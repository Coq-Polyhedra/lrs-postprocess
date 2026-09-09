#!/usr/bin/env python3
"""Certificate pipeline and benchmark driver.

Stages, per instance BASE (a BASE.ine file in the data directory):

  lrs     lrsgmp BASE.ine -> BASE.ext
  cert    lrs-postprocess postprocess --bin -> BASE-cert.bin
  check   the extracted checker on BASE-cert.bin
  rocq    the checker run by vm_compute inside Rocq on BASE-cert.bin
  run     lrs, cert and check in sequence
  report  tabulate the measurements of the selected instances (TSV)
  clean   remove every generated file of the instance
  build   cargo build --release

Every stage records what it measured in BASE-<stage>-timings.json: the
wall-clock time and the peak memory of the process, plus the wall-clock
phase timings the tool itself reports (certificate construction, Rust check
and encoding for cert; loading and the three checks for check; loading,
decoding and the three checks for rocq). Rerunning a stage replaces its
record. The tools' raw output is kept in BASE-<stage>.log.
The report is built from the records alone, one row per instance, and never
appended to.

A stage is skipped when its output (for check, its record) exists and is
newer than its inputs and than the tool producing it, unless --force is
given. Outputs are written to a temporary file and moved into place on
success, so an interrupted run leaves nothing a later run could mistake for
a fresh artifact.

Instances are selected by regular expressions matched against the whole
basename of the .ine files, e.g. 'cube(20|21)' or 'dual_cyclic_d1[5-8]_n.*'.
"""

import argparse
import datetime
import json
import os
import re
import resource
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent

# Stderr records of `lrs-postprocess postprocess`, "<label>: <seconds> s".
TIMING_RECORD = re.compile(r"^(?P<label>.*): (?P<secs>[0-9]+(?:\.[0-9]+)?) s$")
# Stderr records of the extracted checker, "<label>  <seconds> s [<verdict>]".
CHECKER_RECORD = re.compile(r"^(?P<label>.*?)\s+(?P<secs>[0-9]+(?:\.[0-9]+)?) s(?:\s+(?P<verdict>\S+))?$")
# Rocq's `Time` output, "Finished transaction in <wall> secs (<user>u,<sys>s) ...".
ROCQ_TIME = re.compile(r"^Finished transaction in (?P<secs>[0-9]+\.?[0-9]*) secs")
# The template's labels, printed before each timed command, and its `Eval` results.
ROCQ_LABEL = re.compile(r"^(?P<label>[a-z ]+):$")
ROCQ_RESULT = re.compile(r"^\s*= (?P<verdict>true|false)$")

# Report columns: (header, stage, phases summed). A single phase is the tool's
# own timing of that phase; None is the wall time of the whole process. The
# cumulative columns T1-T5, T1-T6 and T1-T7 are the conditions checked by the
# three entry points, vertex containment, vertex equality and graph equality,
# taken in sequence.
CHECKS = ["vertex containment", "vertex equality", "graph equality"]
REPORT_COLUMNS = [
    ("cert", "cert", [None]),
    ("cert_build", "cert", ["Create certificate in memory"]),
    ("rust_check", "cert", ["Check certificate"]),
    ("cert_write", "cert", ["Generate and write binary certificate"]),
    ("load", "check", ["certificate loading"]),
    ("T1-T5", "check", CHECKS[:1]),
    ("T1-T6", "check", CHECKS[:2]),
    ("T1-T7", "check", CHECKS[:3]),
    ("rocq_load", "rocq", ["certificate loading"]),
    ("rocq_decode", "rocq", ["certificate decoding"]),
    ("rocq_T1-T5", "rocq", CHECKS[:1]),
    ("rocq_T1-T6", "rocq", CHECKS[:2]),
    ("rocq_T1-T7", "rocq", CHECKS[:3]),
]


class StageFailed(Exception):
    pass


class Config:
    def __init__(self, args):
        self.data = Path(args.data_dir).resolve()
        self.lrsgmp = tool_path(args.lrsgmp)
        self.bin = tool_path(args.bin)
        self.checker = tool_path(args.checker)
        self.rocq_dir = Path(args.rocq_dir).resolve()
        self.coqtop = tool_path(args.coqtop)
        self.rocq_timeout = args.rocq_timeout
        self.force = args.force

    def ine(self, base): return self.data / f"{base}.ine"
    def ext(self, base): return self.data / f"{base}.ext"
    def cert(self, base): return self.data / f"{base}-cert.bin"
    def log(self, base, stage): return self.data / f"{base}-{stage}.log"
    def record(self, base, stage): return self.data / f"{base}-{stage}-timings.json"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def say(msg=""):
    print(msg, flush=True)


def fresh(output, *inputs):
    """True when OUTPUT exists and is newer than every existing input."""
    if not output.exists():
        return False
    mtime = output.stat().st_mtime
    return all(not p.exists() or p.stat().st_mtime <= mtime for p in inputs)


def read_record(cfg, base, stage):
    try:
        return json.loads(cfg.record(base, stage).read_text())
    except (OSError, ValueError):
        return None


def write_record(cfg, base, stage, timed, **fields):
    record = {
        "stage": stage,
        "date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "wall_s": timed.wall,
        "max_rss_mb": timed.max_rss_mb,
        **fields,
    }
    cfg.record(base, stage).write_text(json.dumps(record, indent=2) + "\n")
    return record


def reference(cfg, base):
    """The lrs wall time of BASE, or None."""
    record = read_record(cfg, base, "lrs")
    return record["wall_s"] if record else None


def with_ratio(secs, lrs_s):
    text = f"{secs:.6f} s"
    if lrs_s:
        text += f" ({secs / lrs_s:.3f}x lrs)"
    return text


def summary(stage, timed):
    return f"    {stage} completed in {timed.wall:.3f} s (max rss {timed.max_rss_mb:.0f} MB)"


def tool_path(name):
    """A command given as a path, or found on the PATH."""
    found = shutil.which(name)
    return Path(found).resolve() if found else Path(name).resolve()


def require(path, what):
    if not path.exists():
        raise StageFailed(f"{what} not found: {path}")


def require_tool(path, what):
    if not (path.is_file() and os.access(path, os.X_OK)):
        raise StageFailed(f"{what} not found or not executable: {path}")


class Timed:
    """Wall-clock time and peak memory of one subprocess, from start to wait()."""

    def __enter__(self):
        self.start = time.perf_counter()
        self.wall = self.max_rss_mb = None
        return self

    def wait(self, proc):
        """Wait for PROC, recording its wall time and peak memory; returns its status."""
        _, status, usage = os.wait4(proc.pid, 0)
        self.wall = time.perf_counter() - self.start
        # ru_maxrss is in bytes on macOS and kilobytes on Linux.
        scale = 1 if sys.platform == "darwin" else 1024
        self.max_rss_mb = usage.ru_maxrss * scale / (1024 * 1024)
        proc.returncode = os.waitstatus_to_exitcode(status)
        return proc.returncode

    def __exit__(self, *exc):
        return False


class Staged:
    """A temporary output that is moved into place only on success."""

    def __init__(self, final):
        self.final = final
        self.tmp = final.with_name(final.name + ".tmp")

    def __enter__(self):
        self.tmp.unlink(missing_ok=True)
        return self.tmp

    def __exit__(self, exc_type, *exc):
        if exc_type is None:
            self.tmp.replace(self.final)
        else:
            self.tmp.unlink(missing_ok=True)
        return False


def stream(proc, log, on_line):
    """Copy PROC's stderr to LOG line by line, handing each line to ON_LINE."""
    for line in proc.stderr:
        line = line.rstrip("\n")
        log.write(line + "\n")
        on_line(line)


# ---------------------------------------------------------------------------
# Stages
# ---------------------------------------------------------------------------

def stage_lrs(cfg, base):
    ine, ext = cfg.ine(base), cfg.ext(base)
    require(ine, "input .ine file")
    if not cfg.force and fresh(ext, ine) and read_record(cfg, base, "lrs"):
        say(f"--- lrs: reusing {ext.name} ({reference(cfg, base):.3f} s)")
        return
    require_tool(cfg.lrsgmp, "lrsgmp")
    say(f"--- lrs: {ine.name} -> {ext.name}")
    with Staged(ext) as tmp, open(cfg.log(base, "lrs"), "w") as log, Timed() as t:
        proc = subprocess.Popen([str(cfg.lrsgmp), str(ine), str(tmp)],
                                stdout=log, stderr=subprocess.STDOUT)
        status = t.wait(proc)
        if status != 0:
            raise StageFailed(f"lrsgmp failed with status {status} (see {cfg.log(base, 'lrs')})")
    write_record(cfg, base, "lrs", t)
    say(summary("lrs", t))


def stage_cert(cfg, base):
    ine, ext, cert = cfg.ine(base), cfg.ext(base), cfg.cert(base)
    require(ext, ".ext file (run the lrs stage first)")
    if not cfg.force and fresh(cert, ine, ext, cfg.bin) and read_record(cfg, base, "cert"):
        say(f"--- cert: reusing {cert.name} ({read_record(cfg, base, 'cert')['wall_s']:.3f} s)")
        return
    require_tool(cfg.bin, "lrs-postprocess binary")
    lrs_s = reference(cfg, base)
    phases, accepted = {}, False

    def on_line(line):
        nonlocal accepted
        m = TIMING_RECORD.match(line)
        if m:
            phases[m["label"]] = float(m["secs"])
            say(f"    {m['label'] + ':':<54}{with_ratio(float(m['secs']), lrs_s)}")
        elif line == "Generated certificate accepted":
            accepted = True
        else:
            say(f"    {line}")

    say(f"--- cert: {ext.name} -> {cert.name}")
    with Staged(cert) as tmp, open(tmp, "wb") as out, open(cfg.log(base, "cert"), "w") as log, \
            Timed() as t:
        proc = subprocess.Popen([str(cfg.bin), "postprocess", "--bin", str(ine), str(ext)],
                                stdout=out, stderr=subprocess.PIPE, text=True)
        stream(proc, log, on_line)
        status = t.wait(proc)
        if status != 0:
            raise StageFailed(f"lrs-postprocess failed with status {status} (see {cfg.log(base, 'cert')})")
        if not accepted:
            raise StageFailed("the Rust check did not report acceptance")
    write_record(cfg, base, "cert", t, phases=phases)
    say(summary("cert", t))


def stage_check(cfg, base):
    cert, record = cfg.cert(base), cfg.record(base, "check")
    require(cert, "certificate (run the cert stage first)")
    if not cfg.force and fresh(record, cert, cfg.checker):
        say(f"--- check: reusing {record.name} ({read_record(cfg, base, 'check')['wall_s']:.3f} s)")
        return
    require_tool(cfg.checker, "extracted checker")
    lrs_s = reference(cfg, base)
    phases, verdicts = {}, {}

    def on_line(line):
        m = CHECKER_RECORD.match(line)
        if m:
            phases[m["label"]] = float(m["secs"])
            if m["verdict"]:
                verdicts[m["label"]] = m["verdict"]
            say(f"    {m['label'] + ':':<54}{with_ratio(float(m['secs']), lrs_s)}"
                + (f"  {m['verdict']}" if m["verdict"] else ""))
        else:
            say(f"    {line}")

    say(f"--- check: {cfg.checker.name} {cert.name}")
    with open(cfg.log(base, "check"), "w") as log, Timed() as t:
        proc = subprocess.Popen([str(cfg.checker), str(cert)], stdout=subprocess.DEVNULL,
                                stderr=subprocess.PIPE, text=True)
        stream(proc, log, on_line)
        status = t.wait(proc)
    write_record(cfg, base, "check", t, phases=phases, verdicts=verdicts, status=status,
                 accepted=status == 0)
    total = sum(phases.get(p, 0.0) for p in CHECKS)
    say(f"    checker total {with_ratio(total, lrs_s)}"
        f" [{'accepted' if status == 0 else f'REJECTED(rc={status})'}]")
    say(summary("check", t))
    if status != 0:
        raise StageFailed(f"the extracted checker rejected the certificate (see {cfg.log(base, 'check')})")


def stage_rocq(cfg, base):
    cert, record = cfg.cert(base), cfg.record(base, "rocq")
    require(cert, "certificate (run the cert stage first)")
    template = cfg.rocq_dir / "src" / "CheckCert.v.in"
    require(template, "Rocq check template")
    checker_vos = [cfg.rocq_dir / "src" / f"{m}.vo" for m in ("LowLevelChecker", "CertificateSchema")]
    if not cfg.force and fresh(record, cert, template, *checker_vos):
        say(f"--- rocq: reusing {record.name} ({read_record(cfg, base, 'rocq')['wall_s']:.3f} s)")
        return
    require_tool(cfg.coqtop, "coqtop")
    for vo in checker_vos:
        require(vo, "compiled checker module (run make in the checker directory)")
    lrs_s = reference(cfg, base)
    phases, verdicts, label = {}, {}, None

    def on_line(line):
        nonlocal label
        if m := ROCQ_LABEL.match(line):
            label = m["label"]
        elif (m := ROCQ_TIME.match(line)) and label:
            phases[label] = float(m["secs"])
            say(f"    {label + ':':<54}{with_ratio(phases[label], lrs_s)}"
                + (f"  {verdicts[label]}" if label in verdicts else ""))
        elif (m := ROCQ_RESULT.match(line)) and label:
            verdicts[label] = m["verdict"]
        elif line.startswith(("Error", "Anomaly", "Toplevel input")):
            say(f"    {line}")

    say(f"--- rocq: {cfg.coqtop.name} on {cert.name}")
    script = cfg.data / f"{base}-rocq.v"
    script.write_text(template.read_text().replace("@CERT@", str(cert)))
    with open(cfg.log(base, "rocq"), "w") as log, Timed() as t:
        proc = subprocess.Popen([str(cfg.coqtop), "-batch", "-Q", str(cfg.rocq_dir / "src"), "Cert",
                                 "-l", str(script)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, cwd=cfg.rocq_dir)
        watchdog = threading.Timer(cfg.rocq_timeout, proc.kill)
        watchdog.start()
        proc.stderr = proc.stdout
        stream(proc, log, on_line)
        status = t.wait(proc)
        timed_out = not watchdog.is_alive()
        watchdog.cancel()
    script.unlink(missing_ok=True)
    accepted = status == 0 and all(verdicts.get(c) == "true" for c in CHECKS)
    write_record(cfg, base, "rocq", t, phases=phases, verdicts=verdicts, status=status,
                 accepted=accepted, timed_out=timed_out)
    total = sum(phases.get(c, 0.0) for c in CHECKS)
    say(f"    rocq checks total {with_ratio(total, lrs_s)} [{'accepted' if accepted else 'REJECTED'}]")
    say(summary("rocq", t))
    if timed_out:
        raise StageFailed(f"coqtop killed after {cfg.rocq_timeout} s")
    if not accepted:
        raise StageFailed(f"the Rocq check did not accept the certificate (see {cfg.log(base, 'rocq')})")


def stage_clean(cfg, base):
    removed = []
    paths = [cfg.ext(base), cfg.cert(base)]
    paths += [cfg.log(base, s) for s in ("lrs", "cert", "check", "rocq")]
    paths += [cfg.record(base, s) for s in ("lrs", "cert", "check", "rocq")]
    paths.append(cfg.data / f"{base}-rocq.v")
    for path in paths:
        for p in (path, path.with_name(path.name + ".tmp")):
            if p.exists():
                p.unlink()
                removed.append(p.name)
    say(f"--- clean: removed {', '.join(removed) if removed else 'nothing'}")


def run_instance(cfg, command, base):
    say()
    say(f"########## {base} ##########")
    if command == "clean":
        stage_clean(cfg, base)
        return
    if command in ("lrs", "run"):
        stage_lrs(cfg, base)
    if command in ("cert", "run"):
        stage_cert(cfg, base)
    if command in ("check", "run"):
        stage_check(cfg, base)
    if command == "rocq":
        stage_rocq(cfg, base)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def report(cfg, bases, out, relative):
    """One row per instance. With RELATIVE, every time column is divided by
    the lrs time of the instance (the lrs column stays in seconds)."""
    unit = "" if relative else "_s"
    header = ["instance", "lrs_s"]
    for name, stage, _ in REPORT_COLUMNS:
        header.append(name + unit)
        if name in ("T1-T7", "rocq_T1-T7"):
            header.append("verdict" if stage == "check" else "rocq_verdict")
    out.write("\t".join(header) + "\n")
    for base in bases:
        records = {stage: read_record(cfg, base, stage) for stage in ("lrs", "cert", "check", "rocq")}
        if not any(records.values()):
            continue
        lrs_s = records["lrs"]["wall_s"] if records["lrs"] else None
        row = [base, lrs_s]
        for name, stage, phases in REPORT_COLUMNS:
            record = records[stage]
            if record is None:
                value = None
            elif phases == [None]:
                value = record["wall_s"]
            else:
                value = sum(record["phases"].get(p, 0.0) for p in phases)
            if relative and value is not None:
                value = value / lrs_s if lrs_s else None
            row.append(value)
            if name in ("T1-T7", "rocq_T1-T7"):
                row.append(None if record is None else
                           "accepted" if record["accepted"] else f"REJECTED(rc={record['status']})")
        out.write("\t".join("" if v is None else v if isinstance(v, str) else f"{v:.6f}" for v in row) + "\n")


# ---------------------------------------------------------------------------
# Instance selection and entry point
# ---------------------------------------------------------------------------

def select_instances(data, patterns):
    if not data.is_dir():
        raise StageFailed(f"data directory not found: {data}")
    bases = sorted(p.stem for p in data.glob("*.ine"))
    if not patterns:
        return bases
    selected = []
    for pattern in patterns:
        pattern = pattern.removesuffix(".ine")
        try:
            regex = re.compile(f"^(?:{pattern})$")
        except re.error as e:
            raise StageFailed(f"invalid regular expression {pattern!r}: {e}")
        matches = [b for b in bases if regex.match(b)]
        if not matches:
            raise StageFailed(f"pattern {pattern!r} matches no .ine file in {data}")
        selected.extend(b for b in matches if b not in selected)
    return selected


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=["lrs", "cert", "check", "rocq", "run", "report", "clean", "build"])
    parser.add_argument("patterns", nargs="*", metavar="PATTERN",
                        help="instance patterns (report: all instances when omitted)")
    parser.add_argument("--force", action="store_true", help="recompute stages whose output is fresh")
    parser.add_argument("--data-dir", default="data", metavar="DIR",
                        help="directory of the .ine inputs and generated files (default: data)")
    parser.add_argument("--lrsgmp", default="lrsgmp", metavar="CMD",
                        help="lrs vertex enumerator (default: lrsgmp, from the PATH)")
    parser.add_argument("--bin", default=str(HERE / "target" / "release" / "lrs-postprocess"), metavar="CMD",
                        help="lrs-postprocess binary (default: target/release/lrs-postprocess)")
    parser.add_argument("--checker", default="homology_checker.exe", metavar="CMD",
                        help="extracted checker (default: homology_checker.exe, from the PATH)")
    parser.add_argument("--rocq-dir", default=str(HERE.parent / "homology-checker"), metavar="DIR",
                        help="checker development, holding src/CheckCert.v.in and the compiled modules "
                             "(default: ../homology-checker)")
    parser.add_argument("--coqtop", default="coqtop", metavar="CMD",
                        help="Rocq toplevel (default: coqtop, from the PATH)")
    parser.add_argument("--rocq-timeout", type=float, default=3600, metavar="SECS",
                        help="kill a Rocq check after this long (default: 3600)")
    parser.add_argument("-o", "--output", metavar="FILE",
                        help="report: write the table to FILE instead of standard output")
    parser.add_argument("--relative", action="store_true",
                        help="report: express every time as a multiple of the instance's lrs time")
    args = parser.parse_intermixed_args()

    if args.command == "build":
        if args.patterns:
            parser.error("build takes no instance pattern")
        sys.exit(subprocess.run(["cargo", "build", "--release"], cwd=HERE).returncode)
    if args.command != "report" and not args.patterns:
        parser.error("expected at least one instance pattern")

    cfg = Config(args)
    failures = []
    try:
        bases = select_instances(cfg.data, args.patterns)
        if args.command == "report":
            if args.output:
                with open(args.output, "w") as out:
                    report(cfg, bases, out, args.relative)
            else:
                report(cfg, bases, sys.stdout, args.relative)
            return
        for base in bases:
            try:
                run_instance(cfg, args.command, base)
            except StageFailed as e:
                say(f"    FAILED: {e}")
                failures.append(base)
    except StageFailed as e:
        print(f"error: {e}", file=sys.stderr)
        sys.exit(2)
    except KeyboardInterrupt:
        print("\ninterrupted", file=sys.stderr)
        sys.exit(130)

    say()
    if failures:
        say(f"failed instances: {', '.join(failures)}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
