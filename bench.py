#!/usr/bin/env python3
"""Certificate pipeline and benchmark driver.

Stages, per instance BASE (a BASE.ine file in the data directory):

  lrs     lrsgmp BASE.ine -> BASE.ext
  cert    lrs-postprocess postprocess --bin -> BASE-cert.bin
  run     lrs and cert in sequence
  report  tabulate the measurements of the selected instances (TSV)
  clean   remove every generated file of the instance
  build   cargo build --release

Both stages record what they measured in BASE-<stage>-timings.json: the
wall-clock time and the peak memory of the process, plus, for
cert, the phase timings lrs-postprocess itself reports (certificate
construction, Rust check, encoding). Rerunning a stage replaces its record.
The tools' raw output is kept in BASE-<stage>.log. The report is built from
the records alone, one row per instance, and never appended to.

A stage is skipped when its output exists and is newer than its inputs (and
than the tool producing it), unless --force is given. Outputs are written to
a temporary file and moved into place on success, so an interrupted run
leaves nothing a later run could mistake for a fresh artifact.

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
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent

# Stderr records of `lrs-postprocess postprocess`, "<label>: <seconds> s".
TIMING_RECORD = re.compile(r"^(?P<label>.*): (?P<secs>[0-9]+(?:\.[0-9]+)?) s$")

# Report columns: (header, stage, phase label or None for the process wall time).
REPORT_COLUMNS = [
    ("lrs_s", "lrs", None),
    ("cert_s", "cert", None),
    ("cert_build_s", "cert", "Create certificate in memory"),
    ("rust_check_s", "cert", "Check certificate"),
    ("cert_write_s", "cert", "Generate and write binary certificate"),
]


class StageFailed(Exception):
    pass


class Config:
    def __init__(self, args):
        self.data = Path(args.data_dir).resolve()
        self.lrsgmp = tool_path(args.lrsgmp)
        self.bin = tool_path(args.bin)
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
    """Wall-clock time and peak memory of a subprocess run."""

    def __enter__(self):
        self.start = time.perf_counter()
        return self

    def __exit__(self, *exc):
        self.wall = time.perf_counter() - self.start
        # ru_maxrss is in bytes on macOS and kilobytes on Linux.
        scale = 1 if sys.platform == "darwin" else 1024
        self.max_rss_mb = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * scale / (1024 * 1024)


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
    return proc.wait()


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
        status = subprocess.run([str(cfg.lrsgmp), str(ine), str(tmp)],
                                stdout=log, stderr=subprocess.STDOUT).returncode
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
        status = stream(proc, log, on_line)
        if status != 0:
            raise StageFailed(f"lrs-postprocess failed with status {status} (see {cfg.log(base, 'cert')})")
        if not accepted:
            raise StageFailed("the Rust check did not report acceptance")
    write_record(cfg, base, "cert", t, phases=phases)
    say(summary("cert", t))


def stage_clean(cfg, base):
    removed = []
    paths = [cfg.ext(base), cfg.cert(base)]
    paths += [cfg.log(base, s) for s in ("lrs", "cert")]
    paths += [cfg.record(base, s) for s in ("lrs", "cert")]
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


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def report(cfg, bases, out):
    header = ["instance", *(c[0] for c in REPORT_COLUMNS)]
    out.write("\t".join(header) + "\n")
    for base in bases:
        records = {stage: read_record(cfg, base, stage) for stage in ("lrs", "cert")}
        if not any(records.values()):
            continue

        def cell(stage, phase):
            record = records[stage]
            if record is None:
                return None
            return record["wall_s"] if phase is None else record.get("phases", {}).get(phase)

        cells = [cell(stage, phase) for _, stage, phase in REPORT_COLUMNS]
        out.write("\t".join([base, *("" if v is None else f"{v:.6f}" for v in cells)]) + "\n")


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
    parser.add_argument("command", choices=["lrs", "cert", "run", "report", "clean", "build"])
    parser.add_argument("patterns", nargs="*", metavar="PATTERN",
                        help="instance patterns (report: all instances when omitted)")
    parser.add_argument("--force", action="store_true", help="recompute stages whose output is fresh")
    parser.add_argument("--data-dir", default="data", metavar="DIR",
                        help="directory of the .ine inputs and generated files (default: data)")
    parser.add_argument("--lrsgmp", default="lrsgmp", metavar="CMD",
                        help="lrs vertex enumerator (default: lrsgmp, from the PATH)")
    parser.add_argument("--bin", default=str(HERE / "target" / "release" / "lrs-postprocess"), metavar="CMD",
                        help="lrs-postprocess binary (default: target/release/lrs-postprocess)")
    parser.add_argument("-o", "--output", metavar="FILE",
                        help="report: write the table to FILE instead of standard output")
    args = parser.parse_args()

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
                    report(cfg, bases, out)
            else:
                report(cfg, bases, sys.stdout)
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
