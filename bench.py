#!/usr/bin/env python3
"""Certificate pipeline driver.

Stages, per instance BASE (a BASE.ine file in the data directory):

  lrs    lrsgmp BASE.ine -> BASE.ext              (reference time: BASE-lrs.time)
  cert   lrs-postprocess postprocess --bin        (BASE-cert.bin, timings in BASE-bin.log)
  check  extracted checker on BASE-cert.bin      (BASE-check.log, one row in the results table)
  run    lrs, cert and check in sequence
  clean  remove every generated file of the instance
  build  cargo build --release

A stage is skipped when its output exists and is newer than its inputs (and
than the tool producing it), unless --force is given. Outputs are written to
a temporary file and moved into place on success, so an interrupted run
leaves nothing a later run could mistake for a fresh artifact.

Instances are selected by regular expressions matched against the whole
basename of the .ine files, e.g. 'cube(20|21)' or 'dual_cyclic_d1[5-8]_n.*'.
"""

import argparse
import datetime
import os
import re
import resource
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent

RESULTS_COLUMNS = [
    "instance", "lrs_s", "cert_s", "cert_build_s", "rust_check_s", "cert_write_s",
    "load_s", "vtx_containment_s", "vtx_equality_s", "graph_equality_s",
    "checker_total_s", "checker_vs_lrs", "verdict", "date",
]

# Stderr records of `lrs-postprocess postprocess`, "<label>: <seconds> s".
TIMING_RECORD = re.compile(r"^(?P<label>.*): (?P<secs>[0-9]+(?:\.[0-9]+)?) s$")
# Stderr records of the extracted checker, "<label>  <seconds> s [<verdict>]".
CHECKER_RECORD = re.compile(r"^(?P<label>.*?)\s+(?P<secs>[0-9]+(?:\.[0-9]+)?) s(?:\s+(?P<verdict>\S+))?$")


class StageFailed(Exception):
    pass


class Config:
    def __init__(self, args):
        self.data = Path(args.data_dir).resolve()
        self.lrsgmp = tool_path(args.lrsgmp)
        self.bin = tool_path(args.bin)
        self.checker = tool_path(args.checker)
        self.results = Path(args.results).resolve() if args.results else self.data / "bench-results.tsv"
        self.force = args.force

    def ine(self, base): return self.data / f"{base}.ine"
    def ext(self, base): return self.data / f"{base}.ext"
    def cert(self, base): return self.data / f"{base}-cert.bin"
    def lrs_time(self, base): return self.data / f"{base}-lrs.time"
    def cert_time(self, base): return self.data / f"{base}-cert.time"
    def ext_log(self, base): return self.data / f"{base}-ext.log"
    def bin_log(self, base): return self.data / f"{base}-bin.log"
    def check_log(self, base): return self.data / f"{base}-check.log"


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


def read_seconds(path):
    try:
        return float(path.read_text().strip())
    except (OSError, ValueError):
        return None


def write_seconds(path, secs):
    path.write_text(f"{secs:.6f}\n")


def fmt(secs):
    return "?" if secs is None else f"{secs:.6f}"


def ratio(secs, reference):
    if secs is None or not reference:
        return "n/a"
    return f"{secs / reference:.3f}"


def with_ratio(secs, reference):
    text = f"{secs:.6f} s"
    if reference:
        text += f" ({secs / reference:.3f}x lrs)"
    return text


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
    """Wall-clock and child resource usage of a subprocess run."""

    def __enter__(self):
        self.start = time.perf_counter()
        self.usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        return self

    def __exit__(self, *exc):
        self.wall = time.perf_counter() - self.start
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        self.user = after.ru_utime - self.usage.ru_utime
        self.system = after.ru_stime - self.usage.ru_stime
        # ru_maxrss is in bytes on macOS and kilobytes on Linux.
        scale = 1 if sys.platform == "darwin" else 1024
        self.max_rss_mb = after.ru_maxrss * scale / (1024 * 1024)


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


# ---------------------------------------------------------------------------
# Stages
# ---------------------------------------------------------------------------

def stage_lrs(cfg, base):
    ine, ext = cfg.ine(base), cfg.ext(base)
    require(ine, "input .ine file")
    if not cfg.force and fresh(ext, ine):
        say(f"--- lrs: reusing {ext.name} ({fmt(read_seconds(cfg.lrs_time(base)))} s)")
        return read_seconds(cfg.lrs_time(base))
    require_tool(cfg.lrsgmp, "lrsgmp")
    say(f"--- lrs: {ine.name} -> {ext.name}")
    with Staged(ext) as tmp, open(cfg.ext_log(base), "w") as log:
        with Timed() as t:
            status = subprocess.run([str(cfg.lrsgmp), str(ine), str(tmp)],
                                    stdout=log, stderr=subprocess.STDOUT).returncode
        if status != 0:
            raise StageFailed(f"lrsgmp failed with status {status} (see {cfg.ext_log(base)})")
    write_seconds(cfg.lrs_time(base), t.wall)
    say(f"    lrs completed in {t.wall:.3f} s (user {t.user:.2f} s, max rss {t.max_rss_mb:.0f} MB)")
    return t.wall


def stage_cert(cfg, base, lrs_s):
    ine, ext, cert = cfg.ine(base), cfg.ext(base), cfg.cert(base)
    require(ext, ".ext file (run the lrs stage first)")
    if not cfg.force and fresh(cert, ine, ext, cfg.bin):
        say(f"--- cert: reusing {cert.name} ({fmt(read_seconds(cfg.cert_time(base)))} s)")
        return read_seconds(cfg.cert_time(base))
    require_tool(cfg.bin, "lrs-postprocess binary")
    say(f"--- cert: {ext.name} -> {cert.name}")
    with Staged(cert) as tmp, open(tmp, "wb") as out, open(cfg.bin_log(base), "w") as log:
        with Timed() as t:
            proc = subprocess.Popen([str(cfg.bin), "postprocess", "--bin", str(ine), str(ext)],
                                    stdout=out, stderr=subprocess.PIPE, text=True)
            accepted = False
            for line in proc.stderr:
                line = line.rstrip("\n")
                log.write(line + "\n")
                m = TIMING_RECORD.match(line)
                if m:
                    say(f"    {m['label'] + ':':<54}{with_ratio(float(m['secs']), lrs_s)}")
                elif line == "Generated certificate accepted":
                    accepted = True
                else:
                    say(f"    {line}")
            status = proc.wait()
        if status != 0:
            raise StageFailed(f"lrs-postprocess failed with status {status} (see {cfg.bin_log(base)})")
        if not accepted:
            raise StageFailed("the Rust check did not report acceptance")
    write_seconds(cfg.cert_time(base), t.wall)
    say(f"    cert completed in {t.wall:.3f} s (user {t.user:.2f} s, max rss {t.max_rss_mb:.0f} MB)")
    return t.wall


def bin_log_totals(cfg, base):
    """The three phase totals recorded by lrs-postprocess, from BASE-bin.log."""
    wanted = {
        "Create certificate in memory": None,
        "Check certificate": None,
        "Generate and write binary certificate": None,
    }
    try:
        for line in cfg.bin_log(base).read_text().splitlines():
            m = TIMING_RECORD.match(line)
            if m and m["label"] in wanted:
                wanted[m["label"]] = float(m["secs"])
    except OSError:
        pass
    return tuple(wanted.values())


def stage_check(cfg, base, lrs_s, cert_s):
    cert = cfg.cert(base)
    require(cert, "certificate (run the cert stage first)")
    require_tool(cfg.checker, "extracted checker")
    say(f"--- check: {cfg.checker.name} {cert.name}")
    with open(cfg.check_log(base), "w") as log, Timed() as t:
        proc = subprocess.Popen([str(cfg.checker), str(cert)], stdout=subprocess.DEVNULL,
                                stderr=subprocess.PIPE, text=True)
        records = {}
        for line in proc.stderr:
            line = line.rstrip("\n")
            log.write(line + "\n")
            m = CHECKER_RECORD.match(line)
            if m:
                records[m["label"]] = (float(m["secs"]), m["verdict"])
                say(f"    {m['label'] + ':':<54}{with_ratio(float(m['secs']), lrs_s)}"
                    + (f"  {m['verdict']}" if m["verdict"] else ""))
            else:
                say(f"    {line}")
        status = proc.wait()

    def secs(label):
        return records[label][0] if label in records else None

    parts = [secs("certificate loading"), secs("vertex containment"),
             secs("vertex equality"), secs("graph equality")]
    total = sum(p for p in parts if p is not None)
    verdict = "accepted" if status == 0 else f"REJECTED(rc={status})"
    build_s, rust_check_s, write_s = bin_log_totals(cfg, base)
    row = [base, fmt(lrs_s), fmt(cert_s), fmt(build_s), fmt(rust_check_s), fmt(write_s),
           *(("?" if p is None else f"{p:.6f}") for p in parts),
           f"{total:.6f}", ratio(total, lrs_s), verdict,
           datetime.datetime.now().strftime("%Y-%m-%d %H:%M")]
    if not cfg.results.exists():
        cfg.results.write_text("\t".join(RESULTS_COLUMNS) + "\n")
    with open(cfg.results, "a") as f:
        f.write("\t".join(row) + "\n")
    say(f"    checker total {total:.6f} s ({ratio(total, lrs_s)}x lrs) [{verdict}]"
        f" (wall {t.wall:.3f} s, max rss {t.max_rss_mb:.0f} MB)")
    if status != 0:
        raise StageFailed(f"the extracted checker rejected the certificate (see {cfg.check_log(base)})")


def stage_clean(cfg, base):
    removed = []
    for path in [cfg.ext(base), cfg.cert(base), cfg.lrs_time(base), cfg.cert_time(base),
                 cfg.ext_log(base), cfg.bin_log(base), cfg.check_log(base)]:
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
    lrs_s = read_seconds(cfg.lrs_time(base))
    cert_s = read_seconds(cfg.cert_time(base))
    if command in ("lrs", "run"):
        lrs_s = stage_lrs(cfg, base)
    if command in ("cert", "run"):
        cert_s = stage_cert(cfg, base, lrs_s)
    if command in ("check", "run"):
        stage_check(cfg, base, lrs_s, cert_s)


# ---------------------------------------------------------------------------
# Instance selection and entry point
# ---------------------------------------------------------------------------

def select_instances(data, patterns):
    if not data.is_dir():
        raise StageFailed(f"data directory not found: {data}")
    bases = sorted(p.stem for p in data.glob("*.ine"))
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
    parser.add_argument("command", choices=["lrs", "cert", "check", "run", "clean", "build"])
    parser.add_argument("patterns", nargs="*", metavar="PATTERN")
    parser.add_argument("--force", action="store_true", help="recompute stages whose output is fresh")
    parser.add_argument("--data-dir", default="data", metavar="DIR",
                        help="directory of the .ine inputs and generated files (default: data)")
    parser.add_argument("--lrsgmp", default="lrsgmp", metavar="CMD",
                        help="lrs vertex enumerator (default: lrsgmp, from the PATH)")
    parser.add_argument("--bin", default=str(HERE / "target" / "release" / "lrs-postprocess"), metavar="CMD",
                        help="lrs-postprocess binary (default: target/release/lrs-postprocess)")
    parser.add_argument("--checker", default="homology_checker.exe", metavar="CMD",
                        help="extracted checker (default: homology_checker.exe, from the PATH)")
    parser.add_argument("--results", metavar="FILE",
                        help="results table appended by the check stage (default: DATA_DIR/bench-results.tsv)")
    args = parser.parse_args()

    if args.command == "build":
        if args.patterns:
            parser.error("build takes no instance pattern")
        sys.exit(subprocess.run(["cargo", "build", "--release"], cwd=HERE).returncode)
    if not args.patterns:
        parser.error("expected at least one instance pattern")

    cfg = Config(args)
    failures = []
    try:
        bases = select_instances(cfg.data, args.patterns)
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
    if args.command in ("check", "run"):
        say(f"results table: {cfg.results}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
