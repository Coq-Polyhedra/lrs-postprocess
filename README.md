# lrs-postprocess

`lrs-postprocess` builds and checks certificates from `lrslib` output. It can
also export certificates in a binary format readable from Coq using
[`coq-binreader`](https://github.com/Coq-Polyhedra/coq-binreader).

The tool trusts the `lrs` enumeration output. Vertex coordinates and rational
numbers are kept as strings in the JSON certificate.

## Build

```bash
make
```

`make test` runs the unit tests and `make clean` removes the build directory.

The binary is then:

```bash
./target/release/lrs-postprocess
```

## Command-line interface

### Generate a JSON certificate

```bash
./target/release/lrs-postprocess postprocess input.ine output.ext --pretty > certificate.json
```

Options:

```text
--pretty       Pretty-print JSON output
--k0 INDEX     Choose the root vertex item index, 0-based
--bin          Output the certificate in binary format for coq-binreader
```

Examples:

```bash
./target/release/lrs-postprocess postprocess data/cross3.ine data/cross3.ext --pretty > data/cross3-cert.json
./target/release/lrs-postprocess postprocess data/cross3.ine data/cross3.ext --pretty --k0 0 > data/cross3-cert.json
```

### Generate a binary certificate

```bash
./target/release/lrs-postprocess postprocess --bin data/cross3.ine data/cross3.ext > data/cross3-cert.bin
```

The binary certificate is intended to be loaded from Coq with `coq-binreader`.

## JSON certificate format

A certificate has the following shape:

```json
{
  "items": [
    {
      "vertex": ["0", "0", "0", "1"],
      "incident": [8, 9, 10, 11, 12, 13, 14, 15],
      "simplices": [[3, 5, 6, 7]]
    }
  ],
  "root": {
    "k0": 0,
    "rows": [8, 9, 10, 11],
    "inverse_rows": [["1", "0", "0", "0"]],
    "same_label_separators": []
  }
}
```

Each simplex is expressed using local indices into the corresponding
`incident` list.

## Binary certificate format

The binary certificate stores the same data as the JSON certificate, encoded
for `coq-binreader`.

The schema is:

```text
certificate := items * roots

items := array item

item :=
  vertex * (incident * simplices)

vertex := array BigQ
incident := array Int63
simplices := array (array Int63)

roots := array root
```

The `roots` array has length `0` if there is no root and length `1` otherwise.

```text
root :=
  k0 * (rows * (inverse_rows * same_label_separators))

k0 := Int63
rows := array Int63
inverse_rows := array (array BigQ)
same_label_separators := array Int63
```

## Pipeline script

`bench.py` drives the whole pipeline, one stage per command:

```text
lrs     lrsgmp BASE.ine -> BASE.ext
cert    lrs-postprocess postprocess --bin -> BASE-cert.bin
check   the extracted checker on BASE-cert.bin
rocq    the checker run by vm_compute inside Rocq on BASE-cert.bin
run     lrs, cert and check in sequence
report  tabulate the measurements of the selected instances (TSV)
clean   remove the generated files of an instance
build   cargo build --release
```

Instances are selected by regular expressions matched against the complete
basename of the `.ine` files in the data directory, so `cross3` selects
exactly `cross3.ine` and `'cube(20|21)'` selects both cubes. A stage is
skipped when its output is newer than its inputs and than the tool producing
it; `--force` recomputes it. Outputs are staged in a `.tmp` file and moved
into place on success.

Every stage writes what it measured to `BASE-<stage>-timings.json`: the
wall-clock time and the peak memory of the process, measured in Python, plus
the wall-clock phase timings reported by the tool itself (certificate
construction, Rust check and encoding for `cert`; loading and the three checks
for `check`; loading, decoding and the three checks for `rocq`, from Rocq's
`Time`). The tool's raw output is kept in `BASE-<stage>.log`. `report` builds
a tab-separated table from these records, one row per instance with blank
cells for stages that have not run, and writes it to standard output or to `-o
FILE`. Its columns are the lrs and certificate-generation wall times, the
three generation phases, and for each of the two checkers its loading time
(plus decoding for Rocq), the cumulative check times T1-T5 (vertex
containment), T1-T6 (plus vertex equality) and T1-T7 (plus graph equality),
and the verdict. With `--relative` every time but the lrs one is divided by
the instance's lrs time. The `rocq` stage instantiates `src/CheckCert.v.in`
from the checker development with the certificate path and runs it through
`coqtop -batch`; each Rocq check decodes the certificate itself, so the
decoding time is included in each of the three check times and is also
reported on its own.

```bash
./bench.py build
./bench.py run cross3
./bench.py --force lrs 'cross(8|9|10)'
./bench.py cert 'dual_cyclic_d1[5-8]_n.*'
./bench.py check cube20
./bench.py rocq cube15
./bench.py report -o results.tsv
./bench.py report --relative 'cube.*'
./bench.py clean cross3
```

Options:

```text
--data-dir DIR   directory of the .ine inputs and generated files (default: data)
--lrsgmp CMD     lrs vertex enumerator (default: lrsgmp, from the PATH)
--bin CMD        lrs-postprocess binary (default: target/release/lrs-postprocess)
--checker CMD    extracted checker (default: homology_checker.exe, from the PATH)
--rocq-dir DIR   checker development with src/CheckCert.v.in and its compiled
                 modules (default: ../homology-checker)
--coqtop CMD     Rocq toplevel (default: coqtop, from the PATH)
--rocq-timeout S kill a Rocq check after S seconds (default: 3600)
-o FILE          report: write the table to FILE instead of standard output
--force          recompute stages whose output is fresh
```

## Notes

The checker verifies the consistency of the certificate with the input
H-representation and the generated certificate data. The binary export is
intended as a compact format for importing the same certificate data into Coq.
