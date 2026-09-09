# lrs-postprocess

`lrs-postprocess` builds and checks certificates from `lrslib` output. It can also export certificates in a binary format readable from Coq using [`coq-binreader`](https://github.com/Coq-Polyhedra/coq-binreader).

The tool trusts the `lrs` enumeration output. Vertex coordinates and rational numbers are kept as strings in the JSON certificate.

## Build

```bash
cargo build --release
```

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

Each simplex is expressed using local indices into the corresponding `incident` list.

## Binary certificate format

The binary certificate stores the same data as the JSON certificate, encoded for `coq-binreader`.

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

The repository contains a helper script:

```bash
run-cert-pipeline.sh
```

It assumes that input and generated files are stored in `data/`.

For a base name `BASE`, the script uses:

```text
data/BASE.ine
data/BASE.ext
data/BASE-cert.json
data/BASE-cert.bin
```

Generated logs are:

```text
data/BASE-ext.log
data/BASE-cert.log
data/BASE-bin.log
```

### Build the Rust binary

```bash
./run-cert-pipeline.sh build
```

### Generate the `.ext` file

```bash
./run-cert-pipeline.sh ext cross3
```

### Generate the JSON certificate

```bash
./run-cert-pipeline.sh cert cross3
```

### Generate the binary certificate

```bash
./run-cert-pipeline.sh bin cross3
```

### Run the whole pipeline

```bash
./run-cert-pipeline.sh all cross3
```

On several files:

```bash
./run-cert-pipeline.sh all cross_9 cross_10 cross_11 cross_12
```

### Clean generated files

```bash
./run-cert-pipeline.sh clean cross3
```

## Pipeline environment variables

The script can be customized with environment variables:

```text
DATA_DIR       directory containing input/output files, default: data
LRS_DIR        directory containing lrsgmp
LRSGMP         full path to lrsgmp
BIN            compiled Rust binary, default: ./target/release/lrs-postprocess
CERT_GEN       JSON certificate command, default: "$BIN postprocess --pretty"
BIN_CERT_GEN   binary certificate command, default: "$BIN postprocess --bin"
```

Example:

```bash
DATA_DIR=data \
LRSGMP=/home/user/lrslib/lrsgmp \
./run-cert-pipeline.sh all cross3
```

## Notes

The checker verifies the consistency of the certificate with the input H-representation and the generated certificate data. The binary export is intended as a compact format for importing the same certificate data into Coq.
