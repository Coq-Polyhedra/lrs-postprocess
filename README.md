# lrs-postprocess

`lrs-postprocess` builds certificates from `lrslib` output and exports them
in a binary format readable from Coq using
[`coq-binreader`](https://github.com/Coq-Polyhedra/coq-binreader). The
certificates are checked by the verified checker of `homology-checker`.

The tool uses the output of the lrslib enumeration algorithm. Vertex coordinates and rational
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
--source INDEX Source vertex of the distance certificate, as a 0-based position
               in the lexicographic order of the vertex coordinates (default 0)
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

## Notes

The tool performs no check of its own: the certificate it writes is validated
by the verified checker (extracted to OCaml, or run inside Rocq), which is
where any inconsistency with the input H-representation is detected.

## License

CeCILL-B, see LICENSE.
