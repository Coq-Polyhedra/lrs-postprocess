# lrs-postprocess

Postprocesses `lrslib` output generated with options such as:

```text
allbases
incidence
printcobasis 1
```

The tool trusts `lrs` and does no rational arithmetic. Vertex coordinates are stored as strings.

For H-representation output, a line such as

```text
V#1 R#0 B#1 h=0 facets  12 14 15 16 : 9 10 11 13 I#8 det= 8
 1  0  0  0  1
```

is interpreted as:

- cobasis simplex: `12 14 15 16`, converted to 0-based indices;
- incident inequalities: `12 14 15 16 9 10 11 13`, converted to 0-based indices;
- vertex: the following output vector after the leading homogeneous coordinate `1`.

Starred indices such as `9*` are stripped.

## Usage

```bash
cargo run --release -- path/to/lrs-output.txt --pretty > triples.json
```

Optional mod-2 toggling of duplicate simplices per vertex:

```bash
cargo run --release -- path/to/lrs-output.txt --pretty --mod2 > triples.json
```

## Output

A JSON list of triples:

```json
[
  {
    "vertex": ["0", "0", "0", "1"],
    "incident": [8, 9, 10, 11, 12, 13, 14, 15],
    "simplices": [[3, 5, 6, 7]]
  }
]
```

Each simplex is expressed using local indices into the `incident` list.
