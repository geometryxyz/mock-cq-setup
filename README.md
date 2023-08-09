# mock-cq-setup

How to run:

1. Generate a table and serialize it to `table.bin`. The expected arguments are
`{k} {random_seed} {table_path}`, where the number of table entries is `1<<k`.
`{random_seed}` is used to generate the table entries; if it isn't provided, the
table will simply be the range `[1..1<<k]`.
```console
cargo run --bin serialize_table 4 LAI serialized/table.bin
```

2. Generate cq mock-srs from `table.bin` and measure time needed to generate it.
The expected arguments are `{k} {toxic_waste} {table_path}`, where the setup size
is `2**k`. **`k` must be the same in steps 1 and 2.**
```console
cargo run --bin run_setup 4 100 serialized/table.bin
```