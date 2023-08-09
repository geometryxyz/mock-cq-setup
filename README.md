# mock-cq-setup

How to run: 

1. cargo run --bin serialize_table 4 LAI serialized/table.bin
table_size = 2**4
seed = LAI
path = serialized/table.bin

2. cargo run --bin run_setup 4 100 serialized/table.bin
table_size = 2**4
toxic_waste = 100
table_path = serialized/table.bin