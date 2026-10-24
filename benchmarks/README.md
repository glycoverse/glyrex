# In-memory performance and parity checks

Run from the package root:

```sh
Rscript benchmarks/functions.R
Rscript benchmarks/parity.R
```

Both scripts load the original baseline (`2a808806`) and the working-tree
implementation into separate, byte-compiled R environments. Dependencies and
input structure parsing are identical. They print to the console and do not
write timings, datasets, reports, or serialized results to disk.

The performance check uses the 26 distinct fixture glycans and a precompiled
`(?<terminal>Hex)-HexNAc` pattern. Every input is unique. It includes warmup,
calibration, nine measured batches, and GC during timing; explicit GC before
each batch is outside timing. Loading and input parsing are outside timing.
The baseline revision can be supplied as the first command-line argument.

The parity check compares function results, missing values, names, captures,
structure representations and raw graph node IDs. Structure objects are
compared through their typed, named IUPAC values rather than igraph pointers.
