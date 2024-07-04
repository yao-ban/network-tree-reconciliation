# Supplementary material for the paper "An efficient algorithm for the reconciliation of a gene network and species tree"


### simulate_species.cpp

Usage:

```
./simulate_species <number of species>
```

This produces a species tree via a pure birth process with birth rate `Srate` (a macro defined in the code, defaulting to 1). The tree is output in the format:

```
<number of nodes> 1
<lines>
```

where each line represents one node in bottom-up order, in the format:

`<node id> <type> <id of parent> <label if a leaf> <time>`

Details:

* node id starts at 0 and increases by 1 up to the number of nodes (no gaps) and is in bottom-up order, so node id of child < node id of parent
* type is 0 = root, 1 = divergence, 2 = reticulation (does not happen for species trees), 3 = leaf
* time starts at 0 for leaf nodes (current day) and increases going backwards in time
* the root node always has no parents and 1 child (not 2 children)


### simulate_gene.cpp

Usage:

```
./simulate_gene <species tree file> [--Drate <D rate>] [--Lrate <L rate>] [--Rrate <R rate>]
```

This produces a gene network via a birth-and-death process inside the species tree given by `<species tree file>`, which must be in the format output by `simulate_species` above.

The birth-and-death process simulates D events at a rate of `<D rate>` per lineage, L events at a rate of `<L rate>` per lineage, and R events at a rate of `<R rate>` per lineage past the first in a species. These rates have default values 0.5.

The network is output in the format:

```
<number of nodes> 0
<lines>
```

where each line represents one node in bottom-up order, in the format:

```
<node id> <type> <id of parent> <id of second parent if a reticulation> <label if a leaf> <breakpoint if a reticulation>
```

Then a reconciliation is output in the format:

```
<lines>
```

where each line represents the mapping of one gene node in id order, in the format:

```
<gene node id> <species node id> <event type>
```

Details:

* event type is C = current, S = speciation, D = duplication, R = reticulation, T = rooT


### perecon.cpp

Usage:

```
./perecon <gene network file> <species tree file>
```

This produces the most parsimonious reconciliation between the gene network given by `<gene network file>` and the species tree given by `<species tree file>`. These must be in the format output by `simulate_gene` and `simulate_species`.

The MPR is calculated with each D event costing `Dcost` (default 1), and each L event costing `Lcost` (default 1), where `Dcost` and `Lcost` are macros defined in the code.

The reconciliation is output in the format:

```
<lines>
```

where each line represents the mapping of one gene node in id order, in the format:

```
<gene node id> <species node id> <event type>
```

Details:

* event type is C = current, S = speciation, D = duplication, R = reticulation, T = rooT
* currently the program also outputs a lot of other information including species tree, gene network, LCA reconciliation, LCA-HCA reconciliation, BCC decomposition, BCC tree-child-ness, intermediate cost calculations, and final MPR cost

### wrapper.py

Usage:

```
python3 wrapper.py [-s <number of species>] [-d <D rate>] [-l <L rate>] [-r <R rate>] [--replicates <number of replicates>]
```

This runs `simulate_species` with `<number of species>` species, then `simulate_gene` with the given rates, then `perecon` to calculate an MPR. It then outputs in one line, separated by `,`:

* number of D events with correct type and location;
* number of D events with correct location only;
* number of D events with correct type only;
* number of D events with neither correct type nor location;
* number of S events with correct type and location;
* number of S events with correct location only;
* number of S events with correct type only;
* number of S events with neither correct type nor location;
* number of R events with correct type and location;
* number of R events with correct location only (always 0);
* number of R events with correct type only;
* number of R events with neither correct type nor location (always 0);
* Proportion of sequence with correct paralogy.

This is then repeated a total of `<number of replicates>` times.

### simulate.py

Usage:

```
python3 simulate.py
```

This runs `wrapper` for the range of parameters shown in the paper.

### process.py

Usage:

```
python3 process.py
```

This opens all `results/results-*.txt` files and amalgamates them into the single file `results/results.csv`. Each line is copied as

```
<line>, <D rate>, <L rate>, <R rate>
```

### analysis.R

This produces (from the `results/results.csv` file) all plots used in the paper.
