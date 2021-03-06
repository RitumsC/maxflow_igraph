# Maxflow for ground-state lattice filling implementation with igraph
This is pure c implementation of maxflow code to calculate ground-state lattice filling in porpus networks. Code mostly was written as a practice.

For igraph library refer to https://igraph.org/.

## Compile with (tested with igraph-0.8.3)
    gcc main.c -I/usr/local/include/igraph -L/usr/local/lib -ligraph -o maxflow

## Usage example

    ./maxflow network_file wettability mu_initial mu_final > output_file

e.g. 

    ./maxflow test_network.dat 2 -20 0 > test.out

The network_file should have following structure
```
Number_of_cells Coordination_number
0 Number_of_neigbouring neighbour_list
1 3 2 3 20
2 4 3 0 6 7
...
21 2 17 18
```
The outputfile will have following structure:
```
mu_i, cells_filled_at_mu_i
mu_2, ...
...
... # will end with mu value before mu_f, but has the same filling as mu_f.
```
for example look refer to test_network.dat and test.out files.
