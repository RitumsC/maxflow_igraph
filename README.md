# Maxflow for ground-state lattice filling implementation with igraph
This is pure c implementation of maxflow code to calculate ground-state lattice filling in porpus networks. Code mostly was written as a practice.

## Compile with (tested with igraph-0.8.3)
    gcc main.c -I/usr/local/include/igraph -L/usr/local/lib -ligraph -o igraph_test -lm
