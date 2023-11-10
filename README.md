# Vertex-critical obstructions for C5-coloring (P_t,F)-free graphs

This repository contains code and data related to the exhaustive generation of vertex-critical obstructions for C5-coloring (P_t,F)-free graphs (see the paper "Vertex-critical obstructions for C5-coloring Pt-free graphs"). The two folders "Code" and "Data" are described below.

### CODE

The program can be compiled by using
```bash
make
```

The program will generate all vertex-critical obstructions for C5-coloring except for the triangle K_3 for efficiency reasons. The program expects 3 parameters: the number of vertices in the path (t), the maximum order of a vertex-critical obstruction (n) and the graph F.

Below are some examples of how the program can be used:
```bash
./generateVertexCriticalGraphsC5Coloring -l6 14 < 30VertexGraph.g6
```

This will generate all vertex-critical obstructions for C5-coloring P_6-free graphs on at most 14 vertices. The graph F is chosen as some arbitrary graph on 30 vertices, and since 30>14, all graphs are automatically F-free (but this input cannot be empty).

One of the output lines says:
```bash
There were 0 non-terminating graphs on level 13
```

This indicates that the algorithm was able to terminate and the produced list of vertex-critical obstructions is guaranteed to be exhaustive.

Another usage example is:
```bash
./generateVertexCriticalGraphsC5Coloring -l7 14 < 30VertexGraph.g6
```
This will generate all vertex-critical obstructions for C5-coloring P_7-free graphs on at most 14 vertices. Note that now one of the output lines says:

```bash
There were 2 non-terminating graphs on level 13
```

This indicates that the algorithm cannot guarantee that there are no vertex-critical obstructions with more than 14 vertices. However, the algorithm still guarantees that the produced list of vertex-critical obstructions is exhaustive over all graphs on at most 14 vertices.

If instead we would use:
```bash
./generateVertexCriticalGraphsC5Coloring -l7 17 < 30VertexGraph.g6
```

one of the output lines would be:
```bash
There were 0 non-terminating graphs on level 16
```

This indicates that the algorithm can guarantee that there are no vertex-critical obstructions with more than 17 vertices and the produced list of vertex-critical obstructions is exhaustive.

### DATA

### PtFreeGraphs

The files "vertCritObstructionsC5ColoringPtFreeGraphs.g6" contain an exhaustive list of all vertex-critical obstructions for C5-coloring Pt-free graphs. There is such a file for t=6, t=7 and t=8.

### SubdividedClawFreeGraphs
The files "vertCritObstructionsC5ColoringS_a_b_cFreeGraphs.g6" contain an exhaustive list of all vertex-critical obstructions for C5-coloring a) S_2_2_1-free graphs and b) S_3_1_1-free graphs, where S_a_b_c represents a subdivided claw (e.g. a claw is the graph S_1_1_1).

