# k-vertex-critical 
This repository contains code for exhaustively generating vertex-critical obstructions for C5-coloring (P_t,F)-free graphs. In fact, it will generate all such obstructions except for the triangle K_3 for efficiency reasons.

The program can be compiled by using
```bash
make
```

The program expects 3 parameters: the number of vertices in the path (t), the maximum order of a vertex-critical obstruction (n) and the graph F.

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
