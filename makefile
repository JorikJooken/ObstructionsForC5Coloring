compiler=gcc
flags=-std=gnu11 -march=native
densenauty32=-DWORDSIZE=32 -DMAXN=WORDSIZE nautyW1.a
densenauty64=-DWORDSIZE=64 -DMAXN=WORDSIZE nautyL1.a

# The 32-bit version of this program takes less memory but only supports graphs up to 32 vertices.
32bit: generateVertexCriticalGraphsC5Coloring.c bitset.h read_graph/readGraph6.c
	$(compiler) -o generateVertexCriticalGraphsC5Coloring generateVertexCriticalGraphsC5Coloring.c planarity.c read_graph/readGraph6.c $(densenauty32) $(flags) -O4

# The 64-bit version of this program.
64bit: generateVertexCriticalGraphsC5Coloring.c bitset.h read_graph/readGraph6.c
	$(compiler) -o generateVertexCriticalGraphsC5Coloring-64 generateVertexCriticalGraphsC5Coloring.c planarity.c read_graph/readGraph6.c $(densenauty64) $(flags) -O4

all: 32bit 64bit

.PHONY: clean
clean:
	rm -f generateVertexCriticalGraphsC5Coloring generateVertexCriticalGraphsC5Coloring-64
