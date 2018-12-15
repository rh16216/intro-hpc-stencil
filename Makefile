skeleton2-headed-plate: skeleton2-heated-plate.c
	mpicc -O3 -funsafe-math-optimizations -std=c99 -Wall $^ -o $@
