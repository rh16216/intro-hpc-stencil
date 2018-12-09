skeleton2-headed-plate: skeleton2-heated-plate.c
	mpicc -O2 -std=c99 -Wall $^ -o $@
