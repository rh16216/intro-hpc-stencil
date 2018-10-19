stencil: stencil.c
	gcc -O3 -funsafe-math-optimizations -std=c99 -Wall $^ -o $@

