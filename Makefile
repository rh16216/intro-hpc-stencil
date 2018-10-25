stencil: stencil.c
	icc -O2 -xHOST -std=c99 -Wall $^ -o $@
