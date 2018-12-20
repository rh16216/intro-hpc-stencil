stencil: stencil.c
	mpiicc -O2 -xHOST -std=c99 -Wall $^ -o $@
