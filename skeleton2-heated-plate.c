/*
** A simple example of halo-exchange for a 2d grid.
** The example is heat diffusion on a heated plate
**
** Boundary conditions for the full grid are:
**
**                      W = 0
**             +--------------------+
**             |                    |
**    W = 100  |                    | W = 100
**             |                    |
**             +--------------------+
**                     W = 100
**
** i.e. 3 sides are held at 100 degress, while the fourth
** is held at 0 degrees.
**
** The grid will be partitioned into 4 subgrids, used by
** each of four ranks:
**
**                       W = 0
**                   |     |     |
**             +-----|-----|-----|-----+
**             |     |     |     |     |
**    W = 100  |     |     |     |     | W = 100
**             |     |     |     |     |
**             +-----|-----|-----|-----+
**                   |     |     |
**                      W = 100
**
** A pattern of communication using only column-based
** halos will be employed, e.g. for 4 ranks:
**
**   +-----+     +-----+     +-----+     +-----+
**   ||   ||     ||   ||     ||   ||     ||   ||
** <-|| 0 || <-> || 1 || <-> || 2 || <-> || 3 || ->
**   ||   ||     ||   ||     ||   ||     ||   ||
**   +-----+     +-----+     +-----+     +-----+
**
**
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mpi.h"


#define EPSILON 0.01
#define MASTER 0

// Define output file name
#define OUTPUT_FILE "stencil.pgm"

/* function prototypes */
int calc_nrows_from_rank(int rank, int size, int nx);
double wtime(void);

// Create the input image
void init_image(const int nx, const int ny, float **  image, float **  tmp_image, int size, int rank, int local_nrows) {
  //Generate entire Checkerboard
  float **checker = (float**)malloc(sizeof(float*) * (nx+2));
  for (int ii=0;ii<nx+2;ii++) {
    checker[ii] = (float*)malloc(sizeof(float) * (ny+2));
  }

  // Zero everything
    for (int i = 0; i < nx+2; ++i) {
      for (int j = 0; j < ny+2; ++j) {
      checker[i][j] = 0.0;
    }
  }

    for (int i = 0; i < 8; ++i) {
      for (int j = 0; j < 8; ++j) {
        for (int ii = (i*nx/8)+1; ii < ((i+1)*nx/8)+1; ++ii) {
          for (int jj = (j*ny/8)+1; jj < ((j+1)*ny/8)+1; ++jj) {
          if ((i+j)%2){
            checker[ii][jj] = 100.0;
          }
        }
      }
    }
  }

   //Copy specific section for current rank
   //printf("Start: %d\n", rank*(ny/size)+1);
   //printf("End: %d\n", rank*(ny/size)+local_ncols);
   //printf("Local n cols: %d\n", local_ncols);
     for (int i = rank*(ny/size)+1; i <= rank*(nx/size)+local_nrows; ++i) {
       for (int j = 0; j < ny+2; ++j) {
       image[i-rank*(ny/size)][j] = checker[i][j];
       tmp_image[i-rank*(ny/size)][j] = checker[i][j];
   }
 }

}

// Routine to output the image in Netpbm grayscale binary image format
void output_image(const char * file_name, const int nx, const int ny, float **image) {

  // Open output file
  FILE *fp = fopen(file_name, "w");
  if (!fp) {
    fprintf(stderr, "Error: Could not open %s\n", OUTPUT_FILE);
    exit(EXIT_FAILURE);
  }

  // Ouptut image header
  fprintf(fp, "P5 %d %d 255\n", nx, ny);

  // Calculate maximum value of image
  // This is used to rescale the values
  // to a range of 0-255 for output
  double maximum = 0.0;
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      if (image[i][j] > maximum)
        maximum = image[i][j];
    }
  }

  // Output image, converting to numbers 0-255
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      fputc((char)(255.0*image[i][j]/maximum), fp);
    }
  }

  // Close the file
  fclose(fp);

}


void stencil(float **u, float **w, float *sendbuf, float *recvbuf,
  int local_nrows, int local_ncols, int rank, int size, MPI_Status status){

  int ii, jj;             /* row and column indices for the grid */
  int tag = 0;            /* scope for adding extra information to a message */
  int up;                 /* the rank of the process above */
  int down;               /* the rank of the process below */

  /*
  ** determine process ranks above and below the rank
  ** respecting periodic boundary conditions
  */
  up = (rank == MASTER) ? (rank + size - 1) : (rank - 1);
  down = (rank + 1) % size;
  /*
  ** halo exchange for the local grids w:
  ** - first send to the top and receive from the bottom,
  ** - then send to the bottom and receive from the top.
  ** for each direction:
  ** - first, pack the send buffer using values from the grid
  ** - exchange using MPI_Sendrecv()
  ** - unpack values from the recieve buffer into the grid
  */

  /* send to the top, receive from bottom */
  if (rank != 0){
    for(ii=0;ii<local_ncols+2;ii++){
      sendbuf[ii] = w[1][ii];
      //printf("up send %6.2f\n", sendbuf[ii]);
    }
    if (rank == size-1){
      MPI_Send(sendbuf, local_ncols+2, MPI_FLOAT,
               up, tag, MPI_COMM_WORLD);
    }
    else{
      //printf("Start send up from rank: %d\n", rank);
      MPI_Sendrecv(sendbuf, local_ncols+2, MPI_FLOAT, up, tag,
       recvbuf, local_ncols+2, MPI_FLOAT, down, tag,
       MPI_COMM_WORLD, &status);
      //printf("End send up from rank: %d\n", rank);
    }
   }
    if (rank != size-1){
      if (rank == 0){
        MPI_Recv(recvbuf, local_ncols+2, MPI_FLOAT,
                 down, tag, MPI_COMM_WORLD, &status);
      }
      //printf("Start recv up to rank: %d\n", rank);
      for(ii=0;ii<local_ncols+2;ii++){
        w[local_nrows + 1][ii] = recvbuf[ii];
        //printf("up recv %6.2f\n", w[ii][local_ncols + 1]);
      }
      //printf("End recv up to rank: %d\n", rank);
    }
  /* send to the bottom, receive from top */
  if (rank != size-1){
    for(ii=0;ii<local_ncols+2;ii++){
      sendbuf[ii] = w[local_nrows][ii];
      //printf("down send %6.2f\n", sendbuf[ii]);
    }
    if (rank == 0){
      MPI_Send(sendbuf, local_ncols+2, MPI_FLOAT,
               down, tag, MPI_COMM_WORLD);
    }
    else{
      //printf("Start send down from rank: %d\n", rank);
      MPI_Sendrecv(sendbuf, local_ncols+2, MPI_FLOAT, down, tag,
       recvbuf, local_ncols+2, MPI_FLOAT, up, tag,
       MPI_COMM_WORLD, &status);
      //printf("End send down from rank: %d\n", rank);
    }
   }
   if (rank != 0){
    if (rank == size-1){
       MPI_Recv(recvbuf, local_ncols+2, MPI_FLOAT,
                up, tag, MPI_COMM_WORLD, &status);
    }
    //printf("Start recv down to rank: %d\n", rank);
    for(ii=0;ii<local_ncols+2;ii++){
      w[0][ii] = recvbuf[ii];
      //printf("down recv %6.2f\n", w[ii][0]);
    }
    //printf("End recv down to rank: %d\n", rank);
  }
  /*
  ** copy the old solution into the u grid
  */
  //for (int q = 0; q < local_nrows+2; q++){
  //  printf("%f\n", w[q][0]);
  //}
  //printf("NO MORE COMMUNICATION\n");
  /*
  ** compute new values of w using u
  ** looping extents depend on rank, as we don't
  ** want to overwrite any boundary conditions
  */
  for(ii=1;ii<local_nrows+1;ii++) {
//     if(rank == 0) {
// start_col = 2;
// end_col = local_ncols;
//     }
//     else if(rank == size -1) {
// start_col = 1;
// end_col = local_ncols - 1;
//     }
//     else {
// start_col = 1;
// end_col = local_ncols;
//     }
    for(jj=1;jj<local_ncols+1;jj++) {
      u[ii][jj] = w[ii][jj-1] * 0.1f;
      u[ii][jj] += w[ii][jj] * 0.6f;
      u[ii][jj] += w[ii][jj+1] * 0.1f;
      u[ii][jj] += w[ii-1][jj] * 0.1f;
      u[ii][jj] += w[ii+1][jj] * 0.1f;
    }
  }
}

int main(int argc, char* argv[])
{
  int ii, jj;             /* row and column indices for the grid */
  int kk;                /* index for looping over ranks */
  //int start_col,end_col; /* rank dependent looping indices */
  int iter;              /* index for timestep iterations */
  int rank;              /* the rank of this process */
  int size;              /* number of processes in the communicator */
  MPI_Status status;     /* struct used by MPI_Recv */
  int tag = 0;           /* scope for adding extra information to a message */
  int local_nrows;       /* number of rows apportioned to this rank */
  int local_ncols;       /* number of columns apportioned to this rank */
  int remote_nrows;      /* number of columns apportioned to a remote rank */
  float **u;            /* local temperature grid at time t - 1 */
  float **w;            /* local temperature grid at time t     */
  float **out;          /* grid for final result      */
  float *sendbuf;       /* buffer to hold values to send */
  float *recvbuf;       /* buffer to hold received values */
  float *printbuf;      /* buffer to hold values for printing */

  // Initiliase problem dimensions from command line arguments
  int nx = atoi(argv[1]);
  int ny = atoi(argv[2]);
  int niters = atoi(argv[3]);

  /* MPI_Init returns once it has started up processes */
  /* get size and rank */
  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );


  /*
  ** determine local grid size
  ** each rank gets all the rows, but a subset of the number of columns
  */
  local_nrows = calc_nrows_from_rank(rank, size, nx);
  local_ncols = ny;

  /*
  ** allocate space for:
  ** - the local grid (2 extra columns added for the halos)
  ** - we'll use local grids for current and previous timesteps
  ** - buffers for message passing
  */
   u = (float**)malloc(sizeof(float*) * (local_nrows + 2));
   for(ii=0;ii<local_nrows+2;ii++) {
     u[ii] = (float*)malloc(sizeof(float) * (local_ncols + 2));
   }
   w = (float**)malloc(sizeof(float*) * (local_nrows + 2));
   for(ii=0;ii<local_nrows+2;ii++) {
     w[ii] = (float*)malloc(sizeof(float) * (local_ncols + 2));
   }

   out = (float**)malloc(sizeof(float*) * nx);
   for(ii=0;ii<nx;ii++) {
     out[ii] = (float*)malloc(sizeof(float) * ny);
   }

   sendbuf = (float*)malloc(sizeof(float) * (local_ncols+2));
   recvbuf = (float*)malloc(sizeof(float) * (local_ncols+2));
  //    printbuf must be big enough to hold this number */
  // /* The last rank has the most rows apportioned.
   remote_nrows = calc_nrows_from_rank(size-1, size, nx);
   printbuf = (float*)malloc(sizeof(float) * (local_ncols + 2));

  /*
  ** initialize the local grid for the present time (w):
  ** - set boundary conditions for any boundaries that occur in the local grid
  ** - initialize inner cells to the average of all boundary cells
  ** note the looping bounds for index jj is modified
  ** to accomodate the extra halo columns
  ** no need to initialise the halo cells at this point
  */

  init_image(nx, ny, u, w, size, rank, local_nrows);

  // for (int i = 0; i < local_nrows+2; i++){
  //   for (int j = 0; j < local_ncols+2; j++){
  //     printf("%6.2f ", w[i][j]);
  //   }
  //     printf("\n");
  // }


  /*
  ** time loop
  */
  double tic = wtime();
  for(iter=0;iter<niters;iter++) {
    stencil(u, w, sendbuf, recvbuf, local_nrows, local_ncols, rank, size, status);
    stencil(w, u, sendbuf, recvbuf, local_nrows, local_ncols, rank, size, status);
   }
   double toc = wtime();
   if (rank == 0){
     // Output
     printf("------------------------------------\n");
     printf(" runtime: %lf s\n", toc-tic);
     printf("------------------------------------\n");
   }
    // if (rank == 0){
    //   for (int i = 0; i < local_nrows+2; i++){
    //     for (int j = 0; j < local_ncols+2; j++){
    //       printf("%6.2f ", w[i][j]);
    //     }
    //     printf("\n");
    //   }
    // }

  /*
  ** at the end, write out the solution.
  ** for each row of the grid:
  ** - rank 0 first prints out its cell values
  ** - then it receives row values sent from the other
  **   ranks in order, and prints them.
  */
  if(rank == MASTER) {
    printf("nx: %d\nny: %d\n",nx,ny);
    printf("Final result after gaussian blur:\n");
  }

  for(ii=1;ii<local_nrows+1;ii++) {
    if(rank == 0) {
      for(jj=1;jj<local_ncols + 1;jj++) {
	       out[ii-1][jj-1] = w[ii][jj];
      }
      for(kk=1;kk<size;kk++) { /* loop over other ranks */
	       MPI_Recv(printbuf,local_ncols + 2,MPI_FLOAT,kk,tag,MPI_COMM_WORLD,&status);
	         for(jj=1;jj<local_ncols + 1;jj++) {
	            out[kk*local_nrows+ii-1][jj-1] = printbuf[jj];
	         }
      }
    }
    else {
      MPI_Send(w[ii],local_ncols + 2,MPI_FLOAT,MASTER,tag,MPI_COMM_WORLD);
    }
  }

  for(ii=local_nrows+1; ii<remote_nrows+1;ii++) {
    if(rank == 0) {
	     MPI_Recv(printbuf,local_ncols + 2,MPI_FLOAT,size-1,tag,MPI_COMM_WORLD,&status);
	     for(jj=1;jj<local_ncols + 1;jj++) {
	        out[(size-1)*local_nrows+ii-1][jj-1] = printbuf[jj];
	     }
    }
    else if(rank == size-1) {
      MPI_Send(w[ii],local_ncols + 2,MPI_FLOAT,MASTER,tag,MPI_COMM_WORLD);
    }
  }


  //if(rank == MASTER){
  //   for (int i = 0; i < nx; i++){
  //     for (int j = 0; j < ny; j++){
  //       printf("%6.2f ", out[i][j]);
  //     }
  //     printf("\n");
  //   }
  // }
  if(rank == MASTER){
    //printf("FILE OUTPUT\n");
    output_image(OUTPUT_FILE, nx, ny, out);
  }
   /* don't forget to tidy up when we're done */
   MPI_Finalize();

  /* free up allocated memory */
  for(ii=0;ii<local_nrows;ii++) {
   free(u[ii]);
   free(w[ii]);
  }
  free(u);
  free(w);
  free(sendbuf);
  free(recvbuf);
  free(printbuf);

  /* and exit the program */
  return EXIT_SUCCESS;
}

int calc_nrows_from_rank(int rank, int size, int nx)
{
  int nrows;

  nrows = nx / size;       /* integer division */
  if ((nx % size) != 0) {  /* if there is a remainder */
    if (rank == size - 1)
      nrows += nx % size;  /* add remainder to last rank */
  }

   return nrows;
}

// Get the current time in seconds since the Epoch
double wtime(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec*1e-6;
}
