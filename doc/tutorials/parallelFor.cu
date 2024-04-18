// Adapted from https://github.com/olcf-tutorials/vector_addition_cuda
#include <stdio.h>

// Large array, 2^20
#define N 1048576

// Kernel that adds the element
// Global = Called on Host, Ran on Device
__global__ void add_vectors(double *a, double *b, double *out) {
  /*
     Indexing within Grids
     ---------------------

     Can use dim3 variables to get the index of a block/thread,
     as well as the size of a grid/block

     dim3 blockIdx - unique
     dim3 threadIdx - unique in own block
     dim3 gridDim - size of grid
     dim3 blockDim - size of block

     A useful indexing command is
     blockDim.x * blockIdx.x + threadIdx.x
     - blockDim.x * blockIdx.x allows going from 0 to block_size one at a time
     - + ThreadIdx.x allows acces to threads in a block
   */

  // Shared Memory Example
  /*
   __shared__ int shared_array[N];
   shared_array[i] = in[i] // Each Thread writes to one element of s_a
   */

  int id = blockDim.x * blockIdx.x + threadIdx.x;
  if (id < N) out[id] = a[id] + b[id];
}

// Main Program
int main() {
  // Time the execution
  cudaEvent_t start, stop;
  float time;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  // Number of bytes to allocate for N Doubles
  size_t bytes = N * sizeof(double);

  // Allocate memory for arrays A, B and Out on host
  // malloc casts to void*, use (double*) to match the pointers
  double *h_a = (double *)malloc(bytes);
  double *h_b = (double *)malloc(bytes);
  double *h_out = (double *)malloc(bytes);

  // Allocate memory for arrays A, B and Out on Device
  double *d_a, *d_b, *d_out;
  cudaMalloc(&d_a, bytes);
  cudaMalloc(&d_b, bytes);
  cudaMalloc(&d_out, bytes);

  // Fill A and B
  for (int i = 0; i < N; i++) {
    h_a[i] = 1.0;
    h_b[i] = 2.0;
  }

  // Copy data from host to device
  cudaMemcpy(d_a, h_a, bytes, cudaMemcpyHostToDevice);
  cudaMemcpy(d_b, h_b, bytes, cudaMemcpyHostToDevice);

  // Grid and Block Size
  int block_size = 256;
  int grid_size = ceil(float(N) / block_size);

  // Start Timer
  cudaEventRecord(start, 0);

  // Launch kernel
  add_vectors<<<grid_size, block_size>>>(d_a, d_b, d_out);

  // End Timer
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);

  // Get Elapsed Time
  cudaEventElapsedTime(&time, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  // Copy data from device to host
  cudaMemcpy(h_out, d_out, bytes, cudaMemcpyDeviceToHost);

  // Free Memory
  free(h_a);
  free(h_b);
  free(h_out);
  cudaFree(d_a);
  cudaFree(d_b);
  cudaFree(d_out);

  // Print if everything works
  printf("\n---------------------------\n");
  printf("SUCCESS\n");
  printf("%d\n", time);
  printf("---------------------------\n");
}
