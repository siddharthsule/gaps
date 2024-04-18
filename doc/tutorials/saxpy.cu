// From https://developer.nvidia.com/blog/easy-introduction-cuda-c-and-c/
#include <stdio.h>

// Large array, 2^20
#define N 1048576

// Kernel function to add the elements of two arrays
__global__ void saxpy(int n, float a, float *x, float *y) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) y[i] = a * x[i] + y[i];
}

int main(void) {
  // Host input vectors
  float *h_x, *h_y, *d_x, *d_y;
  h_x = (float *)malloc(N * sizeof(float));
  h_y = (float *)malloc(N * sizeof(float));

  // Device input vectors
  cudaMalloc(&d_x, N * sizeof(float));
  cudaMalloc(&d_y, N * sizeof(float));

  // Initialize x and y arrays on the host
  for (int i = 0; i < N; i++) {
    h_x[i] = 1.0f;
    h_y[i] = 2.0f;
  }

  // Copy data from host to device
  cudaMemcpy(d_x, x, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_y, y, N * sizeof(float), cudaMemcpyHostToDevice);

  // Perform SAXPY on 1M elements
  saxpy<<<(N + 255) / 256, 256>>>(N, 2.0f, d_x, d_y);

  // Copy data from device to host
  cudaMemcpy(y, d_y, N * sizeof(float), cudaMemcpyDeviceToHost);

  // Check for errors (all values should be 4.0f)
  float maxError = 0.0f;
  for (int i = 0; i < N; i++) maxError = max(maxError, abs(y[i] - 4.0f));
  printf("Max error: %f\n", maxError);

  // Cleanup
  cudaFree(d_x);
  cudaFree(d_y);
  free(h_x);
  free(h_y);
}
