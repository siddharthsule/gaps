int main(void) {
  // DOES NOT RUN

  // Declare Variables
  // h_ = on host, d_ = on device
  int *h_c, d_c;

  // Allocate memory on the device
  // cudaMalloc( Location of the Memory, Size of the Memory )
  cudaMalloc((void**)&d_c, sizeof(int));

  // If h_c initialised, copy info from h_c to d_c
  // cudaMemcpy( destination, host, numBytes, Direction )
  cudaMemcpy(d_c, h_c, sizeof(int) cudaMemcpyHostToDevice);

  // Kernel Configuration Parameters
  dim3 grid_size(1);
  dim3 block_size(1);

  // Launch the Kernel
  kernel<<<grid_size, block_size>>>(...);

  // Copy data back to host
  // cudaMemcpy( destination, device, numBytes, Direction )
  cudaMemcpy(h_c, d_c, sizeof(int), cudaMemcpyDeviceToHost);

  // Deallocate Memory
  cudaFree(d_c);
  free(h_c);

  return 0;
}
