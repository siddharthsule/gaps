#include "base.cuh"

// sync device and check for errors
void sync_gpu_and_check(const char *operation) {
  // synchronize with the device
  cudaDeviceSynchronize();

  // check for an error
  cudaError_t error = cudaGetLastError();
  if (error != cudaSuccess) {
    // print the cuda error message
    std::cerr << "CUDA error @" << operation << ": "
              << cudaGetErrorString(error) << std::endl;

    // abort the program
    std::exit(EXIT_FAILURE);
  }
}

// debug messages
__host__ __device__ void debug_msg(const char *message) {
  if (debug) {
    printf("debug: %s\n", message);
  }
}