#include "base.cuh"

void syncGPUAndCheck(const char *operation) {
  // synchronize with the device
  cudaDeviceSynchronize();

  // check for an error
  cudaError_t error = cudaGetLastError();
  if (error != cudaSuccess) {
    // print the CUDA error message
    std::cerr << "CUDA error @" << operation << ": "
              << cudaGetErrorString(error) << std::endl;

    // abort the program
    std::exit(EXIT_FAILURE);
  }
}