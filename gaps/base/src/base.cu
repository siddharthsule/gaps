#include "base.cuh"

// sync device and check for errors
void sync_gpu_and_check(const char *operation) {
  // synchronize with the device
  cuda_device_synchronize();

  // check for an error
  cuda_error_t error = cuda_get_last_error();
  if (error != cuda_success) {
    // print the cuda error message
    std::cerr << "cuda error @" << operation << ": "
              << cuda_get_error_string(error) << std::endl;

    // abort the program
    std::exit(exit_failure);
  }
}

// debug messages
__host__ __device__ void debug_msg(const char *message) {
  if (debug) {
    printf("debug: %s\n", message);
  }
}