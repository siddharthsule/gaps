#ifndef SPLITTINGS_CUH_
#define SPLITTINGS_CUH_

// #include "qcd.cuh"

// Splitting functions as a class - Here be dragons!

class Kernel {
 public:
  int flavs[3];
  __device__ virtual double Value(double z, double y) = 0;
  virtual ~Kernel() {}
};

class Pqq : public Kernel {
 public:
  __device__ double Value(double z, double y) override {
    return CF * (2. / (1. - z * (1. - y)) - (1. + z));
  }
  __device__ static Pqq* getInstance() {
    __shared__ Pqq instance;
    return &instance;
  }
};

__global__ void initKernels(Kernel** Kernels) {
  if (threadIdx.x == 0 && blockIdx.x == 0) {
    // Quarks
    for (int i = 0; i < 5; i++) {
      Kernels[i] = Pqq::getInstance();
      Kernels[i]->flavs[0] = i + 1;  // id = i + 1
      Kernels[i]->flavs[1] = i + 1;
      Kernels[i]->flavs[2] = 21;
    }

    // Anti-quarks
    for (int i = 5; i < 10; i++) {
      Kernels[i] = Pqq::getInstance();
      Kernels[i]->flavs[0] = i - 11;  // id = i - 11
      Kernels[i]->flavs[1] = i - 11;
      Kernels[i]->flavs[2] = 21;
    }
  }
}

// Kernel to perform computations
__global__ void computeValues(Kernel** Kernels, double* input, double* output,
                              int n) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {
    double temp = 0;
    double max = 0;
    for (int i = 0; i < 10; i++) {
      temp = Kernels[i]->Value(input[idx], input[idx + 1]);
      if (temp > max) {
        max = temp;
      }
    }
    output[idx] = max;
  }
}

int main() {
  Kernel** deviceKernels;
  cudaMalloc(&deviceKernels, 10 * sizeof(Kernel*));

  // Initialize Kernels on device
  initKernels<<<1, 1>>>(deviceKernels);
  cudaDeviceSynchronize();  // Ensure initialization is completed

  // Input data
  double input[2] = {0.1, 0.2};
  double* deviceInput;
  cudaMalloc(&deviceInput, 2 * sizeof(double));
  cudaMemcpy(deviceInput, input, 2 * sizeof(double), cudaMemcpyHostToDevice);

  // Output data
  double output[1];
  double* deviceOutput;
  cudaMalloc(&deviceOutput, 1 * sizeof(double));

  // Compute values
  computeValues<<<1, 1>>>(deviceKernels, deviceInput, deviceOutput, 1);
  cudaDeviceSynchronize();

  // Copy output back to host
  cudaMemcpy(output, deviceOutput, 1 * sizeof(double), cudaMemcpyDeviceToHost);

  // Print output
  printf("Output: %f\n", output[0]);

  // Free memory
  cudaFree(deviceKernels);
  cudaFree(deviceInput);
  cudaFree(deviceOutput);

  return 0;
}

#endif  // SPLITTINGS_CUH_
