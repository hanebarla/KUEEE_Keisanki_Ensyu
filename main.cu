#include <iostream>
// Kernel detenition
__global__ void MatAdd(int N, float *A, float *B, float *C){
    int i = threadIdx.x;
    int j = threadIdx.y;
    C[N*i + j] = A[N*i + j] + B[N*i + j];
}

int main(){
    float *A, *B, *C;
    int N = 100;
    cudaMalloc((void**)&A, N*N*sizeof(float));
    cudaMalloc((void**)&B, N*N*sizeof(float));
    cudaMalloc((void**)&C, N*N*sizeof(float));
    auto *a = malloc(N*N*sizeof(float));
    auto *b = malloc(N*N*sizeof(float));
    auto *c = malloc(N*N*sizeof(float));
    cudaMemcpy(A, a, N*N*sizeof(*A), cudaMemcpyHostToDevice);
    cudaMemcpy(B, b, N*N*sizeof(*B), cudaMemcpyHostToDevice);

    // Kernel invocation with one block of N*N*1thread
    int numBlocks = 1;
    dim3 threadsPerBlock(N, N, 1);

    MatAdd <<< numBlocks, threadsPerBlock >>>(N, (float*)A, (float*)B, (float*)C);

    cudaMemcpy(c, C, N*N*sizeof(*C), cudaMemcpyDeviceToHost);

    cudaFree(A);
    cudaFree(B);
    cudaFree(C);
    std::cout << "Done" << std::endl;
}