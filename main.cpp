#include <iostream>
// Kernel detenition
template<typename N>
__global__ void MatAdd(float A[N][N], float B[N][N], float C[N][N]){
    int i = threadIdx.x;
    int j = threadIdx.y;
    C[i][j] = A[i][j] + B[i][j];
}

int main(){
    float *A, *B, *C;
    int N = 100;
    cudaMalloc((void**)&A, N*N*sizeof(float));
    cudaMalloc((void**)&B, N*N*sizeof(float));
    cudaMalloc((void**)&C, N*N*sizeof(float));
    float *a = malloc(N*N*sizeof(float));
    float *b = malloc(N*N*sizeof(float));
    float *c = malloc(N*N*sizeof(float));
    cudaMemcpy(A, a, N*N*sizeof(*A), cudaMemcpyHostToDevice);
    cudaMemcpy(B, b, N*N*sizeof(*B), cudaMemcpyHostToDevice);

    // Kernel invocation with one block of N*N*1thread
    int numBlocks = 1;
    dim3 threadsPerBlock(N, N);

    MatAdd <<< numBlocks, threadsPerBlocks >>>(A, B, C);

    cudaMemcpy(c, C, N*N*sizeof(*c), cudaMemcpyHostToDevice);

    cudaFree(A);
    cudaFree(B);
    cudaFree(C);
    std::cout << "Done" << std::endl;
}