#include <iostream>

__global__ void fx(int N, double A[N][N], double B[N][N]){
    int i = threadIdx.x;
    int j = threadIdx.y;

}