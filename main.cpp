#include <stdio.h>

#include <cmath>
#include <iostream>
#include <vector>

#include "libs\utils.h"

#define Epochs 30

int main() {
    Matrix<double> A(2, 2, 1);
    A[0][0] = 1;
    A[0][1] = 2;
    A[1][0] = 3;
    A[1][1] = 4;

    //std::vector<double> b = {2, 4};

    auto inv_A = Inv(A);

    //auto x = MatEqu(A, b);

    std::cout << A << std::endl;
    std::cout << inv_A << std::endl;
    //std::cout << x << std::endl;

    return 0;
}