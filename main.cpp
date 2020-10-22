#include <stdio.h>

#include <cmath>
#include <iostream>
#include <vector>

#include "libs\utils.h"

#define Epochs 2000

int main() {
    // initialize
    auto A = Identity<long double>(10);
    std::vector<long double> x(10, 1);
    std::vector<long double> y(10, 0);
    A *= 2.0;
    for (int i = 0; i < A.row; i++) {
        for (int j = 0; j < A.colum; j++) {
            if (abs(i - j) == 1) {
                A[i][j] = -1;
            }
        }
    }

    // repeat
    for (int i=0; i<Epochs; i++){
        y = dot(A, x); // ここがおかしそう
        double l2 = L2norm(y);
        x = y / l2;
        std::cout << "l2: " << l2 << std::endl;
    }

    return 0;
}