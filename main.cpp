#include <cmath>
#include <iostream>
#include <vector>

#include "libs\utils.h"
#include "libs\matrix.h"

int main() {
    Matrix<double> x(100, 100);
    Matrix<double> f(100, 100);

    auto a = multiply(f, x);
    std::cout<< a.value[1] << std::endl;

    return 0;
}