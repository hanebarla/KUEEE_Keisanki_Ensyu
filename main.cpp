#include <cmath>
#include <iostream>
#include <vector>

#include "libs\utils.h"
#include "libs\matrix.h"

int main() {
    Matrix<double> x(11, 3);

    std::cout<< x << std::endl;
    x.T();
    std::cout << x << std::endl;

    return 0;
}