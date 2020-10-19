#include <cmath>
#include <iostream>
#include <vector>

#include "libs\utils.h"
#include "libs\matrix.h"

int main() {
    auto E = Identity<float>(5);
    std::cout << E << std::endl;

    return 0;
}