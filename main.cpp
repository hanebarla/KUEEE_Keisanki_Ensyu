#include <stdio.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>

#include "libs/utils.h"
#include "libs/matrix.h"

int main(){
    Matrix <double>a(2, 2);
    a = {{1, 1}, {1, 1}};
    std::cout << a << std::endl;
    return 0;
}