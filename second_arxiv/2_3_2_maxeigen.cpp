#include <stdio.h>

#include "../libs/matrix.h"

#define RCSize 9
#define REPEAT 100

template <typename Ty>
std::pair<Matrix<Ty>, std::vector<Ty>> Initialize(){
    Matrix<Ty> inita(RCSize, RCSize);
    inita.value = {{ 9,  1,  5,  1,  1,  2, -4,  1,  2},
                   { 1, 16, -1, -4, -5, -2,  1,  2,  3},
                   { 5, -1, 15, -5,  3,  1, -2,  1, -4},
                   { 1, -4, -5, 10,  3, -3, -1,  4,  1},
                   { 1, -5,  3,  3, 11, -1,  4,  1,  1},
                   { 2, -2,  1, -3, -1, 15, -5,  2,  5},
                   {-4,  1, -2, -1,  4, -5, 15,  4, -4},
                   { 1,  2,  1,  4,  1,  2,  4, 11, -1},
                   { 2,  3, -4,  1,  1,  5, -4, -1, 15}};
    std::vector<Ty> initb = {18, 11, 13, 6, 18, 14, 8, 25, 18};

    std::pair<Matrix<Ty>, std::vector<Ty>> initab(inita, initb);

    return initab;
}

int main(){
    int count = 0;
    auto Ab = Initialize<double>();
    std::vector<double> x(RCSize, 1.0);
    std::vector<double> y(RCSize, 0.0);
    double l2norm = 0.0;

    for (int i = 0; i < REPEAT; i++){
        y = dot(Ab.first, x);
        l2norm = L2norm(y);
        std::cout << l2norm << std::endl;
        x = y / l2norm;
    }

    std::cout << "L2Norm: " << l2norm << std::endl;

    return 0;
}