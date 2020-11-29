#include <stdio.h>

#include "../libs/matrix.h"

#define RCSize 9
#define omega 1.5

template <typename Ty>
std::pair<Matrix<Ty>, std::vector<Ty>> Initialize(){
    Matrix<Ty> inita(RCSize, RCSize);
    inita.value = {{12,  1,  5,  1,  1,  2, -4,  1,  2},
                   { 1, 16, -1, -4, -5, -2,  1,  2,  3},
                   { 5, -1, 15, -5,  3,  1, -2,  1, -4},
                   { 1, -4, -5, 10,  3, -3, -1,  4,  1},
                   { 1, -5,  3,  3, 11, -1,  4,  1,  1},
                   { 2, -2,  1, -3, -1, 15, -5,  2,  5},
                   {-4,  1, -2, -1,  4, -5, 15,  4, -4},
                   { 1,  2,  1,  4,  1,  2,  4, 11, -1},
                   { 2,  3, -4,  1,  1,  5, -4,  -1, 15}};
    std::vector<Ty> initb = {21, 11, 13, 6, 18, 14, 8, 25, 18};

    std::pair<Matrix<Ty>, std::vector<Ty>> initab(inita, initb);

    return initab;
}

int main(){
    auto Ab = Initialize<double>();
    std::vector<double> x(RCSize, 0.0);
    std::vector<double> r = Ab.second - dot(Ab.first, x);
    std::vector<double> p(RCSize);
    _CGReturn<double> CG = {x, r, p, 0.0};

    for (int i = 0; i < 10000; i++){
        CG = Conjugate_Gradient(Ab.first, CG, i);
        if(CG.rho == 0.0){
            std::cout << i << std::endl;
            break;
        }
    }

    std::cout << CG.x << std::endl;

    return 0;
}