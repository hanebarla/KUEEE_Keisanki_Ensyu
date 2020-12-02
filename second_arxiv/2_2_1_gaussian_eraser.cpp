#include <stdio.h>

#include "../libs/matrix.h"

template <typename Ty>
std::pair<Matrix<Ty>, std::vector<Ty>> Initialize(){
    Matrix<Ty> inita(9, 9);
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
    auto Ab = Initialize<float>();
    std::vector<int> Exchangememo(Ab.first.row);
    auto fe_Ab = Forward_easure_pair(Ab, Exchangememo);
    auto bs_Ab = Backward_subsitution_pair(fe_Ab, Exchangememo);
    std::cout << "x_idx: " << Exchangememo << std::endl;
    std::cout << "    x: " << bs_Ab.second << std::endl;
    auto resno = Resnorm(Ab, bs_Ab.second);
    std::cout << "Res Norm: " << resno << std::endl;
    std::cout << "ME: " << machine_epsilon() << std::endl;

    return 0;
}