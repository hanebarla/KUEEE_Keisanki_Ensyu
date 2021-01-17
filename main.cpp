#include <stdio.h>

#include "libs/matrix.h"

#define RCSize 6
#define omega 1.5

template <typename Ty>
std::pair<Matrix<Ty>, std::vector<Ty>> Initialize(){
    Matrix<Ty> inita(RCSize, RCSize);
    //phi0, 1, 2, 3, 4, 5
    inita.value = {{-4, 1, 0, 0, 0, 0},
                   {1, -4, 1, 0, 0, 0},
                   {0, 1, -4, 2, 0, 0},
                   {0, 0, 1, -4, 1, 0},
                   {0, 0, 0, 1, -4, 1},
                   {0, 0, 0, 0, 0, 1}};
    std::vector<Ty> initb = {0, 0, -100, -100, 0, 50};

    std::pair<Matrix<Ty>, std::vector<Ty>> initab(inita, initb);

    return initab;
}

int main(){
    int count = 0;
    auto Ab = Initialize<float>();
    std::vector<float> x(RCSize, 0.0);
    std::vector<float> ans(RCSize, 1.0);
    std::vector<float> res_memo = {};
    std::vector<float> ans_memo = {};

    auto res = Resnorm(Ab, x);
    auto l2norm = L2norm(x - ans);

    res_memo.push_back(res);
    ans_memo.push_back(l2norm);

    for (long int i = 0; i < 1000000; i++){
        count++;

        x = Gauss_Seidel_Step_pair(Ab, x);

        res = Resnorm(Ab, x);
        l2norm = L2norm(x - ans);

        res_memo.push_back(res);
        ans_memo.push_back(l2norm);

        if(res < 1e-10){
            break;
        }
    }

    std::cout << " Count: " << count << std::endl;
    std::cout << "L2Norm: " << l2norm << std::endl;
    std::cout << "   res: " << res << std::endl;
    std::cout << "     x: " << x << std::endl;

    return 0;
}