#include <stdio.h>

#include "../libs/matrix.h"

#define RCSize 9
#define REPEAT 2500

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
    double omega = 0.0;
    std::vector<double> ans(RCSize, 1.0);
    std::vector<std::vector<double>> res_memo_omega = {};
    std::vector<std::vector<double>> ans_memo_all = {};

    for (int o = 1; o < 4; o++){
        int count = 0;
        omega += 0.5;
        std::cout << "Omega: " << omega << std::endl;

        auto Ab = Initialize<double>();
        std::vector<double> x(RCSize, 0.0);
        std::vector<double> res_memo = {};
        std::vector<double> ans_memo = {};

        auto l2norm = L2norm(x - ans);
        auto res = Resnorm(Ab, x);

        ans_memo.push_back(l2norm);
        res_memo.push_back(res);

        for (int i = 0; i < 10000; i++){
            count++;

            x = SOR_Step_pair(Ab, x, omega);

            res = Resnorm(Ab, x);
            l2norm = L2norm(x - ans);

            ans_memo.push_back(l2norm);
            res_memo.push_back(res);

            if(res < 1e-10){
                break;
            }
        }

        ans_memo_all.push_back(ans_memo);
        res_memo_omega.push_back(res_memo);

        std::cout << "Count: " << count << std::endl;
        std::cout << "    x: " << x << std::endl;
    }

    // create graph
    FILE* gp;

    gp = _popen("gnuplot -persist", "w");
    fprintf(gp, "unset key\n");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set output \'2_1_2_sor.png\'\n");
    fprintf(gp, "set xrange[0:%d]\n", REPEAT);
    fprintf(gp, "set yrange[1e-15:%lf]\n", 10.0);
    fprintf(gp, "set xlabel \"Time\"\n");
    fprintf(gp, "set ylabel \"Error\"\n");
    fprintf(gp, "set logscale y\n");
    fprintf(gp, "set key right top \n");
    fprintf(gp, "plot '-' with points pt 6 title 'omega=0.5' ,\
                      '-' with points pt 6 lt rgb 'blue' title 'omega=1.0' ,\
                      '-' with points pt 6 lt rgb 'red' title 'omega=1.5' \n");

    int si_0 = res_memo_omega[0].size();
    int si_1 = res_memo_omega[1].size();
    int si_2 = res_memo_omega[2].size();

    for (int i = 0; i < si_0; i++) {
        fprintf(gp, "%d, %g\n", i, res_memo_omega[0][i]);
    }
    fprintf(gp, "e\n");
    for (int i = 0; i < si_1; i++){
        fprintf(gp, "%d, %g\n", i, res_memo_omega[1][i]);
    }
    fprintf(gp, "e\n");
    for (int i = 0; i < si_2; i++){
        fprintf(gp, "%d, %g\n", i, res_memo_omega[2][i]);
    }

    fprintf(gp, "e\n");
    fprintf(gp, "set output\n");
    fprintf(gp, "set output \'2_1_3_sor.png\'\n");
    fprintf(gp, "plot '-' with points pt 6 title 'omega=0.5' ,\
                      '-' with points pt 6 lt rgb 'blue' title 'omega=1.0' ,\
                      '-' with points pt 6 lt rgb 'red' title 'omega=1.5' \n");

    for (int i = 0; i < si_0; i++) {
        fprintf(gp, "%d, %g\n", i, ans_memo_all[0][i]);
    }
    fprintf(gp, "e\n");
    for (int i = 0; i < si_1; i++){
        fprintf(gp, "%d, %g\n", i, ans_memo_all[1][i]);
    }
    fprintf(gp, "e\n");
    for (int i = 0; i < si_2; i++){
        fprintf(gp, "%d, %g\n", i, ans_memo_all[2][i]);
    }

    fprintf(gp, "set output\n");
    fprintf(gp, "set terminal windows\n");
    fflush(gp);
    _pclose(gp);

    return 0;
}