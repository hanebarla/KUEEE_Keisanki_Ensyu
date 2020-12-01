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
    int count = 0;
    auto Ab = Initialize<double>();
    std::vector<double> x(RCSize, 0.0);
    std::vector<double> ans(RCSize, 1.0);
    std::vector<double> res_memo = {};
    std::vector<double> ans_memo = {};

    auto res = Resnorm(Ab, x);
    auto l2norm = L2norm(x - ans);

    res_memo.push_back(res);
    ans_memo.push_back(l2norm);

    for (int i = 0; i < 10000; i++){
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

    std::cout << "Count: " << count << std::endl;
    std::cout <<   "res: " << res << std::endl;
    std::cout << "    x: " << x << std::endl;

        // create graph
    FILE* gp;

    gp = _popen("gnuplot", "w");
    fprintf(gp, "unset key\n");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set output \'2_1_2_gs.png\'\n");
    fprintf(gp, "set xrange[0:%d]\n", REPEAT);
    fprintf(gp, "set yrange[1e-15:%lf]\n", 10.0);
    fprintf(gp, "set xlabel \"Time\"\n");
    fprintf(gp, "set ylabel \"Error\"\n");
    fprintf(gp, "set logscale y\n");
    fprintf(gp, "plot \"-\" with points pt 6 \n");

    int si = res_memo.size();

    for (int i = 0; i < si; i++) {
        fprintf(gp, "%d, %g\n", i, res_memo[i]);
    }


    fprintf(gp, "e\n");
    fprintf(gp, "set output\n");
    fprintf(gp, "set output \'2_1_3_gs.png\'\n");
    fprintf(gp, "plot \"-\" with points pt 6 \n");

    si = ans_memo.size();
    for(int i = 0; i < si; i++){
        fprintf(gp, "%d, %g\n", i, ans_memo[i]);
    }

    fprintf(gp, "e\n");
    fprintf(gp, "set output\n");
    fprintf(gp, "set terminal windows\n");
    fflush(gp);
    _pclose(gp);

    return 0;
}