#include <stdio.h>

#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>

#include "..\libs\utils.h"

#define REPEAT 100

template <typename Ty = double>
inline Ty fx1(Ty x1, Ty x2){
    return x1*x1 + x2*x2 - 2;
}

template <typename Ty = double>
inline Ty fx2(Ty x1, Ty x2){
    return x1 -x2*x2;
}

template <typename Ty>
std::vector<Ty> fx(const std::vector<Ty>& x){
    auto f1 = fx1(x[0], x[1]);
    auto f2 = fx2(x[0], x[1]);
    std::vector<Ty> Ans = {f1, f2};

    return Ans;
}

template <typename Ty = double>
Matrix<Ty> Jac(const std::vector<Ty>& x){
    Matrix<Ty> Jac(2, 2);
    Jac[0][0] = (fx1(x[0]+1e-7, x[1]) - fx1(x[0]-1e-7, x[1])) / 2e-7;
    Jac[0][1] = (fx1(x[0], x[1]+1e-7) - fx1(x[0], x[1]-1e-7)) / 2e-7;
    Jac[1][0] = (fx2(x[0]+1e-7, x[1]) - fx2(x[0]-1e-7, x[1])) / 2e-7;
    Jac[1][1] = (fx2(x[0], x[1]+1e-7) - fx2(x[0], x[1]-1e-7)) / 2e-7;

    return Jac;
}

int main() {
    std::vector<double> x = {sqrt(2.0), sqrt(2.0)};
    std::vector<double> x1_log;
    std::vector<double> x2_log;
    x1_log.push_back(x[0]);
    x2_log.push_back(x[1]);

    for(int i=0; i<REPEAT; i++){
        std::vector<double> f = fx(x);
        Matrix<double> J = Jac(x);
        J = Inv(J);
        x = x - dot(J, f);
        x1_log.push_back(x[0]);
        x2_log.push_back(x[1]);
        if(fabs((x1_log[i-1]-x[0]) < 1e-20) && (fabs(x2_log[i-1]-x[1]) < 1e-20)){
            break;
        }
    }

    std::cout << "x1: " << x1_log << std::endl;
    std::cout << "x2: " << x2_log << std::endl;

    // create graph
    FILE* gp;
    int si = x1_log.size();
    double plotx1 = 0;
    double plotx2 = 0;

    gp = _popen("gnuplot", "w");
    fprintf(gp, "unset key\n");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set output \'1_3_1_jacobi.png\'\n");  // 保存ファイル名
    fprintf(gp, "set xrange[1e-25:100.0]\n");
    fprintf(gp, "set yrange[1e-25:10.0]\n");
    fprintf(gp, "set xlabel \"x1 Error\"\n");
    fprintf(gp, "set ylabel \"x2 Error\"\n");
    for (int i = 0; i < si; i++) {
        plotx1 = fabs(x1_log[i] - x[0]);
        plotx2 = fabs(x2_log[i] - x[1]);
        fprintf(gp, "set label right at %g, %g \"step %d\"\n", plotx1, plotx2, i);
    }
    fprintf(gp, "set logscale x\n");
    fprintf(gp, "set logscale y\n");
    fprintf(gp, "plot \"-\" with points pt 6 \n");

    for (int i = 0; i < si; i++) {
        plotx1 = fabs(x1_log[i] - x[0]);
        plotx2 = fabs(x2_log[i] - x[1]);
        fprintf(gp, "%g, %g\n", plotx1, plotx2);
    }
    std::cout << "Repeat Num: " << si-1 << std::endl;

    fprintf(gp, "e\n");
    fprintf(gp, "set output\n");
    fprintf(gp, "set terminal windows\n");
    fflush(gp);
    _pclose(gp);

    return 0;
}