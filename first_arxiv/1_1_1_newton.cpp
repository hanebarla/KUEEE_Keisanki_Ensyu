#include <stdio.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>

#include "..\libs\matrix.h"

#define REPEAT 30

// f(x)
double func(double x) {
    double tmp = 0.0;
    tmp += pow(x, 5);
    tmp -= 3.0 * pow(x, 4);
    tmp += pow(x, 3);
    tmp += 5.0 * pow(x, 2);
    tmp -= 6 * x;
    tmp += 2;

    return tmp;
}

int main() {
    // initialize
    double x = -1.0;  // x(step 0)
    double xp = 0.0;
    double fx = 0.0;
    double xp_ep = 0.0;
    double xp_em = 0.0;
    double div = 0.0;
    std::vector<double> divAll;
    divAll.push_back(x);

    // Repeat to get answer
    for (int i = 0; i < REPEAT; i++) {
        fx = func(x);
        xp_ep = func(x + 1e-7);
        xp_em = func(x - 1e-7);
        div = (xp_ep - xp_em) / 2e-7;
        xp = x - (fx / div);
        divAll.push_back(xp);

        if (fabs(x - xp) < 1e-15) {
            break;
        }

        x = xp;
    }
    std::cout << divAll << std::endl;
    std::cout << "Answer: " << std::setprecision(20) << x << std::endl;

    // create graph
    FILE* gp;

    gp = _popen("gnuplot", "w");
    fprintf(gp, "unset key\n");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set output \'1_1_newton.png\'\n");
    fprintf(gp, "set xrange[0:%d]\n", REPEAT);
    fprintf(gp, "set yrange[1e-25:%lf]\n", 10.0);
    fprintf(gp, "set xlabel \"Time\"\n");
    fprintf(gp, "set ylabel \"Error\"\n");
    fprintf(gp, "set logscale y\n");
    fprintf(gp, "plot \"-\" with points pt 6 \n");

    int si = divAll.size();
    int count = 0;

    for (int i = 0; i < si; i++) {
        double diff = fabs(divAll[i] - x);
        if(diff > 0){
            count += 1;
        }
        fprintf(gp, "%d, %g\n", i, diff);
    }
    std::cout << "Repeat Num: " << count << std::endl;

    fprintf(gp, "e\n");
    fprintf(gp, "set output\n");
    fprintf(gp, "set terminal windows\n");
    fflush(gp);
    _pclose(gp);

    return 0;
}