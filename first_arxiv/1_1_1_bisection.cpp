#include <stdio.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>

#include "..\libs\utils.h"

#define REPEAT 100


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

// step
std::vector<double> step(double& x1, double& x2, double& x) {
    x = (x1 + x2) / 2;
    double fx1 = func(x1);
    double fx2 = func(x2);
    double fx = func(x);

    if (fx > 0) {
        if (fx1 > 0) {
            x1 = x;
        } else if (fx2 > 0) {
            x2 = x;
        } else {
            std::cout
                << "Method Error: A constraints of Bisection is not completed"
                << std::endl;
        }
    } else if (fx < 0) {
        if (fx1 < 0) {
            x1 = x;
        } else if (fx2 < 0) {
            x2 = x;
        } else {
            std::cout
                << "Method Error: A constraints of Bisection is not completed"
                << std::endl;
        }
    }

    std::vector<double> v = {fx1, fx2};
    return v;
}

int main() {
    // initialize
    double x1 = -2.0;  // x_L (step 0)
    double x2 = 0.0;   // x_R (step 0)
    double x = 0.0;

    double fx = 0.0;
    double fx1 = 0.0;
    double fx2 = 0.0;
    double div = 0.0;
    std::vector<double> is_break_doub = {0.0, 0.0};

    std::vector<double> divAll;
    divAll.push_back(x);

    // Repeat to get answer
    for (int i = 0; i < 100; i++) {
        is_break_doub = step(x1, x2, x);

        if (is_break_doub[0] * is_break_doub[1] >= 0) {
            divAll.push_back(x);
            break;
        }
        if (fabs(is_break_doub[0] - is_break_doub[1]) < 1e-15) {
            divAll.push_back(x);
            break;
        }

        divAll.push_back(x);
    }

    std::cout << divAll << std::endl;
    std::cout << "Answer: " << std::setprecision(20) << x << std::endl;

    // create Graph
    FILE* gp;

    gp = _popen("gnuplot", "w");
    fprintf(gp, "unset key\n");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set output \'1_1_bisection.png\'\n");
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
        fprintf(gp, "%d, %g\n", i, fabs(divAll[i] - x));
    }
    std::cout << "Repeat Num: " << count << std::endl;

    fprintf(gp, "e\n");
    fprintf(gp, "set output\n");
    fprintf(gp, "set terminal windows\n");
    fflush(gp);
    _pclose(gp);

    return 0;
}