#include <stdio.h>

#include <cmath>
#include <iostream>
#include <vector>

#include "libs\utils.h"

#define Epochs 30

int main() {
    // initialize
    auto A = Identity<long double>(10);
    std::vector<long double> x(10, 1);
    std::vector<long double> y(10, 0);
    A *= 2.0;
    for (int i = 0; i < A.row; i++) {
        for (int j = 0; j < A.colum; j++) {
            if (abs(i - j) == 1) {
                A[i][j] = -1;
            }
        }
    }
    double l2 = 0.0;
    std::vector<long double> logger;
    logger.push_back(l2);

    // repeat
    for (int i=0; i<Epochs; i++){
        y = dot(A, x);
        l2 = L2norm(y);
        x = y / l2;
        logger.push_back(l2);
    }

    std::cout << logger << std::endl;
    std::cout << "||y||_2: " << l2 << std::endl;

    // create graph
    FILE* gp;
    gp = _popen("gnuplot", "w");
    fprintf(gp, "unset key\n");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set output \'1_2_1_eigen.png\'\n");
    fprintf(gp, "set xrange[0:%d]\n", Epochs);
    fprintf(gp, "set yrange[0:%g]\n", 4.0);
    fprintf(gp, "set xlabel \"Time\"\n");
    fprintf(gp, "set ylabel \"Error\"\n");
    fprintf(gp, "plot \"-\" with points pt 6 \n");

    int si = logger.size();

    for (int i = 0; i < si; i++) {
        fprintf(gp, "%d, %g\n", i, double(logger[i]));
    }

    fprintf(gp, "e\n");
    fprintf(gp, "set output\n");
    fprintf(gp, "set terminal windows\n");
    fflush(gp);
    _pclose(gp);

    return 0;
}