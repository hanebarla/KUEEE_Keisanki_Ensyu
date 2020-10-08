#include <cmath>
#include <iostream>
#include <vector>

#include "libs\utils.h"

double func(double x) {
    double tmp = 0.0;
    tmp += pow(x, 5);
    tmp -= 3.0 * pow(x, 4);
    tmp += pow(x, 3);
    tmp += 5.0 * pow(x, 2);
    tmp += 2;

    return tmp;
}

int main() {
    double x = -1.0;
    double xp = 0.0;
    double fx = 0.0;
    double xp_ep = 0.0;
    double xp_em = 0.0;
    double div = 0.0;
    std::vector<double> divAll;
    divAll.push_back(x);

    for (int i = 0; i < 100; i++) {
        fx = func(x);
        xp_ep = func(x + 1e-7);
        xp_em = func(x - 1e-7);
        div = (xp_ep - xp_em) / 2e-7;
        xp = x - (fx / div);
        divAll.push_back(xp);

        if (fabs(x - xp) < 1e-8) {
            x = xp;
            break;
        }
        x = xp;
    }
    std::cout << divAll << std::endl;
    std::cout << "Answer: " << x << std::endl;

    return 0;
}