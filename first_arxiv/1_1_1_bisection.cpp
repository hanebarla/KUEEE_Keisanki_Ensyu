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
    double x1 = -2.0;
    double x2 = 0.0;
    double x = 0.0;

    double fx = 0.0;
    double fx1 = 0.0;
    double fx2 = 0.0;
    double div = 0.0;
    std::vector<double> is_break_doub = {0.0, 0.0};

    std::vector<double> divAll;
    divAll.push_back(x);

    for (int i = 0; i < 100; i++) {
        is_break_doub = step(x1, x2, x);

        if (is_break_doub[0] * is_break_doub[1] >= 0) {
            divAll.push_back(x);
            break;
        }
        if (fabs(is_break_doub[0] - is_break_doub[1]) < 1e-7) {
            divAll.push_back(x);
            break;
        }

        divAll.push_back(x);
    }

    std::cout << divAll << std::endl;
    std::cout << "Answer: " << x << std::endl;

    return 0;
}