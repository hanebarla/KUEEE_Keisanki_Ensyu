#include <stdio.h>

#include "../libs/matrix.h"

#define A1 3.0
#define B1 9.0
#define A2 2.0
#define B2 2.0
#define GAMMA 0

// 微分方程式の右辺
template<typename Ty>
std::vector<Ty> Func(std::vector<Ty> var){
    Ty x = var[0];
    Ty y = var[1];
    std::vector<Ty> next_var(2);

    next_var[0] = (A1 - GAMMA*x - B1*y) * x;
    next_var[1] = (-A2 + B2*x)*y;

    return next_var;
}

int main(){
    double end = 20;
    double h = 0.1;
    double t = 0.0;
    std::vector<double> k1 = {0.0, 0.0};
    std::vector<double> k2 = {0.0, 0.0};
    std::vector<double> k3 = {0.0, 0.0};
    std::vector<double> k4 = {0.0, 0.0};

    std::vector<double> var = {2, 1};

    std::vector<std::vector<double>> mem_var;
    mem_var.push_back(var);

    while (t <= end){
        t += h;

        k1 = Func(var);
        k2 = Func(var + (h/2)*k1);
        k3 = Func(var + (h/2)*k2);
        k4 = Func(var + h*k3);
        var = var + (k1 + (2.0 * k2) + (2.0 * k3) + k4) * (h / 6);

        mem_var.push_back(var);
    }

    FILE* gp;

    gp = _popen("gnuplot", "w");
    fprintf(gp, "unset key\n");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set output \'3_2_1_LV_var.png\'\n");
    fprintf(gp, "set xrange[%d:%d]\n", 0, 21);
    fprintf(gp, "set yrange[%d:%d]\n", 0, 5);
    fprintf(gp, "set xlabel \"time\"\n");
    fprintf(gp, "set ylabel \"population\"\n");
    fprintf(gp, "set key right top \n");
    fprintf(gp, "plot \"-\" with points pt 6 title 'Prey' ,\
                      \"-\" with points pt 6 title 'Predator' \n");

    int si = mem_var.size();

    for (int i = 0; i < si; i++) {
        fprintf(gp, "%g, %g\n", i * h, mem_var[i][0]);
    }
    fprintf(gp, "e\n");
    for (int i = 0; i < si; i++) {
        fprintf(gp, "%g, %g\n", i * h, mem_var[i][1]);
    }


    fprintf(gp, "e\n");
    fprintf(gp, "set output\n");
    fprintf(gp, "set terminal windows\n");
    fflush(gp);
    _pclose(gp);

    return 0;
}