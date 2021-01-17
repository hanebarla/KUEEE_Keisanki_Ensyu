#include <stdio.h>

#include "../libs/matrix.h"

#define PI 3.14159265358979323846
#define TAU 2*PI
#define DEL 64.0
#define Q 1.0
#define M 1.0

// 微分方程式の右辺
template<typename Ty>
std::vector<Ty> Func(std::vector<Ty> v, std::vector<Ty> B){
    return Ty(Q / M) * cross(v, B);
}

// 解析解
template<typename Ty>
std::vector<Ty> calcVa(Ty t){
    std::vector<Ty> ans = {-sin(t), -cos(t), 0.0};

    return ans;
}

int main(){
    double end = TAU * 5.0;
    double h = TAU / DEL;
    double t = 0.0;
    std::vector<double> k1 = {0.0, 0.0, 0.0};
    std::vector<double> k2 = {0.0, 0.0, 0.0};
    std::vector<double> k3 = {0.0, 0.0, 0.0};
    std::vector<double> k4 = {0.0, 0.0, 0.0};

    std::vector<double> B = {0, 0, 1};
    std::vector<double> vc = {0, -1, 0};
    std::vector<double> va = calcVa(t);

    std::vector<std::vector<double>> mem_vc;
    mem_vc.push_back(vc);
    std::vector<double> mem_er;
    double er = L2norm(vc - va);
    std::cout << "Start Error: " << std::setprecision(16) << er << std::endl;
    mem_er.push_back(er);

    while (t <= end){
        t += h;

        k1 = Func(vc, B);
        k2 = Func(vc + (h/2)*k1, B);
        k3 = Func(vc + (h/2)*k2, B);
        k4 = Func(vc + h*k3, B);
        vc = vc + (k1 + (2.0 * k2) + (2.0 * k3) + k4) * (h / 6);
        va = calcVa(t);

        mem_vc.push_back(vc);
        er = L2norm(vc - va);
        std::cout << "Error: " << std::setprecision(16) << er << std::endl;
        mem_er.push_back(er);
    }

    FILE* gp;

    gp = _popen("gnuplot", "w");
    fprintf(gp, "unset key\n");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set output \'3_1_4_RK_vc.png\'\n");
    fprintf(gp, "set xrange[%d:%d]\n", -3, 3);
    fprintf(gp, "set yrange[%d:%d]\n", -3, 3);
    fprintf(gp, "set xlabel \"x\"\n");
    fprintf(gp, "set ylabel \"y\"\n");
    fprintf(gp, "plot \"-\" with points pt 6 \n");

    int si = mem_vc.size();

    for (int i = 0; i < si; i++) {
        fprintf(gp, "%g, %g\n", mem_vc[i][0], mem_vc[i][1]);
    }


    fprintf(gp, "e\n");
    fprintf(gp, "set output\n");
    fprintf(gp, "set output \'3_1_4_RK_er.png\'\n");
    fprintf(gp, "set xrange[%d:%d]\n", 0, 35);
    fprintf(gp, "set yrange[%g:%g]\n", 1e-9, 1e-4);
    fprintf(gp, "set xlabel \"time\"\n");
    fprintf(gp, "set ylabel \"Error\"\n");
    fprintf(gp, "set logscale y\n");
    fprintf(gp, "plot \"-\" with points pt 6 \n");

    si = mem_er.size();
    for(int i = 0; i < si; i++){
        fprintf(gp, "%g, %g\n", i * h, mem_er[i]);
    }

    fprintf(gp, "e\n");
    fprintf(gp, "set output\n");
    fprintf(gp, "set terminal windows\n");
    fflush(gp);
    _pclose(gp);

    return 0;
}