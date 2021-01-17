#include <stdio.h>

#include "../libs/matrix.h"

#define PI 3.14159265358979323846
#define TAU 2*PI
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
    double end = TAU;
    std::vector<double> ps = {3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17 ,18};
    std::vector<double> MaxErr = {};
    std::vector<double> Log2ME = {};

    for(auto it = ps.begin(); it != ps.end(); ++it){
        double per = std::pow(2.0, *it);
        double h = TAU / per;
        double t = 0.0;

        std::vector<double> B = {0, 0, 1};
        std::vector<double> vc = {0, -1, 0};
        std::vector<double> va = calcVa(t);

        double er = L2norm(vc - va);
        double max_er = er;

        while (t <= end){
            t += h;
            vc = vc + h*Func(vc, B);
            va = calcVa(t);
            er = L2norm(vc - va);

            if(er > max_er){
                max_er = er;
            }
        }

        MaxErr.push_back(max_er);
        Log2ME.push_back(log2(max_er));
        std::cout << "P: " << *it << ", Max Error: " << std::setprecision(16) << log2(max_er) << std::endl;
    }

    // 最小二乗法
    double sum1 = double(ps.size());
    double sigmax = std::accumulate(ps.begin(), ps.end(), 0);
    double sigmay = std::accumulate(Log2ME.begin(), Log2ME.end(), 0);
    auto xsq = ps * ps;
    double sigmaxsq = std::accumulate(xsq.begin(), xsq.end(), 0);
    auto xy = ps * Log2ME;
    double sigmaxy = std::accumulate(xy.begin(), xy.end(), 0);

    double a = (sigmax * sigmay - sum1 * sigmaxy) / (sigmax * sigmax - sum1 * sigmaxsq);
    double b = (sigmax * sigmaxy - sigmay * sigmaxsq) / (sigmax * sigmax - sum1 * sigmaxsq);
    std::cout << "a: " << a << std::endl;
    std::cout << "b: " << b << std::endl;


    FILE* gp;

    gp = _popen("gnuplot", "w");
    fprintf(gp, "unset key\n");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set output \'3_1_2_euler_p_err.png\'\n");
    fprintf(gp, "set xrange[%d:%d]\n", 1, 20);
    fprintf(gp, "set yrange[%d:%d]\n", -16, 4);
    fprintf(gp, "set xlabel \"p\"\n");
    fprintf(gp, "set ylabel \"log_2(Er)\"\n");
    fprintf(gp, "plot %g*x+%g lc 'red', \"-\" with points pt 6 ps 1.5 linecolor rgb '#3333bb' \n", a, b);

    for (int i=0; i < int(ps.size()); i++) {
        fprintf(gp, "%g, %g\n", ps[i], log2(MaxErr[i]));
    }


    fprintf(gp, "e\n");
    fprintf(gp, "set output\n");
    fprintf(gp, "set terminal windows\n");
    fflush(gp);
    _pclose(gp);

    return 0;
}