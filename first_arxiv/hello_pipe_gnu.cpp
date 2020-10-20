#include <stdio.h>
#include <vector>
#include <cmath> //for sin, cos
#define GNUPLOT_PATH "gnuplot.exe"

int main(){
    int loop = 1000;
    int N = 256;
    int K = 16;

    FILE *gp;
    gp = _popen(GNUPLOT_PATH, "w");
    fprintf(gp, "unset key\n");
    fprintf(gp, "set xrange[0:%d]\n", N);
    fprintf(gp, "set yrange[0:%d]\n", K);
    fprintf(gp, "set zrange[-2:2]\n");
    fprintf(gp, "set xlabel \"Time\"\n");
    fprintf(gp, "set ylabel \"Time\"\n");
    fprintf(gp, "set zlabel \"Amplitude\"\n");

    std::vector<std::vector<double> > x(K, std::vector<double>(N,0));
    int k = 0;
    for(int i = 0; i < loop; i++){
        for(int n = 0; n < N; n++){
            x[k][n] = std::sin(2.0 * std::acos(-1.0) * (n + i * 10) / N);
        }

        fprintf(gp, "splot \"-\" with points \n");
        for(int j = 0; j < K; j++){
            for(int n = 0; n < N; n++){
                fprintf(gp, "%d, %d, %lf\n", n, j, x[(K + k - j)%K][n]);
            }
        }
        fprintf(gp, "e\n");
        fflush(gp);

        k = (k + 1)%K;
    }
    _pclose(gp);
    return 0;
}