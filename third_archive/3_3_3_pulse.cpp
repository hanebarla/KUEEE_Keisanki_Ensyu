#include <stdio.h>

#include "../libs/matrix.h"

#define N 100
#define Delx 1.0/N
#define Delt 1e-3
#define c 5
#define Ts 0.1
#define Td 0.1

template<typename Ty>
Ty Fixedcond(const std::vector<Ty>& um, const std::vector<Ty>& u){
    Ty mu = c*c*Delt*Delt/(Delx*Delx);
    Ty nu_nm_kp = -um[N-1] + mu*u[N-2] + 2*(1-mu)*u[N-1];

    return nu_nm_kp;
}

template<typename Ty>
Ty Freecond(const std::vector<Ty>& um, const std::vector<Ty>& u){
    Ty mu = c*c*Delt*Delt/(Delx*Delx);
    Ty nu_nm_kp = -um[N-1] + mu*u[N-2] + (2-mu)*u[N-1];

    return nu_nm_kp;
}

template<typename Ty>
Ty Smooth(Ty t){
    Ty alpha = 0.01;
    Ty formmer = (tanh((t - Ts)/alpha) + 1)/2;
    Ty latter = (tanh((t - Ts - Td)/alpha) + 1)/2;

    return formmer - latter;
}

template<typename Ty>
Ty Pulse(double t){
    if((t >= Ts) && (t <= (Ts + Td))){
        return Ty(1.0);
    }
    else{
        return Ty(0.0);
    }
}

template<typename Ty>
std::vector<Ty> step(const std::vector<Ty>& um, const std::vector<Ty>& u, Ty t, bool smoothing=true, bool fiexd=true){
    std::vector<Ty> next_u = u;
    std::vector<Ty> u_i_km = um;
    Ty mu = c*c*Delt*Delt/(Delx*Delx);

    auto u_im_k = roll(u, 1);
    auto u_ip_k = roll(u, -1);

    next_u = (mu * u_im_k) + (2*(1 - mu)*u) + (mu*u_ip_k) - u_i_km;

    if(smoothing){
        next_u[0] = Smooth(t);
    }
    else{
        next_u[0] = Pulse<Ty>(t);
    }

    if(fiexd){
        next_u[N-1] = Fixedcond(um, u);
    }
    else{
        next_u[N-1] = Freecond(um, u);
    }

    return next_u;
}


int main(){
    double time = 0.0;
    std::vector<std::vector<double>>mem_u;
    std::vector<double> um(N, 0);
    mem_u.push_back(um);

    time += Delt;
    std::vector<double> u(N, 0);
    u[0] = Smooth(time);
    mem_u.push_back(u);

    std::vector<double> temp_u(N);

    while(time < 3){
        time += Delt;
        temp_u = step(um, u, time, false, false);
        um = u;
        u = temp_u;
        mem_u.push_back(u);
    }

    int flames = mem_u.size();
    FILE* gp;
    gp = _popen("gnuplot", "w");

    for(int i=0; i<flames; i++){
        if((i+1) % 1 != 0){
            continue;
        }

        fprintf(gp, "unset key\n");
        fprintf(gp, "set terminal png\n");
        fprintf(gp, "set output \'3_3_3_free/%d.png\'\n", i);
        fprintf(gp, "set xrange[%d:%d]\n", 0, 1);
        fprintf(gp, "set yrange[%lf:%lf]\n", -3.2, 3.2);
        fprintf(gp, "set xlabel \"x\"\n");
        fprintf(gp, "set ylabel \"u\"\n");
        fprintf(gp, "set title 'time: %g'\n", i*Delt);
        fprintf(gp, "plot \"-\" with linespoints pt 6 ps 0.5 \n");

        int si = mem_u[0].size();

        for (int j = 0; j < si; j++) {
            fprintf(gp, "%g, %g\n", j * Delx, mem_u[i][j]);
        }

        fprintf(gp, "e\n");
        fprintf(gp, "set output\n");
    }
    fprintf(gp, "set terminal windows\n");
    fflush(gp);
    _pclose(gp);

    return 0;
}