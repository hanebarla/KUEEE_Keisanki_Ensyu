#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <numeric>

double machine_epsilon(){
    double e = 1.0;
    double e_m = 1.0;
    while(true){
        e_m /= 2.0;
        if((1.0 + e_m) == 1.0){
            break;
        }
        e = e_m;
    }

    return e;
}

// ベクターの出力
template <typename T>
std::ostream& operator<<(std::ostream& stream, const std::vector<T> v) {
    stream << "[";
    for (const auto& e : v) {
        stream << e << ", ";
    };
    stream << "]";
    return stream;
}

// ベクター同士の四則演算
template <typename T>
std::vector<T> operator+(const std::vector<T>& v, const std::vector<T>& u){
    int vsize = v.size();
    int usize = u.size();

    if (vsize != usize) throw std::range_error("Don't match the shape");

    std::vector<T> Ans(vsize, 0);
    for(int i=0; i<vsize; i++){
        Ans[i] = v[i] + u[i];
    }

    return Ans;
}

template <typename T>
std::vector<T> operator*(const std::vector<T>& v, const std::vector<T>& u){
    int vsize = v.size();
    int usize = u.size();

    if (vsize != usize) throw std::range_error("Don't match the shape");

    std::vector<T> Ans(vsize, 0);
    for(int i=0; i<vsize; i++){
        Ans[i] = v[i] * u[i];
    }

    return Ans;
}

template <typename T>
std::vector<T> operator-(const std::vector<T>& v, const std::vector<T>& u){
    int vsize = v.size();
    int usize = u.size();

    if (vsize != usize) throw std::range_error("Don't match the shape");

    std::vector<T> Ans(vsize, 0);
    for(int i=0; i<vsize; i++){
        Ans[i] = v[i] - u[i];
    }

    return Ans;
}

// ベクターとスカラーの四則演算
template <typename T>
std::vector<T> operator+(const std::vector<T>& v, T n){
    int vsize = v.size();
    std::vector<T> v2(vsize);
    for(int i=0; i<vsize; i++){
        v2[i] = v[i] + n;
    }
    return v2;
}

template <typename T>
std::vector<T> operator+(T n, const std::vector<T>& v){
    return v + n;
}

template <typename T>
std::vector<T> operator*(const std::vector<T>& v, T n){
    int vsize = v.size();
    std::vector<T> v2(vsize);
    for(int i=0; i<vsize; i++){
        v2[i] = v[i] * n;
    }
    return v2;
}

template <typename T>
std::vector<T> operator*(T n, const std::vector<T>& v){
    return v * n;
}

template <typename T>
std::vector<T> operator-(const std::vector<T>& v, T n){
    int vsize = v.size();
    std::vector<T> v2(vsize);
    for(int i=0; i<vsize; i++){
        v2[i] = v[i] - n;
    }
    return v2;
}

template <typename T>
std::vector<T> operator/(const std::vector<T>& v, T n){
    int vsize = v.size();
    std::vector<T> v2(vsize);
    for(int i=0; i<vsize; i++){
        v2[i] = v[i] / n;
    }
    return v2;
}

std::vector<long double> operator/(const std::vector<long double>& v, long double n){
    int vsize = v.size();
    std::vector<long double> v2(vsize);
    for(int i=0; i<vsize; i++){
        v2[i] = v[i] / n;
    }
    return v2;
}

// ベクターの内積
template <typename T>
T dot(const std::vector<T>& x, const std::vector<T>& y){
    T sum = 0;
    for (size_t i = 0; i < x.size(); i++){
        sum += T(x[i] * y[i]);
    }
    return sum;
}

// ベクター同士の外積
template <typename T>
std::vector<T> cross(const std::vector<T>& x, const std::vector<T>& y){
    if (int(x.size()) != 3) throw std::range_error("Not 3-Dimention");

    std::vector<T> ans(3, 0);
    ans[0] = x[1] * y [2] - x[2] * y[1];
    ans[1] = x[2] * y[0] - x[0] * y[2];
    ans[2] = x[0] * y[1] - x[1] * y[0];

    return ans;
}

// ベクターのL2ノルム
template <typename T>
T L2norm(const std::vector<T>& x){
    double sum = 0;
    for(size_t i=0; i<x.size(); i++){
        sum += double(x[i] * x[i]);
    }

    return T(sqrtl(sum));
}

// ベクターのL2ノルム(long double)
long double L2norm(const std::vector<long double>& x){
    double sum = 0;
    for(size_t i=0; i<x.size(); i++){
        sum += x[i] * x[i];
    }

    return sqrtl(sum);
}

#endif