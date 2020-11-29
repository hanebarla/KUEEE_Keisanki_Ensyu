#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

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

// ベクターのL2ノルム
template <typename T>
double L2norm(const std::vector<T>& x){
    double sum = 0;
    for(size_t i=0; i<x.size(); i++){
        sum += double(x[i] * x[i]);
    }

    return sqrtl(sum);
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