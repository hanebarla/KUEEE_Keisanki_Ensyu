#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <vector>
#include <cmath>

#include "matrix.h"

template <typename T>
std::ostream& operator<<(std::ostream& stream, const std::vector<T> v) {
    stream << "[";
    for (const auto& e : v) {
        stream << e << ", ";
    };
    stream << "]";
    return stream;
}

template <typename T>
std::ostream& operator<<(std::ostream& stream, const Matrix<T> v) {
    stream << "[";
    for (int i = 0; i < v.row; i++) {
        if (i != 0) {
            stream << "\n ";
        }
        stream << "[ ";
        for (int j = 0; j < v.colum; j++) {
            stream << v.value[i][j] << ", ";
        }
        stream << "]";
    }
    stream << "]\n";

    return stream;
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

template <typename T>
double L2norm(const std::vector<T>& x){
    double sum = 0;
    for(size_t i=0; i<x.size(); i++){
        sum += double(x[i] * x[i]);
    }

    return sqrtl(sum);
}

long double L2norm(const std::vector<long double>& x){
    double sum = 0;
    for(size_t i=0; i<x.size(); i++){
        sum += x[i] * x[i];
    }

    return sqrtl(sum);
}

#endif