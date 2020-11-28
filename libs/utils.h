#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

// ベクターの出力
template <typename T>
std::ostream& operator<<(std::ostream& stream, const std::vector<T> v);

// ベクター同士の四則演算
template <typename T>
std::vector<T> operator-(const std::vector<T>& v, const std::vector<T>& u);

// ベクターとスカラーの四則演算
template <typename T>
std::vector<T> operator-(const std::vector<T>& v, T n);

template <typename T>
std::vector<T> operator/(const std::vector<T>& v, T n);

std::vector<long double> operator/(const std::vector<long double>& v, long double n);

// ベクターのL2ノルム
template <typename T>
double L2norm(const std::vector<T>& x);

// ベクターのL2ノルム(long double)
long double L2norm(const std::vector<long double>& x);

#endif