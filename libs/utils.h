#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <vector>

#include "matrix.h"

template <typename T>
std::ostream& operator<<(std::ostream& stream, const std::vector<T> v) {
    stream << "Vector: [";
    for (const auto& e : v) {
        stream << e << ", ";
    };
    stream << "], Type: " << typeid(v[0]).name();
    return stream;
}

template <typename T>
std::ostream& operator<<(std::ostream& stream, const Matrix<T> v) {
    stream << "Matrix:: row: " << v.row << ", colum: " << v.colum << "\n [";
    for (int i = 0; i < v.row; i++) {
        if (i != 0) {
            stream << "\n  ";
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

#endif