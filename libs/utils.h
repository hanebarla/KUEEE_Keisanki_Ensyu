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

#endif