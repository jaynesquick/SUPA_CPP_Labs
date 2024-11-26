#ifndef VECTOR_UTILS_H
#define VECTOR_UTILS_H

#include <iostream>
#include <vector>

// Template for printing vectors of any type! :)
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
    os << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        os << vec[i];
        if (i < vec.size() - 1) {
            os << ", ";
        }
    }
    os << "]";
    return os;
}


// Template for adding two vectors
template <typename T>
std::vector<T> operator+(const std::vector<T>& lhs, const std::vector<T>& rhs) {
    // Check that both vectors are of the same size
    if (lhs.size() != rhs.size()) {
        throw std::invalid_argument("Vector shape mismatch.");
    }

    std::vector<T> result(lhs.size());
    for (size_t i = 0; i < lhs.size(); ++i) {
        result[i] = lhs[i] + rhs[i];  // Sum termz
    }
    return result;
}

#endif // VECTOR_UTILS_H HEADER
