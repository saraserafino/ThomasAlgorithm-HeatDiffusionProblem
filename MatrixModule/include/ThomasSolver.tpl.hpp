#include "../include/ThomasSolver.hpp"
#include <iostream>

namespace moduleH {

template<typename MatrixType>
ThomasSolver<MatrixType>::ThomasSolver(const MatrixType& matrix) : matrix(matrix) {};

// The Thomas Algorithm does not assume anything specific about the matrix type,
// it just needs its diagonal vectors and dimension
template<typename MatrixType>
std::vector<double> ThomasSolver<MatrixType>::ThomasAlgorithm(std::vector<double>& f) {
    // Use a non virtual method of Matrix to access its dimension
    unsigned int dimension = matrix.get_dimension();
    // Use virtual methods of MatrixType to access its diagonals
    std::vector<double> a = matrix.get_a();
    std::vector<double> b = matrix.get_b();
    std::vector<double> c = matrix.get_c();

    // Check size compatibility via exception handling
    if (f.size() != dimension) {
        throw std::invalid_argument("Right-hand side vector size f does not match matrix dimensions.");
    }
    // Create the unknown vector
    std::vector<double> u(dimension);
    // Step 1
    for (unsigned int i = 1; i < dimension; ++i) {
        double m = a[i] / b[i - 1];
        b[i] = b[i] - m * c[i - 1];
        f[i] = f[i] - m * f[i - 1];
    }
    // Step 2
    u[dimension - 1] = f[dimension - 1] / b[dimension - 1];
    unsigned int i = dimension - 1;
    while (i >= 1) {
        u[i - 1] = ( f[i - 1] - c[i - 1] * u[i] ) / b[i - 1];
        i -= 1;
    }
    return u;
};
} // end of namespace