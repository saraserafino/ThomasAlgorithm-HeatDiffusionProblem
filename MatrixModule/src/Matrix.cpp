#include "../include/Matrix.hpp"
#include "../include/ThomasSolver.hpp"
#include <iostream>
#include <cmath>
#include <Eigen/Sparse>

namespace moduleH {

// Constructor of the abstract class
Matrix::Matrix(int dimension) : dimension(dimension) {};

unsigned int Matrix::get_dimension() const {
    return dimension;
};

// Derived classes of Matrix

TridiagonalMatrix::TridiagonalMatrix(std::vector<double> a, std::vector<double> b, std::vector<double> c)
    : Matrix(a.size()), a(a), b(b), c(c) {};

// Method for returning the value of the element in (row, col)
double TridiagonalMatrix::operator()(unsigned int row, unsigned int col) const {
    if(row >= this->get_dimension() || col >= this->get_dimension() || row < 0 || col < 0) {
        throw std::out_of_range ("Indexes out of range");
    }
    
    if (row == col) { // diagonal
        return b[row];
    } else if (row == col - 1) {  // superdiagonal
        return c[row];
    } else if (row == col + 1) {  // subdiagonal
        return a[row];
    } else { // not a diagonal value so it's empty
        return 0;
    }
};

void TridiagonalMatrix::printMatrix() const {
    for (unsigned int i = 0; i < dimension; ++i) {
        for (unsigned int j = 0; j < dimension; ++j) {
            // Use the operator() to access values
            std::cout << (*this)(i, j) << "\t";
        }
        // At the end of the row, break line
        std::cout << std::endl;
    }
};

std::vector<double> TridiagonalMatrix::get_a() const {
    return a;
};

std::vector<double> TridiagonalMatrix::get_b() const {
    return b;
};

std::vector<double> TridiagonalMatrix::get_c() const {
    return c;
};

std::vector<double> TridiagonalMatrix::solve(std::vector<double> &f) {
    // Create a ThomasSolver instance for solving it, where *this indicates the matrix
    ThomasSolver<TridiagonalMatrix> thomasSolver(*this);
    // Return its solution
    return thomasSolver.ThomasAlgorithm(f);
};


EigenMatrix::EigenMatrix(std::vector<double> a, std::vector<double> b, std::vector<double> c)
    : Matrix(a.size()), a(a), b(b), c(c)
    {
        Eigen::SparseMatrix<double> sparseMatrix(b.size(), b.size());
        // Insert diagonal values
        // Since a[0] and c[n-1] are undefined, the for-cycle doesn't implement i=0 and i=n
        sparseMatrix.insert(0, 0) = b[0];
        sparseMatrix.insert(0, 1) = c[0];
        for (unsigned int i = 1; i < a.size() - 1; ++i) {
            sparseMatrix.insert(i, i - 1) = a[i]; // subdiagonal
            sparseMatrix.insert(i, i) = b[i]; // diagonal
            sparseMatrix.insert(i, i + 1) = c[i]; // superdiagonal
        }
        sparseMatrix.insert(a.size() - 1, a.size() - 2) = a[a.size() - 1]; // (n,n-1) i.e. a[n-1]
        sparseMatrix.insert(b.size() - 1, b.size() - 1) = b[b.size() - 1]; // (n,n) i.e. b[n-1]
        this->sparseMatrix = sparseMatrix;
    };

double EigenMatrix::operator()(unsigned int row, unsigned int col) const {
    // Implement access to the elements of the sparse matrix using an already existent method of Eigen
    return sparseMatrix.coeff(row, col);
};

void EigenMatrix::printMatrix() const {
    for (int i = 0; i < sparseMatrix.rows(); ++i) {
        for (int j = 0; j < sparseMatrix.cols(); ++j) {
            // Use an already existent method of Eigen for printing the values
            std::cout << sparseMatrix.coeff(i, j) << "\t";
        }
        // At the end of the row, break line
        std::cout << std::endl;
    }
};

std::vector<double> EigenMatrix::get_a() const {
    return a;
};

std::vector<double> EigenMatrix::get_b() const {
    return b;
};

std::vector<double> EigenMatrix::get_c() const {
    return c;
};

std::vector<double> EigenMatrix::solve(std::vector<double> &f) {
    // Create a ThomasSolver instance for solving it, where *this indicates the matrix
    ThomasSolver<EigenMatrix> thomasSolver(*this);
    // Return its solution
    return thomasSolver.ThomasAlgorithm(f);
};
} // end of namespace