#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <iostream>
#include <vector>
#include <Eigen/Sparse>

namespace moduleH {

// Abstract class
class Matrix {
public:
    // Constructor
    Matrix(int dimension);

    // Method that returns the value of (row,col)
    virtual double operator()(unsigned int row, unsigned int col) const = 0;

    virtual void printMatrix() const = 0;
    
    // Methods for getting the diagonals for each matrix
    virtual std::vector<double> get_a() const = 0;
    virtual std::vector<double> get_b() const = 0;
    virtual std::vector<double> get_c() const = 0;

    // Method for getting the dimension of the matrix
    unsigned int get_dimension() const;

    // Method that solves a linear system given a right hand side f and returns its solution
    virtual std::vector<double> solve(std::vector<double> &f) = 0;

    virtual ~Matrix() {}; // Virtual destructor

protected:
int dimension;
};

// Derived classes from Matrix

class TridiagonalMatrix : public Matrix {
public:
    TridiagonalMatrix(std::vector<double> a, std::vector<double> b, std::vector<double> c);

    double operator()(unsigned int row, unsigned int col) const override;
    void printMatrix() const override;
    std::vector<double> get_a() const override;
    std::vector<double> get_b() const override;
    std::vector<double> get_c() const override;

    // Override it to use the Thomas algorithm
    std::vector<double> solve(std::vector<double> &f) override;

    ~TridiagonalMatrix() {}; // Default destructor

protected:
std::vector<double> a; // subdiagonal
std::vector<double> b; // diagonal
std::vector<double> c; // superdiagonal
};

// Derived class using Eigen's SparseMatrixXd
class EigenMatrix : public Matrix {
public:
    EigenMatrix(std::vector<double> a, std::vector<double> b, std::vector<double> c);

    double operator()(unsigned int row, unsigned int col) const override;
    void printMatrix() const override;
    std::vector<double> get_a() const override;
    std::vector<double> get_b() const override;
    std::vector<double> get_c() const override;

    // Override it to use the Thomas algorithm
    std::vector<double> solve(std::vector<double> &f) override;
    
    ~EigenMatrix() {}; // Default destructor

protected:
Eigen::SparseMatrix<double> sparseMatrix;
std::vector<double> a; // subdiagonal
std::vector<double> b; // diagonal
std::vector<double> c; // superdiagonal
};
} // end of namespace
#endif // MATRIX_HPP_
