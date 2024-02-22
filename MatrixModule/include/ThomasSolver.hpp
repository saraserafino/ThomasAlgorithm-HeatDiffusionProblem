#ifndef THOMAS_SOLVER_HPP_
#define THOMAS_SOLVER_HPP_

#include <iostream>
#include <vector>

namespace moduleH {

// template class to allow it to work with any matrix type
template<typename MatrixType>
class ThomasSolver {
public:
    // Constructor templated over the matrix type
    ThomasSolver(const MatrixType& matrix);

    std::vector<double> ThomasAlgorithm(std::vector<double>& f);

    ~ThomasSolver() {}; // Default destructor

protected:
    const MatrixType& matrix;
};
} // end of namespace
#include "ThomasSolver.tpl.hpp"
#endif // THOMAS_SOLVER_HPP_