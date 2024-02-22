#ifndef HEAT_DIFFUSION_HPP_
#define HEAT_DIFFUSION_HPP_

#include "Matrix.hpp"
#include "mpParser.h" // for using muParserX
#include <iostream>
#include <vector>

namespace moduleH {

class HeatDiffusion {
public:
    // Constructor
    HeatDiffusion(int dimension, double L, const double boundaryCondition1, const double boundaryCondition2);

    // Function to evaluate in point y
    mup::Value evaluate(double y, mup::ParserX& parser);

    // Set the diagonal values of the heat matrix and the heat source term (i.e. f)
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> setHeatMatrix();
    std::vector<double> computeHeatSource(const std::string &function);

    // Validate the solution against the exact solution
    double validate_solution(const std::vector<double> solution, const std::string &exactSol);

    ~HeatDiffusion() {}; // Default destructor

protected:
int dimension;
std::vector<double> a; // subdiagonal
std::vector<double> b; // diagonal
std::vector<double> c; // superdiagonal
double L; // Upper boundary of domain
const double boundaryCondition1; // Boundary condition for u(0)
const double boundaryCondition2; // Boundary condition for u(L)
std::vector<double> HeatSource; // Heat source term i.e. f(x)
};
} // end of namespace
#endif // HEAT_DIFFUSION_HPP_