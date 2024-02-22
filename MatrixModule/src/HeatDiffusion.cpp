#include "../include/HeatDiffusion.hpp"

#include <numeric> // for std::inner_product

namespace moduleH {

HeatDiffusion::HeatDiffusion(int dimension, double L, const double boundaryCondition1, const double boundaryCondition2)
    : dimension(dimension + 2), L(L),
    boundaryCondition1(boundaryCondition1), boundaryCondition2(boundaryCondition2) {};
// The vectors u (u_0, ..., u_{N+1}) and f (alpha, f(x_1)h^2, ..., f(x_N)h^2, beta) have size "dimension + 2"
// and not just "dimension" because they are computed over N+2 points. It's cleaner changing the dimension
// here in the constuctor since it affects also the matrix A (that must match their dimensions)

// Define the function to evaluate in point y
mup::Value HeatDiffusion::evaluate(double y, mup::ParserX& parser) {
    // Remove any existing variable x to ensure a clean slate
    parser.RemoveVar("x");
    // Create variable x with value a
	mup::Value x(y);
    parser.DefineVar("x", mup::Variable(&x));
    // Evaluate the expression
    return parser.Eval();
};

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> HeatDiffusion::setHeatMatrix() {
    // subdiagonal
    a.resize(dimension, 1);
    a[dimension - 1] = 0;
    // diagonal
    b.resize(dimension, -2);
    b[0] = 1;
    b[dimension - 1] = 1;
    // superdiagonal
    c.resize(dimension, 1);
    c[0] = 0;
    return make_tuple(a, b, c);
};

std::vector<double> HeatDiffusion::computeHeatSource(const std::string &function) {
    // Create the parser instance for the function
    mup::ParserX parser;
    // Set the expression of the function
    parser.SetExpr(function);
    
    // Fill the vector domain (i.e. x) with the equally spaced temperatures
    std::vector<double> domain(dimension);
    // Calculate number of spatial steps
    double h = L / (dimension - 1);
    for (unsigned int i = 0; i < dimension; ++i) {
        domain[i] = i * h;
    }

    // Compute the heat source term i.e. f
    HeatSource.resize(dimension, h * h);
    HeatSource[0] = boundaryCondition1;
    HeatSource[dimension - 1] = boundaryCondition2;
    for (unsigned int i = 1; i < dimension - 1; ++i) {
        // Evaluate the function in each point of the domain
        double ev = evaluate(domain[i], parser);
        HeatSource[i] *= ev;
    }
    return HeatSource;
};

// Validate the solution against the exact solution computing the error ||solution - exactSolution||
// in Euclidean norm to assess the correctness of the implementation
double HeatDiffusion::validate_solution(const std::vector<double> solution, const std::string &exactf) {
    // Create the parser instance for the exact solution exactf
    mup::ParserX parser;
    // Set the expression of the exact solution exactf
    parser.SetExpr(exactf);

    // Calculate number of spatial steps for the domain (i.e. x)
    double h = L / (dimension + 1);
    std::vector<double> exactSol(dimension);
    for (unsigned int i = 0; i < dimension; ++i) {
        // Evaluate the exact solution in each point of the domain (which is i*h equally spaced)
        exactSol[i] = evaluate(i * h, parser);
    }

    // Compute the difference between the two vectors
    std::vector<double> difference(dimension);
    std::transform(solution.begin(), solution.end(), exactSol.begin(), difference.begin(),
                   [](double x, double y) { return x - y; });
    // Compute the Euclidean norm doing the square root of std::inner_product(first1, last1, first2, init)
    // where first1 and first2 coincide because Euclidean norm = sqrt( (x_1)^2 + ... + (x_n)^2 )
    // so you have to multiply the vector difference for itself. 0.0 is the initial value of this sum
    return std::sqrt(std::inner_product(difference.begin(), difference.end(), difference.begin(), 0.0));
};
} // end of namespace