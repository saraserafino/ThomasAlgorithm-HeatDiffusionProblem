#include "../include/Matrix.hpp"
#include "../include/HeatDiffusion.hpp"
#include "../include/ThomasSolver.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// ----------------------------------
// Python interface - trampoline
// ----------------------------------

using namespace moduleH;

class PyMatrix : public Matrix {
public:
    // Inherit the constructors
    using Matrix::Matrix;

    // Trampoline (one for each virtual function).
    virtual double operator()(unsigned int input_row_idx, unsigned int input_col_idx) const override {
        PYBIND11_OVERRIDE_PURE(
            double,         /* Return type. */
            Matrix,         /* Parent class. */
            operator(),     /* Name of function in C++ (must match Python name). */
            input_row_idx,  /* Argument(s). */
            input_col_idx
        );
    }

    virtual void printMatrix() const override {
        PYBIND11_OVERRIDE_PURE(
            void,           /* Return type. */
            Matrix,         /* Parent class. */
            printMatrix    /* Name of function in C++ (must match Python name). */
                            /* Argument(s). */
        );
    }

    std::vector<double> solve(std::vector<double> &f) override {
        PYBIND11_OVERRIDE_PURE(
            std::vector<double>,    /* Return type. */
            Matrix,                 /* Parent class. */
            solve                   /* Name of function in C++ (must match Python name). */
            f                       /* Argument(s). */
        );
    }

    virtual std::vector<double> get_a() const override {
        PYBIND11_OVERRIDE_PURE(
            std::vector<double>,    /* Return type. */
            Matrix,                 /* Parent class. */
            get_a                   /* Name of function in C++ (must match Python name). */
                                    /* Argument(s). */
        );
    }

    virtual std::vector<double> get_b() const override {
        PYBIND11_OVERRIDE_PURE(
            std::vector<double>,    /* Return type. */
            Matrix,                 /* Parent class. */
            get_b                   /* Name of function in C++ (must match Python name). */
                                    /* Argument(s). */
        );
    }

    virtual std::vector<double> get_c() const override {
        PYBIND11_OVERRIDE_PURE(
            std::vector<double>,    /* Return type. */
            Matrix,                 /* Parent class. */
            get_c                   /* Name of function in C++ (must match Python name). */
                                    /* Argument(s). */
        );
    }
};


namespace py = pybind11;

// Wrap as Python module

PYBIND11_MODULE(moduleH, m) {
    m.doc() = "pybind11 moduleH plugin";

    py::class_<Matrix, PyMatrix>(m, "Matrix")
        .def(py::init<int>(),
            py::arg("dimension"))
        .def("get_dimension", &Matrix::get_dimension);

    py::class_<TridiagonalMatrix, Matrix>(m, "TridiagonalMatrix")
        .def(py::init<std::vector<double>, std::vector<double>, std::vector<double>>(),
            py::arg("a"), py::arg("b"), py::arg("c"))
        .def("solveTridiagonal", &TridiagonalMatrix::solve, py::arg("f"))
        .def("readTridiagonalEntry", &TridiagonalMatrix::operator(), py::arg("row"), py::arg("col"))
        .def("printTridiagonalMatrix", &TridiagonalMatrix::printMatrix)
        .def("get_aTridiagonal", &TridiagonalMatrix::get_a)
        .def("get_bTridiagonal", &TridiagonalMatrix::get_b)
        .def("get_cTridiagonal", &TridiagonalMatrix::get_c);
    
        
    py::class_<EigenMatrix, Matrix>(m, "EigenMatrix")
        .def(py::init<std::vector<double>, std::vector<double>, std::vector<double>>(),
            py::arg("a"), py::arg("b"), py::arg("c"))
        .def("solveEigen", &EigenMatrix::solve, py::arg("f"))
        .def("readEigenEntry", &EigenMatrix::operator(), py::arg("row"), py::arg("col"))
        .def("printEigenMatrix", &EigenMatrix::printMatrix)
        .def("get_aEigen", &EigenMatrix::get_a)
        .def("get_bEigen", &EigenMatrix::get_b)
        .def("get_cEigen", &EigenMatrix::get_c);

    py::class_<HeatDiffusion>(m, "HeatDiffusion")
        .def(py::init<int, double, const double, const double>(),
            py::arg("dimension"), py::arg("L"),
            py::arg("boundaryCondition1"), py::arg("boundaryCondition2"))
        .def("evaluate", &HeatDiffusion::evaluate, py::arg("y"), py::arg("parser"))
        .def("setHeatMatrix", &HeatDiffusion::setHeatMatrix)
        .def("computeHeatSource", &HeatDiffusion::computeHeatSource, py::arg("function"))
        .def("validate_solution", &HeatDiffusion::validate_solution, py::arg("solution"), py::arg("exactSol"));

    py::class_<ThomasSolver<TridiagonalMatrix>>(m, "ThomasSolverTridiagonal")
        .def(py::init<TridiagonalMatrix>(), py::arg("matrix"))
        .def("ThomasAlgorithmTridiagonal", &ThomasSolver<TridiagonalMatrix>::ThomasAlgorithm, py::arg("f"));

    py::class_<ThomasSolver<EigenMatrix>>(m, "ThomasSolverEigen")
        .def(py::init<EigenMatrix>(), py::arg("matrix"))
        .def("ThomasAlgorithmEigen", &ThomasSolver<EigenMatrix>::ThomasAlgorithm, py::arg("f"));
}