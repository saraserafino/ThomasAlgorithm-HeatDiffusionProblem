# Import the module created with pybind11
import moduleH

from tabulate import tabulate # for a cute table of results
import time # for the wrapper execution_time
import numpy as np
import matplotlib.pyplot as plt # for plot
import os
# for validating heat solution with NumPy and Scipy
from numpy.linalg import solve
from scipy.linalg import solve_banded

# Decorator for computing the execution time
def execution_time(func):
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        executionTime = time.time() - start
        return result, executionTime
    return wrapper

# --------------------------------------------------
# definition of functions for heat diffusion problem
# --------------------------------------------------

# Solve the heat diffusion problem using the Tridiagonal Matrix implemented in C++
@execution_time
def solveHeatTridiagonal_cpp(HDproblem, function):
    # Set the diagonal values of the Heat matrix
    a, b, c = HDproblem.setHeatMatrix()
    # Create a Tridiagonal Matrix with them
    matrix = moduleH.TridiagonalMatrix(a, b, c)

    # In case you want to print the matrix and the diagonal vectors to check them
    #matrix.printTridiagonalMatrix()
    #print(f"Subdiagonal a: {matrix.get_aTridiagonal()}")
    #print(f"Diagonal b: {matrix.get_bTridiagonal()}")
    #print(f"Superdiagonal c: {matrix.get_cTridiagonal()}")

    # Compute the heat source term f starting from f(x)
    f = HDproblem.computeHeatSource(function)
    # Return the solution
    return matrix.solveTridiagonal(f)

# Solve the heat diffusion problem using the Eigen Matrix implemented in C++
@execution_time
def solveHeatEigen_cpp(HDproblem, function):
    # Set the diagonal values of the Heat matrix
    a, b, c = HDproblem.setHeatMatrix()
    # Create a Eigen Matrix with them
    matrix = moduleH.EigenMatrix(a, b, c)

    # In case you want to print the matrix and the diagonal vectors to check them
    #matrix.printEigenMatrix()
    #print(f"Subdiagonal a: {matrix.get_aEigen()}")
    #print(f"Diagonal b: {matrix.get_bEigen()}")
    #print(f"Superdiagonal c: {matrix.get_cEigen()}")

    # Compute the heat source term f starting from f(x)
    f = HDproblem.computeHeatSource(function)
    # Return the solution
    return matrix.solveEigen(f)

# Solve the heat diffusion problem using NumPy
@execution_time
def solveHeat_numpy(dimension, L, boundaryCondition1, boundaryCondition2, f):
    # Define matrix A
    A = np.diag(- 2 * np.ones(dimension + 2)) + np.diag(np.ones(dimension + 1), 1) + np.diag(np.ones(dimension + 1), -1)
    # subdiagonal
    A[-1, -2] = 0 #a[-1]=0
    # diagonal
    A[0, 0] = 1
    A[-1, -1] = 1
    # superdiagonal
    A[0, 1] = 0 # c[0]=0

    # Calculate number of spatial steps
    h = L / (dimension + 1)
    # Compute points of the domain
    x = np.arange(0, dimension + 2) * h

    # Construct heat source term f
    f_vec = np.zeros(dimension + 2)
    f_vec[0] = boundaryCondition1
    f_vec[1:-1] = f(x[1:-1]) * h**2
    f_vec[-1] = boundaryCondition2

    # Solve the system using NumPy
    return np.linalg.solve(A, f_vec)

# Solve the heat diffusion problem using SciPy
@execution_time
def solveHeat_scipy(dimension, L, boundaryCondition1, boundaryCondition2, f):
    # Calculate number of spatial steps
    h = L / (dimension + 1)
    # Compute points of the domain
    x = np.arange(0, dimension + 2) * h

    # Construct heat source term f
    f_vec = np.zeros(dimension + 2)
    f_vec[0] = boundaryCondition1
    f_vec[1:-1] = f(x[1:-1]) * h**2
    f_vec[-1] = boundaryCondition2

    # Define matrix A in banded format (must be defined from subdiagonal to superdiagonal)
    diagonals = np.array([np.ones(dimension + 2), -2 * np.ones(dimension + 2), np.ones(dimension + 2)])

    # subdiagonal
    diagonals[-1][-2] = 0 #a[-1]=0
    # diagonal
    diagonals[0][0] = 1
    diagonals[-1][-1] = 1
    # superdiagonal
    diagonals[0][1] = 0 # c[0]=0

    # Solve the system using SciPy
    return solve_banded((1, 1), diagonals, f_vec)

# Validate Python solution against the exact solution
def validate_solution_py(solution, exactSol):
    # Compute the difference between the two vectors
    difference = np.array(solution) - np.array(exactSol)
    # Compute the Euclidean norm = sqrt( (x_1)^2 + ... + (x_n)^2 )
    return np.sqrt(np.sum(difference**2))

# Create a folder in which to store plot images
if not os.path.exists("images"):
    os.mkdir("images")

# -------------------------------
# main for Heat Diffusion problem
# -------------------------------

# Initial data of the heat diffusion problem
dimension = 10
L = 1
boundaryCondition = 0 # both alpha and beta are 0

# Create an instances of the heat diffusion problem
HDproblem1 = moduleH.HeatDiffusion(dimension, L, boundaryCondition, boundaryCondition)

# Define the heat source term and the exact analytical solution (muParseX is used in the C++ code)
f_cpp = "-sin(pi*x)"
exactf_cpp = "(sin(pi*x))/(pi*pi)"

# Compute solution with the Tridiagonal Matrix implemented in C++
solutionTridiag_cpp, timesolTridiag_cpp = solveHeatTridiagonal_cpp(HDproblem1, f_cpp)
# Validate it against the exact solution
errorTridiag_cpp = HDproblem1.validate_solution(solutionTridiag_cpp, exactf_cpp)

# Compute solution with the Eigen Matrix implemented in C++
HDproblem2 = moduleH.HeatDiffusion(dimension, L, boundaryCondition, boundaryCondition)
solutionEigen_cpp, timesolEigen_cpp = solveHeatEigen_cpp(HDproblem2, f_cpp)
# Validate it against the exact solution
errorEigen_cpp = HDproblem2.validate_solution(solutionEigen_cpp, exactf_cpp)


# Validate results against the solve function from numpy.linalg

# Define the heat source term and the exact analytical solution
f_py = lambda x: -np.sin(np.pi * x)

exactf_py = lambda x: np.sin(np.pi * x) / np.pi**2

# Compute solution with NumPy and SciPy
solution_numpy, timesol_numpy = solveHeat_numpy(dimension, L, boundaryCondition, boundaryCondition, f_py)
solution_scipy, timesol_scipy = solveHeat_scipy(dimension, L, boundaryCondition, boundaryCondition, f_py)

# Compute exact solution

# Calculate number of spatial steps for the domain (i.e. x)
h = L / (dimension + 1)
# Evaluate the exact solution in each point of the domain (which is i*h equally spaced)
exact_solution = [exactf_py(i * h) for i in range(dimension + 2)]

# Validate the solutions with the exact one
error_numpy = validate_solution_py(solution_numpy, exact_solution)
error_scipy = validate_solution_py(solution_scipy, exact_solution)

# Prepare the results to be printed in a more uniform way with tabulate
header = ["", "Exact", "Tridiagonal", "Eigen", "NumPy", "SciPy"]
data = [
        ["Solution[0]", exact_solution[0], solutionTridiag_cpp[0], solutionEigen_cpp[0], solution_numpy[0], solution_scipy[0]],
        ["Solution[1]", exact_solution[1], solutionTridiag_cpp[1], solutionEigen_cpp[1], solution_numpy[1], solution_scipy[1]],
        ["Solution[2]", exact_solution[2], solutionTridiag_cpp[2], solutionEigen_cpp[2], solution_numpy[2], solution_scipy[2]],
        ["Solution[3]", exact_solution[3], solutionTridiag_cpp[3], solutionEigen_cpp[3], solution_numpy[3], solution_scipy[3]],
        ["Solution[4]", exact_solution[4], solutionTridiag_cpp[4], solutionEigen_cpp[4], solution_numpy[4], solution_scipy[4]],
        ["Solution[n]", round(exact_solution[-1], 3), solutionTridiag_cpp[-1], solutionEigen_cpp[-1], solution_numpy[-1], solution_scipy[-1]],
        ["Error against exact solution", None, round(errorTridiag_cpp, 5), round(errorEigen_cpp, 5), round(error_numpy, 5), round(error_scipy, 5)],
        ["Execution time (s)", None, round(timesolTridiag_cpp, 5), round(timesolEigen_cpp, 5), round(timesol_numpy, 5), round(timesol_scipy, 5)]
        ]
print(f"Solution of the heat diffusion problem with heat source term {f_cpp} and exact solution {exactf_cpp}:")
print(tabulate(data, header, tablefmt = "fancy_grid"))


# Plot the numerical and exact solutions vs x for visual comparison

# Generate x values
x_values = [i * h for i in range(dimension + 2)]

# Plot exact solution
plt.plot(x_values, exact_solution, label = 'Exact Solution', color = 'blue')

# Plot solutions computed with each method
plt.plot(x_values, solutionTridiag_cpp, label = 'Tridiagonal Solution', linestyle = '--', marker = 'o', color = 'red')
plt.plot(x_values, solutionEigen_cpp, label = 'Eigen Solution', linestyle = ':', marker = 's', color = 'cyan')
plt.plot(x_values, solution_numpy, label = 'NumPy Solution', linestyle = ':', marker = 's', color = 'green')
plt.plot(x_values, solution_scipy, label = 'SciPy Solution', linestyle = '--', marker = 'o', color = 'orange')

plt.xlabel('x')
plt.ylabel('Temperature')
plt.title('Numerical Solutions VS Exact Solution')
plt.legend()
plt.grid(True)
# Save the plot
plt.savefig('images/NumericalVSExactSolution.png')
plt.show()

print("Thomas Algorithm not only solves the heat diffusion problem, but systems of equations with a tridiagonal matrix")

# ------------------------------------------------------
# definition of functions for tridiagonal linear systems
# ------------------------------------------------------

# Solve a tridiagonal linear system using the Tridiagonal Matrix implemented in C++
@execution_time
def solveTridiagonal_cpp(a, b, c, f):
    # Create a Tridiagonal Matrix with them
    matrix = moduleH.TridiagonalMatrix(a, b, c)

    # Return the solution
    return matrix.solveTridiagonal(f)

# Solve a tridiagonal linear system using the Eigen Matrix implemented in C++
@execution_time
def solveEigen_cpp(a, b, c, f):
    # Create a Eigen Matrix with them
    matrix = moduleH.EigenMatrix(a, b, c)

    # Return the solution
    return matrix.solveEigen(f)

# Solve a tridiagonal linear system using NumPy
@execution_time
def solve_numpy(a, b, c, f):
    # Construct the tridiagonal matrix
    A = np.zeros((len(b), len(b)))
    # Fill it with the diagonal values
    np.fill_diagonal(A[1:, :], a[1:]) # a[0] is undefined
    np.fill_diagonal(A, b)
    np.fill_diagonal(A[:, 1:], c)

    # Solve the system using NumPy
    return np.linalg.solve(A, f)

# -----------------------------------
# main for tridiagonal linear systems
# -----------------------------------

print("Solution of the system with")
a = [0, 5, 4, 1, 8]
b = [3, 2, 5, -1, 1]
c = [-1, 2, 1, 9, 0]
print("Matrix A: ")
trimatrix = moduleH.TridiagonalMatrix(a,b,c)
trimatrix.printTridiagonalMatrix()
# Of course printing the other matrix is the same
#eimatrix = moduleH.EigenMatrix(a,b,c)
#eimatrix.printEigenMatrix()

# If you want to check an entry you can use the overload of the C++ operator()
# which with Pybind is called readMatrixEntry with Matrix the type of matrix
#print(trimatrix.readTridiagonalEntry(1,0))
#print(eimatrix.readEigenEntry(1,0))

f = [-29, -27, 21, 5, 74]
print(f"Right-hand side vector f: {f}")

solutionTri, timeTri = solveTridiagonal_cpp(a, b, c, f)
solutionEi, timeEi = solveEigen_cpp(a, b, c, f)
solutionNumpy, timeNumpy = solve_numpy(a, b, c, f)
exactsolution = [-7, 8, -4, 9, 2]

# Prepare the results to be printed in a more uniform way with tabulate
header = ["", "Exact", "Tridiagonal", "Eigen", "NumPy"]
data = [
        ["Solution[0]", exactsolution[0], solutionTri[0], solutionEi[0], solutionNumpy[0]],
        ["Solution[1]", exactsolution[1], solutionTri[1], solutionEi[1], solutionNumpy[1]],
        ["Solution[2]", exactsolution[2], solutionTri[2], solutionEi[2], solutionNumpy[2]],
        ["Solution[3]", exactsolution[3], solutionTri[3], solutionEi[3], solutionNumpy[3]],
        ["Solution[4]", exactsolution[-1], solutionTri[-1], solutionEi[-1], solutionNumpy[-1]],
        ["Execution time (s)", None, round(timeTri, 5), round(timeEi, 5), round(timeNumpy, 5)]
        ]
print(tabulate(data, header, tablefmt = "fancy_grid"))