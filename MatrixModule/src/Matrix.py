# ----------------
# Python interface
# ----------------

class Matrix:
    def __init__(self, dimension):
        self.dimension = dimension

# When defining a custom constructor in a derived Python class,
# you must explicitly call the bound C++ constructor using __init__

class TridiagonalMatrix(Matrix):
    def __init__(self, a, b, c):
        Matrix.__init__(self) # Without this, a TypeError is raised
        self.a = a
        self.b = b
        self.c = c

class EigenMatrix(Matrix):
    def __init__(self, a, b, c):
        Matrix.__init__(self) # Without this, a TypeError is raised
        self.a = a
        self.b = b
        self.c = c

class ThomasSolver:
    def __init__(self, matrix):
        self.matrix = matrix

class HeatDiffusion:
    def __init__(self, dimension, L, boundaryCondition1, boundaryCondition2):
        self.dimension = dimension
        self.L = L
        self.boundaryCondition1 = boundaryCondition1
        self.boundaryCondition2 = boundaryCondition2
