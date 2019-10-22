import sympy
import numpy as np
from numpy import *
from sympy import Matrix, Symbol, lambdify, solve_linear_system, Transpose
from sympy import pprint
from sympy import symbols
import sys
sys.displayhook = pprint
from sympy.matrices import *
from sympy import *
from mpmath import *
from sympy import sympify

#-F(x) = J(x0)s

def limit_diff(M_self, X, i, j, x_zero): #i e j são as coordenadas atuais (Fj(x) e dxi)
        #X = dx, M_self = F(x)
        h = 0.01
        n = 0
        f_xis = M_self
        for x in range(0, X.shape[0]):
            f_xis = f_xis.subs(X[n], x_zero[n])
            n = n + 1

        f_xis_h = M_self
        x_zero_h = x_zero
        x_zero_h[i] = x_zero_h[i] + h

        n = 0
        for x in range(0, X.shape[0]):
            f_xis_h = f_xis_h.subs(X[n], x_zero_h[n])
            n = n + 1

        return ((f_xis_h - f_xis)/h)

def jcb_direct_diff(self, X, x_zero):
        if not isinstance(X, MatrixBase):
            X = self._new(X)
        if self.shape[0] == 1:
            m = self.shape[1]
        elif self.shape[1] == 1:
            m = self.shape[0]
        else:
            raise TypeError("``self`` must be a row or a column matrix")
        if X.shape[0] == 1:
            n = X.shape[1]
        elif X.shape[1] == 1:
            n = X.shape[0]
        else:
            raise TypeError("X must be a row or a column matrix")
        return self._new(m, n, lambda j, i: limit_diff(self[j], X, i, j, x_zero))

def matrix_z(M_self, X, x_zero, matriz_um):
        f_xis = M_self
        n = 0
        for x in range(0, X.shape[0]):
            f_xis = f_xis.subs(X[n], x_zero[n])
            n = n + 1
        f_xis = -1*f_xis
        f_final = matriz_um.col_insert(X.shape[0], f_xis)
        return solve_linear_system(f_final, *X)

def check_diff(maz):
        n = 0
        for x in range(0, matriz_um.shape[0]):
            if abs(maz[n]) > 1e-5:
                return 0
            n = n + 1
        return 1

n_exp_str = input("digite o numero de expressões (variáveis devem ser x1, x2,...xn): ")
n_exp = int(n_exp_str)
M = Matrix(n_exp,1,lambda i,j: i+j)
x_zero = Matrix(n_exp,1,lambda i,j: i+j)

n = 0
for x in range(0, n_exp):
    ss = input("enter a matrix (0): ")
    ff = sympify(ss)
    x_zero[n] = ff
    n = n + 1

n = 0
for x in range(0, n_exp):
    s = input("enter a formula: ")
    f = sympify(s)
    M[n] = f
    n = n + 1


flag = 0
while(1):
    n = 0
    matrix_X = Matrix(n_exp, 1, lambda i,j: i+j)
    for x in range(0, n_exp):
        matrix_X[n] = Symbol("x" + str(n+1))
        n = n + 1

    Jac = Matrix([matrix_X])

    if(flag == 0):
        matriz_um = jcb_direct_diff(M, Jac, x_zero)

    matriz_dois = matrix_z(M, Jac, x_zero, matriz_um)
    res = np.array([list(matriz_dois.values()) for item in matriz_dois.values()])
    res = sympy.Matrix(res)
    maz2_true = res.row(0)
    maz2_true = sympy.Transpose(maz2_true)
    maz = (maz2_true - maz2_true) + maz2_true
    matriz_dois = maz
    var = check_diff(matriz_dois)


    if var > 0:
        break
    else:
        x_zero = (x_zero + matriz_dois)
        flag = flag + 1

print(x_zero)
