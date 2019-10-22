import matplotlib
from matplotlib import pyplot as plt
import sympy
import numpy
from numpy import *
from sympy import Symbol, lambdify, diff
import sys
from sympy.matrices import *
from sympy import *
from mpmath import *
from sympy import sympify

def deriv_v(Va, p, n, R, T, a, b):
    V = Symbol("V")
    func_f = str((p*V**2 + a*n**2)*(V - n*b) - n*R*T*V**2)
    ff = sympify(func_f)
    ff = ff.diff(V)
    return ff.subs(V, Va)

def func_phi(Va, p, n, R, T, a, b):
    V = Symbol("V")
    func_f = str((p*V**2 + a*n**2)*(V - n*b) - n*R*T*V**2)
    ff = sympify(func_f)
    return ff.subs(V, Va)

def calc_Vk(Vk, P, n, R, T, a, b):
    Vkant = 0
    while(abs(Vkant-Vk) > 1e-10):
        print(Vk)
        Vkant = Vk
        Vk = Vk - func_phi(Vk, P, n, R, T, a, b)/(deriv_v(Vk, P, n, R, T, a, b))
    return Vk

T = 734
T = float(T)

p = 90.79694054
R = 0.0820574587
n = 1.0
a = 0.3656
b = 0.00004283

Vk_list = list()
N = 0.1
P = p * N

for x in range(0, 100):
    Vk = n*R*T/P
    Vk = calc_Vk(Vk, P, n, R, T, a, b)
    Vk_list.append(Vk)
    N = N + 0.1
    P = p * N


Z_vetor = list()
x_plot = list()
N = 0.1
P = p * N

for x in range(0,100):
    x_plot.append(N * p)
    N = N + 0.1

N = 0.1
P = p * N
Z = P * Vk_list[0]/(n*R*T)
const = 0

for x in range(0, 100):
    print(const)
    print(P)
    const = Vk_list[x]/(n*R*T)
    Z = const * P
    print(Z)
    Z_vetor.append(Z)
    N = N + 0.1
    P = p * N

plt.plot(x_plot,Z_vetor)
plt.grid(linestyle='--', linewidth=0.5)

plt.show()
