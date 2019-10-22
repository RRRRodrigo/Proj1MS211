def deriv_v(V, p, n, r, T, a, b):
    return 3*V**2*p + (-2*n*R*T - 2*b*n*p)*V + a*n**2

def func_phi(V, p, n, R, T, a, b):
    return ((p*V**2 + a*n**2)*(V - n*b) - n*R*T*V**2)

T = input("temperatura media: ")
p = input("pressão atmosferica: ")
T = float(T)
p = float(p)
R = 8.3144621
n = 1.0
a = 0.3656
b = 0.00004283

Vo = n*R*T/p
flag = 0
Vk = Vo - func_phi(Vo, p, n, R, T, a, b)/(deriv_v(Vo, p, n, R, T, a, b))
Vkant = 0
while(abs(Vkant-Vk) > 1e-10):
        Vkant = Vk
        Vk = Vo - func_phi(Vk, p, n, R, T, a, b)/(deriv_v(Vk, p, n, R, T, a, b))
        Vk = Vk - func_phi(Vk, p, n, R, T, a, b)/(deriv_v(Vk, p, n, R, T, a, b))
        print(Vk)
