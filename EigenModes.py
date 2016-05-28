import numpy as np
import matplotlib.pyplot as plt
from numpy import pi, exp
import scipy as sp


# def degeneracy(omega, xyz):
#
#
#     def isgood(n1, n2, n3):
#         x, y, z = xyz
#         return abs(((n1 / float(x)) ** 2 + (n2 / float(y)) ** 2 + (n3 / float(z)) ** 2)/(omega*(x*y*z)**(1./3)/pi)**2 - 1.) <= 1e-2
#
#     x, y, z = xyz
#     a = (x * y * z)**(1./3)
#     n1, n2, n3 = 0, 0, 0
#     g = 0
#
#     while (n1/x) ** 2 + (n2/y) ** 2 + (n3/z) ** 2 < (omega*(x*y*z)**(1./3)/pi)**2:
#
#         if not isgood(n1, n2, n3):
#             pass
#         elif n1 + n2 + n3 >= 1:
#             g = 2 - (n1 == 0 + n2 == 0 + n3 == 0)
#             return g
#         if not isgood(n1+1, n2, n3):
#             pass
#         elif n1 + n2 + n3 >= 1:
#             g += 2 - (n1 == 0 + n2 == 0 + n3 == 0)
#         if not isgood(n1, n2+1, n3):
#             pass
#         elif n1 + n2 + n3 >= 1:
#             g += 2 - (n1 == 0 + n2 == 0 + n3 == 0)
#         if not isgood(n1, n2, n3+1):
#             pass
#         elif n1 + n2 + n3 >= 2:
#             g += 2 - (n1 == 0 + n2 == 0 + n3 == 0)
#         if not isgood(n1+1, n2+1, n3):
#             pass
#         elif n1 + n2 + n3 >= 2:
#             g += 2 - (n1 == 0 + n2 == 0 + n3 == 0)
#         if not isgood(n1+1, n2, n3+1):
#             pass
#         elif n1 + n2 + n3 >= 2:
#             g += 2 - (n1 == 0 + n2 == 0 + n3 == 0)
#         if not isgood(n1, n2+1, n3+1):
#             pass
#         elif n1 + n2 + n3 >= 2:
#             g += 2 - (n1 == 0 + n2 == 0 + n3 == 0)
#
#         n1 += 1
#         n2 += 1
#         n3 += 1
#
#     return g


def degeneracy2(omega, xyz):

    x, y, z = xyz
    a = (x * y * z)**(1./3)
    cutoff_n = int(1e5)
    g = 0

    def isgood(n1, n2, n3):
        x, y, z = xyz
        return abs(np.linalg.norm([n1/x, n2/y, n3/z])/(omega*a/pi) - 1.) < 1

    for n1 in range(0, cutoff_n):
        if n1/float(x) > omega*a/pi:
            break
        for n2 in range(0, cutoff_n):
            if np.linalg.norm([n1/float(x), n2/float(y)]) > omega*a/pi:
                break
            for n3 in range(0, cutoff_n):
                if np.linalg.norm([n1/float(x), n2/float(y), n3/float(z)]) > omega*a/pi:
                    break
                elif isgood(n1, n2, n3) and n1 + n2 + n3 > 1:
                    g += 2 - (n1 == 0 + n2 == 0 + n3 == 0)

    return g


def cuboid(v, alpha, beta):

    z = (alpha * beta * v) ** (1./3)
    x, y = z/alpha, z/beta

    return [x, y, z]


def dos(omega, domega, xyz):

    x, y, z = xyz
    omegas = np.arange(omega, omega + domega, domega)
    gs = np.array([degeneracy2(o, xyz) for o in omegas])
    a = (x * y * z)**(1./3)

    E = np.sum(gs * omegas / (a * (exp(omegas) - 1)))/domega

    return E


omegas = np.linspace(0.1, 10, 500)
domega = 50 * float(omegas[1] - omegas[0])
cube = cuboid(1, 1e-2, 1e-2)
print [degeneracy2(1., cuboid(v, 1, 1)) for v in range(10, 50)]
plt.plot(omegas, [dos(o, domega, cube) for o in omegas], ".")
plt.show()