from sympy import *
from upg_kinetic import *
from upg_jacobian import *

# Position dans la base flottante / dans le référentiel absolu
l, m, n, xx0, yy0, zz0, xx1, yy1, zz1, xx2, yy2, zz2, xx3, yy3, zz3, xx4, yy4, zz4 = symbols('l m n xx0 yy0 zz0 xx1 yy1 zz1 xx2 yy2 zz2 xx3 yy3 zz3 xx4 yy4 zz4')
# Positions dans le référentiel du robot
x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4 = symbols('x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4')

init_printing(use_unicode=True)

# Matrices de rotation
R_l = Matrix([[cos(l), -sin(l), 0],
       [sin(l), cos(l), 0],
       [0, 0, 1]])

R_m = Matrix([[1, 0, 0],
       [0, cos(m), -sin(m)],
       [0, sin(m), cos(m)]])

R_n = Matrix([[cos(n), -sin(n), 0],
       [sin(n), cos(n), 0],
       [0, 0, 1]])

R = R_l * R_m * R_n

PP = [Matrix([xx1, yy1, zz1]), Matrix([xx2, yy2, zz2]), Matrix([xx3, yy3, zz3]), Matrix([xx4, yy4, zz4])]
P = [Matrix([x1, y1, z1]), Matrix([x2, y2, z2]), Matrix([x3, y3, z3]), Matrix([x4, y4, z4])]
O = Matrix([xx0, yy0, zz0])
Omega = Matrix([l, m, n])

HG = []
for i in range(4):
       expr1 = PP[i] - O - R * P[i]
       J11 = diff(expr1, PP[i][0]).row_join(diff(expr1, PP[i][1])).row_join(diff(expr1, PP[i][2]))
       J12 = diff(expr1, P[i][0]).row_join(diff(expr1, P[i][1])).row_join(diff(expr1, P[i][2]))
       J13 = diff(expr1, O[0]).row_join(diff(expr1, O[1])).row_join(diff(expr1, O[2]))
       J14 = diff(expr1, Omega[0]).row_join(diff(expr1, Omega[1])).row_join(diff(expr1, Omega[2]))

       M1 = (R * Matrix([1, 0, 0])).row_join(R * Matrix([0, 1, 0])).row_join(PP[i])
       M2 = (R * Matrix([1, 0, 0])).row_join(PP[i]).row_join(R * Matrix([0, 1, 0]))
       M3 = PP[i].row_join(R * Matrix([1, 0, 0])).row_join(R * Matrix([0, 1, 0]))
       expr2 = Matrix([M1.det(), M2.det(), M3.det()])
       J21 = diff(expr2, PP[i][0]).row_join(diff(expr2, PP[i][1])).row_join(diff(expr2, PP[i][2]))
       J22 = diff(expr2, P[i][0]).row_join(diff(expr2, P[i][1])).row_join(diff(expr2, P[i][2]))
       J23 = diff(expr2, O[0]).row_join(diff(expr2, O[1])).row_join(diff(expr2, O[2]))
       J24 = diff(expr2, Omega[0]).row_join(diff(expr2, Omega[1])).row_join(diff(expr2, Omega[2]))

       G = J11.row_join(J12).col_join(J21.row_join(J22))
       H = J13.row_join(J14).col_join(J23.row_join(J24))
       HG.append((H, G))