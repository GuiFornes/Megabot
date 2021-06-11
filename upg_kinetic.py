# L'ensemble des distances sont exprimées en mm : segments de patte et élongations des vérins

import numpy as np
from numpy.linalg import inv
from qpsolvers import solve_qp
import kinetic as k

x_A = - k.LEGS[k.FL]['lengths']['ao']
y_A = 0
x_B = k.LEGS[k.FL]['lengths']['bo']
y_B = 0

AB = k.LEGS[k.FL]['lengths']['ao'] + k.LEGS[k.FL]['lengths']['bo']
AC = np.sqrt((AB + k.LEGS[k.FL]['lengths']['bcx'])**2 + k.LEGS[k.FL]['lengths']['bcy']**2)
BC = np.sqrt(k.LEGS[k.FL]['lengths']['bcx']**2 + k.LEGS[k.FL]['lengths']['bcy']**2)
AE = k.LEGS[k.FL]['lengths']['ae']
AD = AE - k.LEGS[k.FL]['lengths']['de']
BF = k.LEGS[k.FL]['lengths']['bf']
FH = k.LEGS[k.FL]['lengths']['fh']
BH = BF - FH
FG = k.LEGS[k.FL]['lengths']['fg']
EF = k.LEGS[k.FL]['lengths']['ef']
EG = EF + FG
GI = k.LEGS[k.FL]['lengths']['gi']
GJ = k.LEGS[k.FL]['lengths']['gj']

KO = k.LEGS[k.FL]['lengths']['yaw_c']
LM = k.LEGS[k.FL]['lengths']['yaw_b']
MO = k.LEGS[k.FL]['lengths']['yaw_a']
LO = np.sqrt(LM**2 + MO**2)

ori = np.array([[1, 1, -1, -1], [1, -1, -1, 1]]) #[[oritentation selon x][orientation selon y]]

V = 460, 565, 500


def solve_indirect(x, y, z, x0, y0, z0, v1, v2, v3, leg_id, pts):
  X, Y, calpha = d3_to_d2(x, y, z)
  Xt = np.array([X/1000, Y/1000])
  #dv3 = (- np.sqrt(1-calpha**2) / v3 ) * np.arctan(x/y)-np.arctan(x0-y0)
  new_v3 = cos_angle_to_v3(calpha)

  X, Y, calpha = d3_to_d2(x0, y0, z0)
  X0 = np.array([X/1000, Y/1000])

  J = gen_jacob(pts, v1/1000, v2/1000)

  P = 2 * J.T @ J 
  q = J.T @ (X0 - Xt)
  G = np.array([
    [1, 1],
    [1, 1],
  ])
  h = np.array([0.01, 0.01])
  lb = np.array([450.0 - v1, 450.0 - v2])
  ub = np.array([650.0 - v1, 650.0 - v2])
  dV = solve_qp(P, q, G, h, lb=lb, ub=ub)

  return v1+dV[0]*1000, v2+dV[1]*1000, new_v3


  
def d3_to_d2(x, y, z):
  X = np.sqrt(x**2 + y**2)
  Y = z
  calpha = y/X
  return X, Y, calpha

def d2_to_d3(X, Y, calpha):
  x = X * np.sqrt(1-calpha**2)
  y = X * calpha
  z = Y
  return x, y, z



'''
Retourne la Jacobienne correspondant au modèle cinématique indirect
Prend en argument la position des points de la patte et l'élongation des verrins en m
'''
def gen_jacob(pts, v1, v2):
  x_A, y_A = pts['A']
  x_B, y_B = pts['B']
  x_C, y_C = pts['C']
  x_D, y_D = pts['D']
  x_E, y_E = pts['E']
  x_F, y_F = pts['F']
  x_G, y_G = pts['G']
  x_H, y_H = pts['H']
  x_I, y_I = pts['I']

  A = np.array([
    [2*(x_D - x_A), 2*(y_D - y_A)],
    [2*(x_D - x_C), 2*(y_D - y_C)]])
  B = np.array([0, 2*v1])
  M_D = inv(A) @ B

  M_E = (AE/AD) * M_D

  A = np.array([
    [2*(x_F - x_E), 2*(y_F - y_E)], 
    [2*(x_F - x_B), 2*(y_F - y_B)]])
  B = np.array([
    [2*(x_F - x_E), 2*(y_F - y_E)], 
    [0, 0]])
  M_F = (inv(A) @ B) @ M_E

  M_G = (EG/EF) * M_F

  M_H = (BH/BF) * M_F

  A = np.array([
    [2*(x_I - x_G), 2*(y_I - y_G)],
    [2*(x_I - x_H), 2*(y_I - y_H)]])
  B = np.array([
    [2*(x_I - x_G), 2*(y_I - y_G)],
    [0, 0]])
  C = np.array([
    [0, 0],
    [2*(x_I - x_H), 2*(y_I - y_H)]])
  D = np.array([0, 2*v2])
  V1 = inv(A) @ (B @ M_G + C @ M_H)
  V2 = inv(A) @ D 
  M_I = np.array([
    [V1[0], V2[0]],
    [V1[1], V2[1]]
  ])
  Jacob = inv((GJ/GI) * M_I)

  return Jacob



'''
Fonction auxiliaire de move_xyz
Retourne l'élongation de v3 en fonction de l'angle de la patte au chassis
'''
def cos_angle_to_v3(cangle):
  #print(angle)
  #print(beta)
  return np.sqrt(LO**2 + KO**2 - 2*LO*KO * (cangle * MO/LO + np.sqrt(1-cangle**2)*np.sqrt(1-(MO/LO)**2)))

'''
Fonction auxiliaire de move_xyz
Retourne l'angle de la patte au chassis en fonction de l'élongation de v3
'''
def v3_to_cos_angle(v3):
  #print(v3)
  #print(beta)
  return (KO**2 + LO**2 - v3**2)/(2*KO*LO) * MO/LO - np.sqrt(1 - ((KO**2 + LO**2 - v3**2)/(2*KO*LO))**2) * np.sqrt(1-(MO/LO)**2)



'''
Retourne la liste des élongations successives des verrins permettant de placer le bout de la patte en (x,y,z)
'''
def move_xyz(x, y, z, v1, v2, v3, dstep, p, eps, leg_id):
  L = []
  c = 0

  pts = k.get_leg_points_V1_V2(v1/1000, v2/1000, k.LEGS[k.FL]['lengths'])
  X, Y = pts['J'][0]*1000, pts['J'][1]*1000
  x0, y0, z0 = d2_to_d3(X, Y, v3_to_cos_angle(v3))
  dist = distance(x - x0, y - y0, z - z0)
  while dist > eps and c < 300:
    c += 1
    U = np.array([(x - x0), (y - y0), (z - z0)])
    U = dstep / 100 * U # / np.linalg.norm(U)**p
    dx, dy, dz = U[0], U[1], U[2]

    v1, v2, v3 = solve_indirect(x0+dx, y0+dy, z0+dz, x0, y0, z0, v1, v2, v3, leg_id, pts)
    L.append((v1, v2, v3))
    
    pts = k.get_leg_points_V1_V2(v1/1000, v2/1000, k.LEGS[k.FL]['lengths'])
    X, Y = pts['J'][0]*1000, pts['J'][1]*1000
    x0, y0, z0 = d2_to_d3(X, Y, v3_to_cos_angle(v3))
    dist = distance(x - x0, y - y0, z - z0)

  return L



'''
Retourne les positions x, y, z du bout de la patte en fonctions des v1, v2, v3
'''
def direct_xyz(v1 ,v2, v3, leg_id):
  X, Y = k.get_leg_points_V1_V2(v1/1000, v2/1000, k.LEGS[k.FL]['lengths'])['J']
  #print (X, Y)
  calpha = v3_to_cos_angle(v3)
  #alpha = k.v3_to_delta_angle(v3, k.LEGS[k.FL]['lengths'])
  x = ori[leg_id][0] * X * np.sqrt(1-calpha**2) * 1000
  y = ori[leg_id][1] * X * calpha * 1000
  z = Y * 1000
  return x, y, z



'''
Calcule la distance euclidienne dans l'espace en 3 dimensions (ou 2D selon le nombre de coordonnées passées en paramètre)
'''
def distance(x, y, z=0):
  return np.sqrt(x**2+y**2+z**2)



################################ DEPRECATED ################################

def distance_euclidienne(xi, yi, xj, yj):
  return np.sqrt((xj - xi)**2 + (yj - yi)**2)

def al_kashi_longueur(a,b,alpha):
  return np.sqrt(a**2 + b**2 + 2*a*b*np.cos(alpha))

def al_kashi_angle(a, b, c):
  return np.arccos((a**2 + b**2 - c**2) / (2*a*b))

alpha = al_kashi_angle(AB, AC, BC)

def direct_v1(v1, v2):
  theta1 = alpha + al_kashi_angle(AD, AC, v1)
  x_E = x_A + AE * np.cos(theta1)
  y_E = y_A + AE * np.sin(theta1)
  EB = distance_euclidienne(x_E, y_E, x_B, y_B)

  theta2 = al_kashi_angle(EF, EB, BF)
  theta3 = al_kashi_angle(AE, EB, AB)

  beta = theta2 + theta3 - (np.pi - theta1)
  x_F = x_E + EF * np.cos(beta)
  y_F = y_E + EF * np.sin(beta)
  x_G = x_E + EG * np.cos(beta)
  y_G = y_E + EG * np.sin(beta)

  x_H = (FH * x_B + BH * x_F) / (FH + BH)
  y_H = (FH * y_B + BH * y_F) / (FH + BH)

  GH = distance_euclidienne(x_G, x_H, y_G, y_H)
  
  theta4 = al_kashi_angle(FG, GH, FH)
  theta5 = al_kashi_angle(GI, GH, v2)
  theta6 = np.pi - (theta4 + theta5 + beta)
  
  return x_G + GJ * np.cos(theta6), y_G - GJ * np.sin(theta6)

'''
Fonction auxiliaire de move_xyz
Retourne les valeurs de dx, dy et dz en fonction des angles theta1 et theta2 et de la distance d'un pas
'''
def deltas(theta1, theta2, dstep):
  deltaX = dstep * np.cos(theta2)
  dz = dstep * np.sin(theta2)
  dx = deltaX * np.cos(theta1)
  dy = deltaX * np.sin(theta1)
  return dx, dy, dz

############################################################################
if __name__ == "__main__":
    import doctest
    doctest.testmod()