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

def solve_indirect_cyl(x, y, z, x0, y0, z0, v1, v2, v3, leg_id, pts):
  X, Z, calpha = d3_to_d2(x, y, z)
  Xt = np.array([X/1000, Z/1000])
  new_v3 = cos_angle_to_v3(calpha)

  X, Z, calpha = d3_to_d2(x0, y0, z0)
  X0 = np.array([X/1000, Z/1000])

  J = gen_jacob_plan(pts, v1/1000, v2/1000)

  P = 2 * J.T @ J  
  q = J.T @ (X0 - Xt)
  lb = np.array([(450.0 - v1)/1000, (450.0 - v2)/1000])
  ub = np.array([(650.0 - v1)/1000, (650.0 - v2)/1000])
  dV = solve_qp(P, q, lb=lb, ub=ub)

  return v1+dV[0]*1000, v2+dV[1]*1000, new_v3



def solve_indirect_cart(x, y, z, x0, y0, z0, v1, v2, v3, leg_id, pts):
  Xt = np.array([x, y, z]) / 1000
  X0 = np.array([x0, y0, z0]) / 1000

  J = gen_matrix_leg(v1, v2, v3, np.arctan(x0/y0), k.FL) # MARCHE PAS POUR AUTRE QUE k.FL

  P = 2 * J.T @ J  
  q = J.T @ (X0 - Xt)
  lb = np.array([450.0 - v1, 450.0 - v2, 450.0 - v3]) / 1000
  ub = np.array([650.0 - v1, 650.0 - v2, 650.0 - v3]) / 1000
  dV = solve_qp(P, q, lb=lb, ub=ub)

  return v1+dV[0]*1000, v2+dV[1]*1000, v3+dV[2]*1000
  


def d3_to_d2(x, y, z):
  X = np.sqrt(x**2 + y**2)
  Z = z
  calpha = y/X
  return X, Z, calpha

def d2_to_d3(X, Z, calpha):
  x = X * np.sqrt(1-calpha**2)
  y = X * calpha
  z = Z
  return x, y, z

def distance_3_points(A, B, C):
  '''
  Retourne la matrice des distances de A à B et à C
  '''
  return 2 * np.array([[A[0] - B[0], A[1] - B[1]], [A[0] - C[0], A[1] - C[1]]])

def gen_jacob_plan(pts, v1, v2):
  '''
  Retourne la Jacobienne correspondant au modèle cinématique indirect dans le plan de la patte 
  Prend en argument la position des points de la patte et l'élongation des verrins en m

  >>> gen_jacob_plan(k.get_leg_points_V1_V2(0.495, 0.585, k.LEGS[k.FL]['lengths']), 0.495, 0.585) @ np.array([0, 0])
  array([0., 0.])

  >>> gen_jacob_plan(k.get_leg_points_V1_V2(0.495, 0.585, k.LEGS[k.FL]['lengths']), 0.495, 0.585) @ np.array([1, 0])
  array([1, 0])
  
  >>> gen_jacob_plan(k.get_leg_points_V1_V2(0.495, 0.585, k.LEGS[k.FL]['lengths']), 0.495, 0.585) @ np.array([1, 0])

  '''
  x_E, y_E = pts['E']
  x_F, y_F = pts['F']
  x_G, y_G = pts['G']
  x_H, y_H = pts['H']
  x_I, y_I = pts['I']

  A = distance_3_points(pts['D'], pts['A'], pts['C'])
  B = np.array([0, 2*v1])
  M_D = inv(A) @ B

  M_E = (k.LEGS[k.FL]['lengths']['ae']/(k.LEGS[k.FL]['lengths']['ae'] - k.LEGS[k.FL]['lengths']['de'])) * M_D

  A = distance_3_points(pts['F'], pts['E'], pts['B'])
  B = np.array([
    [2*(x_F - x_E), 2*(y_F - y_E)], 
    [0, 0]])
  M_F = (inv(A) @ B) @ M_E

  M_G =((k.LEGS[k.FL]['lengths']['ef']+k.LEGS[k.FL]['lengths']['fg']) / k.LEGS[k.FL]['lengths']['ef']) * M_F

  M_H = (BH/BF) * M_F

  A = distance_3_points(pts['I'], pts['G'], pts['H'])
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
    [V1[1], V2[1]]])

  Jacob = inv((k.LEGS[k.FL]['lengths']['gj']/k.LEGS[k.FL]['lengths']['gi']) * M_I)

  return Jacob

def gen_matrix_leg(v1, v2, v3, alpha, leg_id): # ATTENTION A ORI
  '''
  Génère la matrice de passage de (dx, dy, dz) vers (dv1, dv2, dv3)
  '''
  pts = k.get_leg_points_V1_V2(v1/1000, v2/1000, k.LEGS[k.FL]['lengths'])

  A = mat_A(pts, v1, v2, v3, alpha)
  B = mat_B(pts, alpha, leg_id)
  
  return A @ inv(B)

def mat_A(pts, v1, v2, v3, alpha):  
  Jacob = gen_jacob_plan(pts, v1/1000, v2/1000)

  A = np.array([
  [Jacob[0][0], Jacob[0][1], 0],
  [Jacob[1][0], Jacob[1][1], 0],
  [0, 0, -np.sin(alpha)/v3]])

  return A

def mat_B(pts, alpha, leg_id):
  X = pts['J'][0]*1000
  Z = pts['J'][1]*1000

  B = np.array([
  [np.cos(alpha), 0, -X * np.sin(alpha)],
  [np.sin(alpha), 0, Z * np.cos(alpha)],
  [0, 1, 0]])

  return B

def cos_angle_to_v3(cangle):
  '''
  Fonction auxiliaire de move_xyz
  Retourne l'élongation de v3 en fonction de l'angle de la patte au chassis  

  >>> cos_angle_to_v3(np.cos(np.pi/4)) - al_kashi_longueur(KO, LO, np.pi/4 - np.arccos(MO/LO))
  0.0
  >>> v3_to_cos_angle(cos_angle_to_v3(np.cos(np.pi/4))) - np.cos(np.pi/4) < 0.0000001 
  True
  '''
  return np.sqrt(LO**2 + KO**2 - 2*LO*KO * (cangle * MO/LO + np.sqrt(1-cangle**2)*np.sqrt(1-(MO/LO)**2)))

def v3_to_cos_angle(v3):
  '''
  Fonction auxiliaire de move_xyz
  Retourne l'angle de la patte au chassis en fonction de l'élongation de v3
      
  >>> v3_to_cos_angle(500) - np.cos(al_kashi_angle(LO, KO, 500) + np.arccos(MO/LO))
  0.0
  >>> cos_angle_to_v3(v3_to_cos_angle(500)) - 500 < 0.0000001
  True
  '''
  return (KO**2 + LO**2 - v3**2)/(2*KO*LO) * MO/LO - np.sqrt(1 - ((KO**2 + LO**2 - v3**2)/(2*KO*LO))**2) * np.sqrt(1-(MO/LO)**2)




def move_xyz(x, y, z, v1, v2, v3, dstep, p, eps, leg_id):
  '''
  Retourne la liste des élongations successives des verrins permettant de placer le bout de la patte en (x,y,z)
  '''
  L = []
  c = 0

  pts = k.get_leg_points_V1_V2(v1/1000, v2/1000, k.LEGS[k.FL]['lengths'])
  X, Z = pts['J'][0]*1000, pts['J'][1]*1000
  x0, y0, z0 = d2_to_d3(X, Z, v3_to_cos_angle(v3))
  dist = distance(x - x0, y - y0, z - z0)
  while dist > eps and c < 300:
    c += 1
    U = np.array([(x - x0), (y - y0), (z - z0)])
    U = dstep / 100 * U # / np.linalg.norm(U)**p
    dx, dy, dz = U[0], U[1], U[2]

    v1, v2, v3 = solve_indirect_cyl(x0+dx, y0+dy, z0+dz, x0, y0, z0, v1, v2, v3, leg_id, pts)
    L.append((v1, v2, v3))
    
    pts = k.get_leg_points_V1_V2(v1/1000, v2/1000, k.LEGS[k.FL]['lengths'])
    X, Z = pts['J'][0]*1000, pts['J'][1]*1000
    x1, y1, z1 = x0+dx, y0+dy, z0+dz
    x0, y0, z0 = d2_to_d3(X, Z, v3_to_cos_angle(v3))

    dist = distance(x - x0, y - y0, z - z0)

  return L



def direct_xyz(v1 ,v2, v3, leg_id):
  '''
  Retourne les positions x, y, z du bout de la patte en fonctions des v1, v2, v3
  '''
  X, Z = k.get_leg_points_V1_V2(v1/1000, v2/1000, k.LEGS[k.FL]['lengths'])['J']
  calpha = v3_to_cos_angle(v3)
  x = ori[leg_id][0] * X * np.sqrt(1-calpha**2) * 1000
  y = ori[leg_id][1] * X * calpha * 1000
  z = Z * 1000
  return x, y, z



def distance(x, y, z=0):
  '''
  Calcule la distance euclidienne dans l'espace en 3 dimensions (ou 2D selon le nombre de coordonnées passées en paramètre)
  
  >>> distance(0, 0, 0)
  0.0
  >>> distance(1, 0, 0)
  1.0
  >>> distance(0, 1, 0)
  1.0
  >>> distance(0, 1, 0)
  1.0
  >>> distance(0, 0, 1)
  1.0
  >>> np.abs(distance(1, 1, 0) - np.sqrt(2)) < 0.0000001 
  True
  '''
  return np.sqrt(x**2+y**2+z**2)



def al_kashi_longueur(a, b, alpha):
  return np.sqrt(a**2 + b**2 - 2*a*b*np.cos(alpha))

def al_kashi_angle(a, b, c):
  return np.arccos((a**2 + b**2 - c**2) / (2*a*b))



################################ DEPRECATED ################################

def distance_euclidienne(xi, yi, xj, yj):
  return np.sqrt((xj - xi)**2 + (yj - yi)**2)

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

beta = np.arccos(MO/LO)

def angle_to_v3(angle):
  '''
  Fonction auxiliaire de move_xyz
  Retourne l'élongation de v3 en fonction de l'angle de la patte au chassis
    
  >>> angle_to_v3(np.pi/4) - al_kashi_longueur(KO, LO, np.pi/4 - beta)
  0.0
  >>> v3_to_angle(angle_to_v3(np.pi/4)) - np.pi/4 < 0.0000001 
  True
  '''
  return np.sqrt(LO**2 + KO**2 - 2*LO*KO*np.cos(angle - beta))

def v3_to_angle(v3):
  '''
  Fonction auxiliaire de move_xyz
  Retourne l'angle de la patte au chassis en fonction de l'élongation de v3
      
  >>> v3_to_angle(500) - (al_kashi_angle(LO, KO, 500) + beta)
  0.0
  >>> angle_to_v3(v3_to_angle(500)) - 500 < 0.0000001 
  True
  '''
  return np.arccos((KO**2 + LO**2 - v3**2)/(2*KO*LO)) + beta

############################################################################

if __name__ == "__main__":
    import doctest
    doctest.testmod()

def gen_jacob_direct(pts, v1, v2):
  '''
  Retourne la Jacobienne correspondant au modèle cinématique indirect dans le plan de la patte
  Prend en argument la position des points de la patte et l'élongation des verrins en m

  >>> gen_jacob_plan(k.get_leg_points_V1_V2(0.495, 0.585, k.LEGS[k.FL]['lengths']), 0.495, 0.585) @ np.array([0, 0])
  array([0., 0.])
  '''
  x_E, y_E = pts['E']
  x_F, y_F = pts['F']
  x_G, y_G = pts['G']
  x_H, y_H = pts['H']
  x_I, y_I = pts['I']

  A = distance_3_points(pts['D'], pts['A'], pts['C'])
  B = np.array([0, 2*v1])
  M_D = inv(A) @ B

  M_E = (k.LEGS[k.FL]['lengths']['ae']/(k.LEGS[k.FL]['lengths']['ae'] - k.LEGS[k.FL]['lengths']['de'])) * M_D

  A = distance_3_points(pts['F'], pts['E'], pts['B'])
  B = np.array([
    [2*(x_F - x_E), 2*(y_F - y_E)],
    [0, 0]])
  M_F = (inv(A) @ B) @ M_E

  M_G =((k.LEGS[k.FL]['lengths']['ef']+k.LEGS[k.FL]['lengths']['fg']) / k.LEGS[k.FL]['lengths']['ef']) * M_F

  M_H = (BH/BF) * M_F

  A = distance_3_points(pts['I'], pts['G'], pts['H'])
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
    [V1[1], V2[1]]])

  return (k.LEGS[k.FL]['lengths']['gj']/k.LEGS[k.FL]['lengths']['gi']) * M_I

