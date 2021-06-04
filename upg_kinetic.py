import numpy as np
from numpy.linalg import inv
from numpy import dot

import kinetic as k

x_A = - k.LEG_PARTS_LENGTHS['ao']
y_A = 0
x_B = k.LEG_PARTS_LENGTHS['bo']
y_B = 0

AB = k.LEG_PARTS_LENGTHS['ao'] + k.LEG_PARTS_LENGTHS['bo']
AC = np.sqrt((AB + k.LEG_PARTS_LENGTHS['bcx'])**2 + k.LEG_PARTS_LENGTHS['bcy']**2)
BC = np.sqrt(k.LEG_PARTS_LENGTHS['bcx']**2 + k.LEG_PARTS_LENGTHS['bcy']**2)
AE = k.LEG_PARTS_LENGTHS['ae']
AD = AE - k.LEG_PARTS_LENGTHS['de']
BF = k.LEG_PARTS_LENGTHS['bf']
FH = k.LEG_PARTS_LENGTHS['fh']
BH = BF - FH
FG = k.LEG_PARTS_LENGTHS['fg']
EF = k.LEG_PARTS_LENGTHS['ef']
EG = EF + FG
GI = k.LEG_PARTS_LENGTHS['gi']
GJ = k.LEG_PARTS_LENGTHS['gj']



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
  M_D = dot(inv(A), B)

  M_E = (AE/AD) * M_D

  A = np.array([
    [2*(x_F - x_E), 2*(y_F - y_E)], 
    [2*(x_F - x_B), 2*(y_F - y_B)]])
  B = np.array([
    [2*(x_F - x_E), 2*(y_F - y_E)], 
    [0, 0]])
  M_F = dot(dot(inv(A), B), M_E)

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
  V1 = dot(inv(A),(dot(B, M_G) + dot(C, M_H)))
  V2 = dot(inv(A), D) 
  M_I = np.array([[V1[0], V2[0]], [V1[1], V2[1]]])

  Jacob = inv((GJ/GI) * M_I)

  return Jacob



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

'''
Fonction auxiliaire de move_xyz
Retourne l'élongation de v3 en fonction de l'angle de la patte au chassis
'''
def angle_to_v3(angle):
  return np.sqrt(LO**2 + KO**2 - 2*LO*KO*np.cos(angle))

'''
Fonction auxiliaire de move_xyz
Retourne l'angle de la patte au chassis en fonction de l'élongation de v3
'''
def v3_to_angle(v3):
  return np.arccos((KO**2 + LO**2 - v3**2)/(2*KO*LO))

'''
Retourne la liste des élongations successives des verrins permettant de placer le bout de la patte en (x,y,z)
'''
def move_xyz(x, y, z, dstep, eps):
  L = []
  X = np.sqrt(x**2 + y**2)
  v1, v2, v3 = get_ver()
  pts = k.get_leg_points_V1_V2(v1, v2, k.LEG_PARTS_LENGTHS)
  cur_X, cur_Y = pts['J']
  cur_angle = v3_to_angle(v3)
  cur_x, cur_y = cur_X * np.cos(cur_angle), cur_X * np.sin(cur_angle)
  dist = distance(x - cur_x, y - cur_y, z - cur_Y)
  while dist > eps:
    # on détermine dx, dy, dz tq conforme à dstep
    theta1 = np.arccos((x - cur_x)/(X - cur_X))
    theta2 = np.arccos((X - cur_X)/(z - cur_Y))
    dx, dy, dz = deltas(theta1, theta2, dstep)

    # on calcule les nouveaux v1, v2, v3 (et on les ajoute à L)
    dv3 = angle_to_v3(np.tan(dy/dx))
    dX = np.sqrt(dx**2 + dy**2)
    deltaXY = np.array([dX, dz])
    Jacob = gen_jacob(pts, v1, v2)
    dv1v2 = np.dot(Jacob, deltaXY)
    v1 += dv1v2[0]
    v2 += dv1v2[1]
    v3 += dv3
    L.append((v1, v2, v3))

    # on recalcule les nouvelles positions courantes et la distance à l'objectif
    pts = k.get_leg_points_V1_V2(v1, v2, k.LEG_PARTS_LENGTHS)
    cur_X, cur_Y = pts['J']
    cur_angle = v3_to_angle(v3)
    cur_x, cur_y = cur_X * np.cos(cur_angle), cur_X * np.sin(cur_angle)
    dist = distance(x - cur_x, y - cur_y, z - cur_Y)
  return L



'''
Retourne les positions x, y, z du bout de la patte en fonctions des v1, v2, v3
'''
def direct(v1 ,v2, v3):
  X, Y = k.get_leg_points_V1_V2(v1, v2, k.LEG_PARTS_LENGTHS)['J']
  alpha = v3_to_angle(v3)
  x = X * np.cos(alpha)
  y = Y * np.sin(alpha)
  return x, y, Y 



'''
Calcule la distance euclidienne dans l'espace en 3 dimensions (ou 2D selon le nombre de coordonnées passées en paramètre)
'''
def distance(x, y, z=0):
  return np.sqrt(x**2+y**2+z**2)

'''
Retourne les positions des verrins
'''
def get_ver():
  #return 485, 545, 515
  return (controlers[l].la[1]['position']+450.0)/1000.0, (controlers[l].la[0]['position']+450.0)/1000.0, (controlers[l].la[2]['position']+450.0)/1000.0



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

############################################################################
