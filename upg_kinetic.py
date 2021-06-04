# L'ensemble des distances sont exprimées en mm : segments de patte et élongations des vérins

import numpy as np
from numpy.linalg import inv
from numpy import dot

import kinetic as k

x_A = - k.LEGS[k.FR]['lengths']['ao']
y_A = 0
x_B = k.LEGS[k.FR]['lengths']['bo']
y_B = 0

AB = k.LEGS[k.FR]['lengths']['ao'] + k.LEGS[k.FR]['lengths']['bo']
AC = np.sqrt((AB + k.LEGS[k.FR]['lengths']['bcx'])**2 + k.LEGS[k.FR]['lengths']['bcy']**2)
BC = np.sqrt(k.LEGS[k.FR]['lengths']['bcx']**2 + k.LEGS[k.FR]['lengths']['bcy']**2)
AE = k.LEGS[k.FR]['lengths']['ae']
AD = AE - k.LEGS[k.FR]['lengths']['de']
BF = k.LEGS[k.FR]['lengths']['bf']
FH = k.LEGS[k.FR]['lengths']['fh']
BH = BF - FH
FG = k.LEGS[k.FR]['lengths']['fg']
EF = k.LEGS[k.FR]['lengths']['ef']
EG = EF + FG
GI = k.LEGS[k.FR]['lengths']['gi']
GJ = k.LEGS[k.FR]['lengths']['gj']

KO = k.LEGS[k.FR]['lengths']['yaw_c']
LM = k.LEGS[k.FR]['lengths']['yaw_b']
MO = k.LEGS[k.FR]['lengths']['yaw_a']
LO = np.sqrt(LM**2 + MO**2)
beta = np.arctan(LM/MO)

V = 460, 565, 500

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
  #print(angle)
  return np.sqrt(LO**2 + KO**2 - 2*LO*KO*np.cos(angle - beta))

'''
Fonction auxiliaire de move_xyz
Retourne l'angle de la patte au chassis en fonction de l'élongation de v3
'''
def v3_to_angle(v3):
  #print(v3)
  return np.arccos((KO**2 + LO**2 - v3**2)/(2*KO*LO)) + beta

'''
Retourne la liste des élongations successives des verrins permettant de placer le bout de la patte en (x,y,z)
'''
def move_xyz(x, y, z, dstep, eps):
  L = []
  c = 0
  X = np.sqrt(x**2 + y**2)
  v1, v2, v3 = get_ver()
  pts = k.get_leg_points_V1_V2(v1/1000, v2/1000, k.LEGS[k.FR]['lengths'])
  cur_X, cur_Y = pts['J'][0]*1000, pts['J'][1]*1000
  cur_angle = v3_to_angle(v3)
  cur_x, cur_y = cur_X * np.cos(cur_angle), - cur_X * np.sin(cur_angle)
  dist = distance(x - cur_x, y - cur_y, z - cur_Y)
  while dist > eps and c < 1000:
    c += 1
    # on détermine dx, dy, dz tq conforme à dstep
    theta1 = np.arctan((y - cur_y)/(x - cur_x))
    theta2 = np.arctan((z - cur_Y)/(X - cur_X))
    print(theta1, theta2)
    dx, dy, dz = deltas(theta1, theta2, dstep)
    #print(dx, dy, dz)
    # on calcule les nouveaux v1, v2, v3 (et on les ajoute à L)
    new_v3 = angle_to_v3(np.arctan((cur_y + dy)/(cur_x + dx))) 
    dX = np.sqrt(dx**2 + dy**2)
    #print(dz)
    deltaXY = np.array([dX/1000, dz/1000])
    Jacob = gen_jacob(pts, v1/1000, v2/1000)
    #print(Jacob)
    dv1v2 = Jacob @ deltaXY
    #print(dv1v2)
    v1 += dv1v2[0]*1000
    v2 += dv1v2[1]*1000
    v3 = new_v3
    #print(v1, v2, v3)
    L.append((v1, v2, v3))

    # on recalcule les nouvelles positions courantes et la distance à l'objectif
    pts = k.get_leg_points_V1_V2(v1/1000, v2/1000, k.LEGS[k.FR]['lengths'])
    cur_X, cur_Y = pts['J'][0]*1000, pts['J'][1]*1000
    cur_angle = v3_to_angle(v3)
    cur_x, cur_y = cur_X * np.cos(cur_angle), cur_X * np.sin(cur_angle)
    dist = distance(x - cur_x, y - cur_y, z - cur_Y)
  return L



'''
Retourne les positions x, y, z du bout de la patte en fonctions des v1, v2, v3
'''
def direct_xyz(v1 ,v2, v3):
  X, Y = k.get_leg_points_V1_V2(v1/1000, v2/1000, k.LEGS[k.FR]['lengths'])['J']
  #print (X, Y)
  alpha = v3_to_angle(v3)
  x = X * np.cos(alpha) * 1000
  y = -X * np.sin(alpha) * 1000
  z = Y * 1000
  return x, y, z



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
  #return (controlers[l].la[1]['position']+450.0)/1000.0, (controlers[l].la[0]['position']+450.0)/1000.0, (controlers[l].la[2]['position']+450.0)/1000.0
  return V

def set_ver(v1, v2, v3):
  global V
  V = v1, v2, v3



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
