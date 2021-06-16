# L'ensemble des distances sont exprimées en mm : segments de patte et élongations des vérins

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import inv
from qpsolvers import solve_qp
import kinetic as kin
from upg_kinetic import *

x_A = - kin.LEGS[kin.FL]['lengths']['ao']
y_A = 0
x_B = kin.LEGS[kin.FL]['lengths']['bo']
y_B = 0

AB = kin.LEGS[kin.FL]['lengths']['ao'] + kin.LEGS[kin.FL]['lengths']['bo']
AC = np.sqrt((AB + kin.LEGS[kin.FL]['lengths']['bcx'])**2 + kin.LEGS[kin.FL]['lengths']['bcy']**2)
BC = np.sqrt(kin.LEGS[kin.FL]['lengths']['bcx']**2 + kin.LEGS[kin.FL]['lengths']['bcy']**2)
AE = kin.LEGS[kin.FL]['lengths']['ae']
AD = AE - kin.LEGS[kin.FL]['lengths']['de']
BF = kin.LEGS[kin.FL]['lengths']['bf']
FH = kin.LEGS[kin.FL]['lengths']['fh']
BH = BF - FH
FG = kin.LEGS[kin.FL]['lengths']['fg']
EF = kin.LEGS[kin.FL]['lengths']['ef']
EG = EF + FG
GI = kin.LEGS[kin.FL]['lengths']['gi']
GJ = kin.LEGS[kin.FL]['lengths']['gj']

KO = kin.LEGS[kin.FL]['lengths']['yaw_c']
LM = kin.LEGS[kin.FL]['lengths']['yaw_b']
MO = kin.LEGS[kin.FL]['lengths']['yaw_a']
LO = np.sqrt(LM**2 + MO**2)

ori = np.array([[1, 1, -1, -1], [1, -1, -1, 1]]) #[[oritentation selon x][orientation selon y]]

V = 460, 565, 500



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