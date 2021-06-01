import numpy as np

alpha = 1

x_A = 1
y_A = 1
x_B = 1
y_B = 1

AB = 1
AC = 1
AD = 1
AE = 1
BH = 1
BF = 1
EF = 1
EG = 1
FH = 1
FG = 1
GI = 1
GJ = 1


def distance_euclidienne(xi, yi, xj, yj):
  return np.sqrt((xj - xi)**2 + (yj - yi)**2)

def al_kashi_longueur(a,b,alpha):
  return np.sqrt(a**2 + b**2 + 2*a*b*np.cos(alpha))

def al_kashi_angle(a, b, c):
  return np.arccos((a**2 + b**2 - c**2) / (2*a*b))

def direct(v1, v2):
  theta1 = alpha + al_kashi_angle(AD, AC, v1)
  x_E = AE * np.cos(theta1) + x_A
  y_E = AE * np.sin(theta1) + y_A
  EB = distance_euclidienne(x_E, y_E, x_B, y_B)

  theta2 = al_kashi_angle(EF, EB, BF)
  theta3 = al_kashi_angle(AE, EB, AB)

  beta = theta2 + theta3 - (np.pi - theta1)
  x_F = x_E + EF * np.cos(beta)
  y_F = y_E + EF * np.sin(beta)
  x_G = x_E + EG * np.cos(beta)
  y_G = y_E + EG * np.sin(beta)
  
  x_H = BH / BF * (x_F - x_B)
  y_H = BH / BF * (y_F - y_B)
  GH = distance_euclidienne(x_G, x_H, y_G, y_H)
  
  theta4 = al_kashi_angle(FG, GH, FH)
  theta5 = al_kashi_angle(GI, GH, v2)
  theta6 = np.pi - (theta4 + theta5 + beta)
  
  return x_G + GJ * np.cos(theta6), y_G - GJ * np.sin(theta6)
  