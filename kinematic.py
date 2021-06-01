import numpy as np
import kinetic as k


x_A = k.LEG_PARTS_LENGTHS['ao']
y_A = 0
x_B = k.LEG_PARTS_LENGTHS['bo']
y_B = 0

AB = k.LEG_PARTS_LENGTHS['ao'] + k.LEG_PARTS_LENGTHS['bo']
AC = np.sqrt((AB + k.LEG_PARTS_LENGTHS['bcx'])**2 + k.LEG_PARTS_LENGTHS['bcy']**2)
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

def distance_euclidienne(xi, yi, xj, yj):
  return np.sqrt((xj - xi)**2 + (yj - yi)**2)

def al_kashi_longueur(a,b,alpha):
  return np.sqrt(a**2 + b**2 + 2*a*b*np.cos(alpha))

def al_kashi_angle(a, b, c):
  print(a,b,c)
  return np.arccos((a**2 + b**2 - c**2) / (2*a*b))

alpha = al_kashi_angle(AB, AC, BC)
print(alpha *180/np.pi)

def direct(v1, v2):
  theta1 = alpha + al_kashi_angle(AD, AC, v1)
  x_E = AE * np.cos(theta1) + x_A
  y_E = AE * np.sin(theta1) + y_A
  EB = distance_euclidienne(x_E, y_E, x_B, y_B)

  theta2 = al_kashi_angle(EF, EB, BF)
  theta3 = al_kashi_angle(AE, EB, AB)

  beta = theta2 + theta3 - (np.pi - theta1)
  print(beta)
  x_F = x_E + EF * np.cos(beta)
  y_F = y_E + EF * np.sin(beta)
  x_G = x_E + EG * np.cos(beta)
  y_G = y_E + EG * np.sin(beta)

  x_H = (BH / BF * (x_F - x_B)) + x_B
  y_H = (BH / BF * (y_F - y_B)) + y_B
  print(x_B, y_B)
  print(x_F, y_F)
  print(x_H, y_H)

  GH = distance_euclidienne(x_G, x_H, y_G, y_H)
  
  theta4 = al_kashi_angle(FG, GH, FH)
  theta5 = al_kashi_angle(GI, GH, v2)
  theta6 = np.pi - (theta4 + theta5 + beta)
  
  return x_G + GJ * np.cos(theta6), y_G - GJ * np.sin(theta6)

print(direct(470,580))
