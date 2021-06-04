import numpy as np
from numpy.linalg import inv
from numpy import dot
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d  # Fonction pour la 3D
import time

import kinetic as k
from upg_kinetic import *

'''
  Retourne l'erreur relative en x et en y de l'application de la Jacobienne du modèle indirect
  Prend en argument l'élongation des vérins et la distance de déplacement selon x et y en m
'''
def test_jacob(v1, v2, dstep_x, dstep_y):
   
  pts = k.get_leg_points_V1_V2(v1, v2, k.LEG_PARTS_LENGTHS)
  x, y = pts['J']
  Jacob = gen_jacob(pts, v1, v2)
  Dpos = np.array([dstep_x, dstep_y])
  DV = np.dot(Jacob, Dpos)
  
  pts = k.get_leg_points_V1_V2(v1 + DV[0], v2 + DV[1], k.LEG_PARTS_LENGTHS)
  print("position initiale : ", x, y)
  
  new_x, new_y = pts['J']
  print("nouvelle position : ", new_x, new_y)
  print("on s'est déplacé de ", new_x-x, " mm en x et de ", new_y-y, " mm en y")

  err_pos = np.array([new_x - (x + dstep_x), new_y - (y + dstep_y)])
  print("erreur de position en x : ", err_pos[0] * 1000, " mm")
  print("erreur de position en y : ", err_pos[1] * 1000, " mm")  
  err_rel = np.array([None, None])
  if dstep_x != 0 : err_rel[0] = err_pos[0] / dstep_x
  if dstep_y != 0 : err_rel[1] = err_pos[1] / dstep_y
  if dstep_x != 0 : print("erreur relative en x : ", err_rel[0])
  if dstep_y != 0 : print("erreur relative en y : ", err_rel[1])
  print("\n")

  return err_rel

'''
Trace la trajectoire de la patte du robot entre 2 points de l'espace en suivant move_xyz
'''
def test_move_xyz(x0, y0, z0, x, y, z, dstep, eps):
  # Tableau pour les 3 axes
  X = [x0, x]
  Y = [y0, y]
  Z = [z0, z]

  # Tracé du résultat en 3D
  fig = plt.figure()
  ax = fig.gca(projection='3d')  # Affichage en 3D
  ax.plot(X, Y, Z, label='Courbe')  # Tracé de la courbe 3D
  plt.title("Courbe 3D")
  ax.set_xlabel('X')
  ax.set_ylabel('Y')
  ax.set_zlabel('Z')
  plt.tight_layout()
  plt.show()



################################# TESTS ####################################

t_jacob = 0
# test de l'erreur relative de position en appliquant la Jacobienne
if t_jacob == 1:
  test_jacob(495/1000, 585/1000, 10/1000, 0)
  test_jacob(505/1000, 585/1000, 10/1000, 0)
  test_jacob(515/1000, 585/1000, 10/1000, 0)

  test_jacob(495/1000, 585/1000, -10/1000, 0)
  test_jacob(505/1000, 585/1000, -10/1000, 0)
  test_jacob(515/1000, 585/1000, -10/1000, 0)

  test_jacob(495/1000, 585/1000, 0, 10/1000)
  test_jacob(505/1000, 585/1000, 0, 10/1000)
  test_jacob(515/1000, 585/1000, 0, 10/1000)

  test_jacob(495/1000, 585/1000, 0, -10/1000)
  test_jacob(505/1000, 585/1000, 0, -10/1000)
  test_jacob(515/1000, 585/1000, 0, -10/1000)

  test_jacob(495/1000, 585/1000, 10/1000, 10/1000)
  test_jacob(505/1000, 585/1000, 10/1000, 10/1000)
  test_jacob(515/1000, 585/1000, 10/1000, 10/1000)

  test_jacob(495/1000, 585/1000, 10/1000, -10/1000)
  test_jacob(505/1000, 585/1000, 10/1000, -10/1000)
  test_jacob(515/1000, 585/1000, 10/1000, -10/1000)

  test_jacob(495/1000, 585/1000, -10/1000, 10/1000)
  test_jacob(505/1000, 585/1000, -10/1000, 10/1000)
  test_jacob(515/1000, 585/1000, -10/1000, 10/1000)

  test_jacob(495/1000, 585/1000, -10/1000, -10/1000)
  test_jacob(505/1000, 585/1000, -10/1000, -10/1000)
  test_jacob(515/1000, 585/1000, -10/1000, -10/1000)

############################################################################



################################ DEPRECATED ################################

t_comp = 0
# comparaison en temps de nos modèles directs (v1 / Julien)
if t_comp == 1:
  t = time.time()
  for i in range (10000):
    direct_v1(495, 585)
  t1 = time.time() - t

  t = time.time()
  for i in range (10000):
    k.get_leg_points_V1_V2(495/1000, 585/1000, k.LEG_PARTS_LENGTHS)['J']
  t2 = time.time() - t

  print("direct_v1 prend ", t1*100, " us")
  print("L'algo de Julien prend ", t2*100, " us")

############################################################################