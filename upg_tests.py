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

    pts = k.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, k.LEGS[k.FR]['lengths'])
    x, y = pts['J'][0], pts['J'][1] 
    Jacob = gen_jacob(pts, v1 / 1000, v2 / 1000)

    Dpos = np.array([dstep_x / 1000, dstep_y / 1000])
    DV = Jacob @ Dpos

    pts = k.get_leg_points_V1_V2(v1 / 1000 + DV[0], v2 / 1000 + DV[1], k.LEGS[k.FR]['lengths'])
    #print("position initiale : ", x * 1000, y * 1000)

    new_x, new_y = pts['J']
    #print("nouvelle position : ", new_x * 1000, new_y * 1000)
    print("on s'est déplacé de ", (new_x - x) * 1000, " mm en x et de ", (new_y - y) * 1000, " mm en y")

    err_pos = np.array([new_x - (x + dstep_x/1000), new_y - (y + dstep_y/1000)])
    #print("erreur de position en x : ", err_pos[0] * 1000, " mm")
    #print("erreur de position en y : ", err_pos[1] * 1000, " mm")
    err_rel = np.array([None, None])
    if dstep_x != 0: err_rel[0] = err_pos[0]*1000 / dstep_x
    if dstep_y != 0: err_rel[1] = err_pos[1]*1000 / dstep_y
    if dstep_x != 0: print("erreur relative en x : ", err_rel[0])
    if dstep_y != 0: print("erreur relative en y : ", err_rel[1])
    print("\n")

    return err_rel


'''
Trace la trajectoire de la patte du robot entre 2 points de l'espace en suivant move_xyz
'''
def test_move_xyz(x0, y0, z0, x, y, z, dstep, eps, leg_id):
    # Théorique
    Xt = [x0, x]
    Yt = [y0, y]
    Zt = [z0, z]

    #Positions des vérins
    L = move_xyz(x, y, z, dstep, eps, leg_id)
    #print (V, L[0])
    Pos = list(map(direct_xyz, [v[0] for v in L], [v[1] for v in L], [v[2] for v in L], [k.FL for v in L]))
    #print((x0, y0, z0), Pos[0])
    X = [p[0] for p in Pos]
    Y = [p[1] for p in Pos]
    Z = [p[2] for p in Pos]

    # Tracé du résultat en 3D
    fig = plt.figure()
    ax = fig.gca(projection='3d')  # Affichage en 3D
    ax.plot(Xt, Yt, Zt, label='Théorique')  # Tracé de la courbe 3D
    ax.scatter(X, Y, Z, label='Positions',
               marker='d')  # Tracé des points 3D
    plt.title("Trajectoire patte FR avec dstep = " + str(dstep))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.tight_layout()
    plt.show()


################################# TESTS ####################################

t_move = 1
# test de move_xyz
if t_move:
    x0, y0, z0 = direct_xyz(V[0], V[1], V[2], k.FL)
    print(x0, y0, z0)
    # test_move_xyz(x0, y0, z0, 750, y0, z0, 10, 10, k.FL)
    # test_move_xyz(x0, y0, z0, x0, 600, z0, 10, 10, k.FL)
    # test_move_xyz(x0, y0, z0, x0, y0, -550, 10, 10, k.FL)
    test_move_xyz(x0, y0, z0, 750, 600, -400, 10, 5, k.FL)
    test_move_xyz(x0, y0, z0, 750, 600, -550, 20, 5, k.FL)
    test_move_xyz(x0, y0, z0, 750, 600, -550, 5, 5, k.FL)
    # test_move_xyz(x0, y0, z0, 750, 600, -550, 2, 10, k.FL)
    # test_move_xyz(x0, y0, z0, 750, 600, -550, 1, 10, k.FL)

t_jacob = 0
# test de l'erreur relative de position en appliquant la Jacobienne
if t_jacob:
    test_jacob(495, 585, 0, 0)
    test_jacob(505, 585, 0, 0)
    test_jacob(515, 585, 0, 0)

    test_jacob(495, 585, 10, 0)
    test_jacob(505, 585, 10, 0)
    test_jacob(515, 585, 10, 0)

    # test_jacob(495, 585, -10, 0)
    # test_jacob(505, 585, -10, 0)
    # test_jacob(515, 585, -10, 0)

    test_jacob(495, 585, 0, 10)
    test_jacob(505, 585, 0, 10)
    test_jacob(515, 585, 0, 10)

    # test_jacob(495, 585, 0, -10)
    # test_jacob(505, 585, 0, -10)
    # test_jacob(515, 585, 0, -10)

    # test_jacob(495, 585, 10, 10)
    # test_jacob(505, 585, 10, 10)
    # test_jacob(515, 585, 10, 10)

    # test_jacob(495, 585, 10, -10)
    # test_jacob(505, 585, 10, -10)
    # test_jacob(515, 585, 10, -10)

    # test_jacob(495, 585, -10, 10)
    # test_jacob(505, 585, -10, 10)
    # test_jacob(515, 585, -10, 10)

    # test_jacob(495, 585, -10, -10)
    # test_jacob(505, 585, -10, -10)
    # test_jacob(515, 585, -10, -10)

############################################################################


################################ DEPRECATED ################################

t_comp = 0
# comparaison en temps de nos modèles directs (v1 / Julien)
if t_comp:
    t = time.time()
    for i in range(10000):
        direct_v1(495, 585)
    t1 = time.time() - t

    t = time.time()
    for i in range(10000):
        k.get_leg_points_V1_V2(495 / 1000, 585 / 1000 ,
                               k.LEGS[k.FR]['lengths'])['J']
    t2 = time.time() - t

    print("direct_v1 prend ", t1 * 100, " us")
    print("L'algo de Julien prend ", t2 * 100, " us")

############################################################################