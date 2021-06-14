import numpy as np
from numpy.linalg import inv
from numpy import dot
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d  # Fonction pour la 3D
import time

import kinetic as k
from upg_kinetic import *

############################# FONCTIONS DE TEST ################################
def test_jacob(v1, v2, dstep_x, dstep_y):
    '''
      Retourne l'erreur relative en x et en y de l'application de la Jacobienne du modèle indirect
      Prend en argument l'élongation des vérins et la distance de déplacement selon x et y en m
    '''

    pts = k.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, k.LEGS[k.FR]['lengths'])
    x, y = pts['J'][0], pts['J'][1] 
    Jacob = gen_jacob_plan(pts, v1 / 1000, v2 / 1000)

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


def test_move_xyz(x0, y0, z0, x, y, z, dstep, p, eps, leg_id):
    '''
    Trace la trajectoire de la patte du robot entre 2 points de l'espace en suivant move_xyz
    '''
    # Théorique
    Xt = [x0, x]
    Yt = [y0, y]
    Zt = [z0, z]

    # Positions des vérins
    L = move_xyz(x, y, z, V[0], V[1], V[2], dstep, p, eps, leg_id)
    V1 = [v[0] for v in L]
    V2 = [v[1] for v in L]
    V3 = [v[2] for v in L]
    T = np.linspace(0, 1, len(L))

    # Positions du bout de la patte
    Pos = list(map(direct_xyz, [v[0] for v in L], [v[1] for v in L], [v[2] for v in L], [k.FL for v in L]))
    Xp = [p[0] for p in Pos]
    Yp = [p[1] for p in Pos]
    Zp = [p[2] for p in Pos]

    # Tracé de la position du bout de la patte
    fig = plt.figure()
    ax = fig.gca(projection='3d')  # Affichage en 3D
    ax.plot(Xt, Yt, Zt, label='Théorique')  # Tracé de la courbe 3D
    ax.scatter(Xp, Yp, Zp, label='Positions',
               marker='d')  # Tracé des points 3D
    plt.title("Trajectoire patte FL avec (dstep, p, eps) = (" + str(dstep) + ", " + str(p) + ", " + str(eps) + ")")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.tight_layout()
    plt.show()

    # Tracé des positions des vérins
    plt.plot(T, V1, label='V1' )
    plt.plot(T, V2, label='V2' )
    plt.plot(T, V3, label='V3' )
    plt.title("Elongations des vérins dans le mouvement")
    # plt.plot.set_xlabel('L (en mm)')
    # plt.plot.set_ylabel('T')
    plt.show()

    # fig = plt.figure()
    # fig.subplots_adjust(top=0.8)
    # ax = fig.add_subplot(221)




    # ax1.set_ylabel('volts')
    # ax1.set_title('a sine wave')

    # t = np.arange(0.0, 1.0, 0.01)
    # s = np.sin(2 * np.pi * t)
    # line, = ax1.plot(t, s, lw=2)

    # ax2 = fig.add_subplot(212)
    # ax2.set_xlabel('time (s)')


def test_A(dX, dZ, dalpha, v1, v2, v3):
  alpha = v3_to_angle(v3)
  print (alpha)
  pts = k.get_leg_points_V1_V2(v1/1000, v2/1000, k.LEGS[k.FL]['lengths'])
  A = mat_A(pts, v1, v2, v3, alpha)

  dPos = np.array([dX/1000, dZ/1000, dalpha])
  dV = A @ dPos

  new_v1, new_v2 = v1 + dV[0]*1000, v2 + dV[1]*1000
  new_pts = k.get_leg_points_V1_V2(new_v1/1000, new_v2/1000, k.LEGS[k.FL]['lengths'])

  new_alpha = v3_to_angle(v3 + dV[2]*1000)
  print (dV)
  print(new_alpha)
  err_rel = (abs(new_pts['J'][0] - pts['J'][0])*1000 - dX)/dX, (abs(new_pts['J'][1] - pts['J'][1])*1000 - dZ)/dZ, (abs(new_alpha - alpha) - dalpha)/dalpha

  return err_rel


def test_jacob_direct():
    """
    Comparaison des resultats en J du modèle direct utilisant une jacobienne et de celui utilisant la géométrie
    """
    leg_id = k.FL
    V = np.array([495, 585, 515])
    x0, y0, z0 = direct_xyz(V[0], V[1], V[2], leg_id)
    print("position de départ :", x0, y0, z0)
    deltaV = np.array([5, -5, 5])
    V += deltaV
    deltaV2 = np.array([deltaV[0], deltaV[1]])
    Jacob = gen_jacob_direct(k.get_leg_points_V1_V2(V[0]/1000, V[1]/1000, k.LEGS[leg_id]['lengths']), V[0]/1000, V[1]/1000)
    deltaX = Jacob @ deltaV2
    calpha = v3_to_cos_angle(V[2])
    deltax = np.array([
        deltaX[0] * np.sqrt(1-calpha**2),
        deltaX[0] * calpha,
        deltaX[1]
    ])
    print("position selon jacob_direct :", x0+deltax[0], y0+deltax[1], z0+deltax[2])
    print("position selon direct julien :", direct_xyz(V[0], V[1], V[2], leg_id))
    print("")


def test_precision_jacobienne_en_direct(v1, v2, dv1, dv2, leg_id):
    """
    Test des précisions pas à pas pour chaque point de la patte en suivant une méthode géométrique exacte et une methode
    approchée utilisant la jacobienne
    @param v1: position initiale du vérin 1
    @param v2: position initiale du vérin 2
    @param dv1: déplacement du vérin 1
    @param dv2: déplacement du vérin 2
    @param leg_id: id de la patte
    @return:
    """
    pts = k.get_leg_points_V1_V2(v1, v2, k.LEGS[leg_id]['lengths'])
    x_E, y_E = pts['E']
    x_F, y_F = pts['F']
    x_G, y_G = pts['G']
    x_H, y_H = pts['H']
    x_I, y_I = pts['I']
    print('initial D :', pts['D'])
    pts_end = k.get_leg_points_V1_V2((v1+dv1), (v2+dv2), k.LEGS[leg_id]['lengths'])
    A = distance_3_points(pts['D'], pts['A'], pts['C'])
    B = np.array([0, 2 * v1])
    M_D = inv(A) @ B
    deltaD = M_D * dv1
    D = pts['D'] + deltaD
    print("real D :", pts_end['D'])
    print("jacob D :", D, "\n")

    M_E = (k.LEGS[k.FL]['lengths']['ae']/(k.LEGS[k.FL]['lengths']['ae'] - k.LEGS[k.FL]['lengths']['de'])) * M_D
    deltaE = M_E * dv1
    E = pts['E'] + deltaE
    print('initial E :', pts['E'])
    print("real E :", pts_end['E'])
    print("jacob E :", E, "\n")

    A = distance_3_points(pts['F'], pts['E'], pts['B'])
    B = np.array([
        [2 * (x_F - x_E), 2 * (y_F - y_E)],
        [0, 0]])
    M_F = (inv(A) @ B) @ M_E
    deltaF = M_F * dv1
    F = pts['F'] + deltaF
    print('initial F:', pts['F'])
    print("real F:", pts_end['F'])
    print("jacob F", F, "\n")

    M_G = ((k.LEGS[k.FL]['lengths']['ef'] + k.LEGS[k.FL]['lengths']['fg']) / k.LEGS[k.FL]['lengths']['ef']) * M_F
    deltaG = M_G * dv1
    G = pts['G'] + deltaG
    print('initial G:', pts['G'])
    print("real G:", pts_end['G'])
    print("jacob G", G, "\n")

    M_H = (BH / BF) * M_F
    deltaH = M_H *  dv1
    H = pts['H'] + deltaH
    print('initial H:', pts['H'])
    print("real H:", pts_end['H'])
    print("jacob H", H, "\n")

    A = distance_3_points(pts['I'], pts['G'], pts['H'])
    B = np.array([
        [2 * (x_I - x_G), 2 * (y_I - y_G)],
        [0, 0]])
    C = np.array([
        [0, 0],
        [2 * (x_I - x_H), 2 * (y_I - y_H)]])
    D = np.array([0, 2 * v2])
    V1 = inv(A) @ (B @ M_G + C @ M_H)
    V2 = inv(A) @ D
    M_I = np.array([
        [V1[0], V2[0]],
        [V1[1], V2[1]]])
    deltaI = M_I @ np.array([dv1, dv2])
    I = pts['I'] + deltaI
    print('initial I:', pts['I'])
    print("real I:", pts_end['I'])
    print("jacob I", I, "\n")

    M_J = (k.LEGS[k.FL]['lengths']['gj']/k.LEGS[k.FL]['lengths']['gi']) * M_I
    deltaJ = M_J @ np.array([dv1, dv2])
    J = pts['J'] + deltaJ
    print('initial J:', pts['J'])
    print("real J:", pts_end['J'])
    print("jacob J", J, "\n")

def test_comparaison_minimize_vs_jacob_indirect():
    """
    test de comparaison des méthodes de cinématique indirecte
    """
    v1, v2 = 0.495, 0.575
    print("initialement v1, v2 =", v1, v2)
    leg_id = k.FL
    lpl = k.LEGS[leg_id]['lengths']
    pts = k.get_leg_points_V1_V2(v1, v2, lpl)
    x0, y0 = pts['J']
    dX = np.array([0.010, -0.005])
    x, y = x0 + dX[0], y0 + dX[1]

    r = k.inverse_kinetic(x, y, lpl)
    v1_minimize, v2_minimize = r[0], r[1]
    print("selon minimize, V =", r[0], r[1])

    J = gen_jacob_plan(pts, v1, v2)
    dV = J @ dX
    print("selon la methode jacobienne, V =", V[0]/1000, V[1]/1000)

################################# TESTS ####################################

t_move = 1
# test de move_xyz
if t_move:
    x0, y0, z0 = direct_xyz(V[0], V[1], V[2], k.FL)
    print(x0, y0, z0)
    test_move_xyz(x0, y0, z0, 750, y0, z0, 1, 1, 5, k.FL)
    test_move_xyz(x0, y0, z0, x0, 700, z0, 1, 1, 5, k.FL)
    test_move_xyz(x0, y0, z0, x0, y0, -550, 1, 1, 5, k.FL)
    test_move_xyz(x0, y0, z0, 750, 600, -400, 1, 1, 5, k.FL)

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

t_a = 0
if t_a:
    print (test_A(1, -1, 0.1, V[0], V[1], V[2]))

t_jacob_direct = 0
if t_jacob_direct:
    test_jacob_direct()
    test_precision_jacobienne_en_direct(0.495, 0.555, 0.015, -0.015, k.FL)

test_comp_indirect = 0
if test_comp_indirect:
    test_comparaison_minimize_vs_jacob_indirect()
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
        k.get_leg_points_V1_V2(495 / 1000, 585 / 1000 , k.LEGS[k.FR]['lengths'])['J']
    t2 = time.time() - t

    print("direct_v1 prend ", t1 * 100, " us")
    print("L'algo de Julien prend ", t2 * 100, " us")

    print("direct v1 retourne : ", direct_v1(495, 585))
    print("L'algo de Julien retourne :", k.get_leg_points_V1_V2(495 / 1000, 585 / 1000 , k.LEGS[k.FR]['lengths'])['J'])

############################################################################

