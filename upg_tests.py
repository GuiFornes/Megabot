import numpy as np
from numpy.linalg import inv
from numpy import dot
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d  # Fonction pour la 3D
import time

import kinetic as kin
from upg_kinetic import *
from upg_deprecated import *

############################# FONCTIONS DE TEST ################################

def test_jacob_2(v1, v2, dstep_x, dstep_y):
    """
    Retourne l'erreur relative en x et en y de l'application de la Jacobienne du modèle indirect plan
    Prend en argument l'élongation des vérins et la distance de déplacement selon x et y en m
    """
    pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, kin.LEGS[kin.FR]['lengths'])
    x, y = pts['J'][0], pts['J'][1] 
    Jacob = gen_jacob_2(pts, v1 / 1000, v2 / 1000)

    Dpos = np.array([dstep_x / 1000, dstep_y / 1000])
    DV = Jacob @ Dpos

    pts = kin.get_leg_points_V1_V2(v1 / 1000 + DV[0], v2 / 1000 + DV[1], kin.LEGS[kin.FR]['lengths'])
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


def test_normalized_move_xyz(x0, y0, z0, x, y, z, dstep, p, eps, leg_id):
    """
    Trace la trajectoire de la patte du robot entre 2 points de l'espace en suivant normalized_move_xyz
    """
    # Théorique
    Xt = [x0, x]
    Yt = [y0, y]
    Zt = [z0, z]

    # Positions des vérins
    L = normalized_move_xyz(x, y, z, V[0], V[1], V[2], dstep, p, eps, leg_id)
    V1 = [v[0] for v in L]
    V2 = [v[1] for v in L]
    V3 = [v[2] for v in L]
    T = np.linspace(0, 1, len(L))

    # Positions du bout de la patte
    Pos = list(map(direct_xyz, [v[0] for v in L], [v[1] for v in L], [v[2] for v in L], [kin.FL for v in L]))
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
    ax.set_xbound(400, 800)
    ax.set_ybound(400, 800)
    ax.set_zbound(-300, -700)
    plt.tight_layout()
    plt.show()

    # Tracé des positions des vérins
    plt.plot(T, V1, label='V1' )
    plt.plot(T, V2, label='V2' )
    plt.plot(T, V3, label='V3' )
    plt.title("Elongations des vérins dans le mouvement")
    plt.show()

def test_A(dX, dZ, dalpha, v1, v2, v3):
  lpl = kin.LEGS[kin.FL]['lengths']
  # alpha = np.cos(v3_to_cos_angle(v3))
  # print (alpha)
  # pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
  # A = mat_A(pts, v1, v2, v3, alpha)
  #
  # dPos = np.array([dX/1000, dZ/1000, dalpha])
  # dV = A @ dPos
  #
  # new_v1, new_v2 = v1 + dV[0]*1000, v2 + dV[1]*1000
  # new_pts = kin.get_leg_points_V1_V2(new_v1 / 1000, new_v2 / 1000, lpl)
  #
  # new_alpha = np.arccos(v3_to_cos_angle(v3 + dV[2]*1000))
  # print (dV)
  # print(new_alpha)
  # err_rel = (abs(new_pts['J'][0] - pts['J'][0])*1000 - dX)/dX, (abs(new_pts['J'][1] - pts['J'][1])*1000 - dZ)/dZ, \
  #           (abs(new_alpha - alpha) - dalpha)/dalpha


  alpha = np.cos(v3_to_cos_angle(v3))
  pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
  dPos = np.array([dX/1000, dZ/1000, dalpha])

  # en utilisant mat A :
  A = mat_A(pts, v1, v2, v3, alpha)
  dVa = A @ dPos
  v1a = v1 + dVa[0]*1000
  v2a = v2 + dVa[1]*1000
  v3a = v3 + dVa[2]*1000
  print(direct_xyz(v1a, v2a, v3a, kin.FL))

  # en séparant jacob 2D et v3
  J = gen_jacob_2(pts, v1, v2)
  dPos2d = np.array([dPos[0], dPos[1]])
  dVj = J @ dPos2d
  v1j = v1 + dVj[0]*1000
  v2j = v2 + dVj[1]*1000
  v3j = cos_angle_to_v3(np.cos(alpha + dalpha))
  print(direct_xyz(v1j, v2j, v3j, kin.FL))


def test_jacob_2_direct():
    """
    Comparaison des resultats en J du modèle direct plan utilisant une jacobienne et de celui utilisant la géométrie
    """
    leg_id = kin.FL
    V = np.array([495, 585, 515])
    x0, y0, z0 = direct_xyz(V[0], V[1], V[2], leg_id)
    print("position de départ :", x0, y0, z0)
    deltaV = np.array([5, -5, 5])
    V += deltaV
    deltaV2 = np.array([deltaV[0], deltaV[1]])
    Jacob = gen_jacob_direct(kin.get_leg_points_V1_V2(V[0] / 1000, V[1] / 1000, kin.LEGS[leg_id]['lengths']), V[0] / 1000, V[1] / 1000)
    deltaX = Jacob @ deltaV2
    calpha = v3_to_cos_angle(V[2])
    deltax = np.array([
        deltaX[0] * np.sqrt(1-calpha**2),
        deltaX[0] * calpha,
        deltaX[1]
    ])
    print("position selon jacob_direct :", x0+deltax[0], y0+deltax[1], z0+deltax[2])
    print("position selon direct Julien :", direct_xyz(V[0], V[1], V[2], leg_id))
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
    pts = kin.get_leg_points_V1_V2(v1, v2, kin.LEGS[leg_id]['lengths'])
    x_E, y_E = pts['E']
    x_F, y_F = pts['F']
    x_G, y_G = pts['G']
    x_H, y_H = pts['H']
    x_I, y_I = pts['I']
    print('initial D :', pts['D'])
    pts_end = kin.get_leg_points_V1_V2((v1 + dv1), (v2 + dv2), kin.LEGS[leg_id]['lengths'])
    A = distance_3_points(pts['D'], pts['A'], pts['C'])
    B = np.array([0, 2 * v1])
    M_D = inv(A) @ B
    deltaD = M_D * dv1
    D = pts['D'] + deltaD
    print("real D :", pts_end['D'])
    print("jacob D :", D, "\n")

    M_E = (kin.LEGS[kin.FL]['lengths']['ae'] / (kin.LEGS[kin.FL]['lengths']['ae'] - kin.LEGS[kin.FL]['lengths']['de'])) * M_D
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

    M_G = ((kin.LEGS[kin.FL]['lengths']['ef'] + kin.LEGS[kin.FL]['lengths']['fg']) / kin.LEGS[kin.FL]['lengths']['ef']) * M_F
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

    M_J = (kin.LEGS[kin.FL]['lengths']['gj'] / kin.LEGS[kin.FL]['lengths']['gi']) * M_I
    deltaJ = M_J @ np.array([dv1, dv2])
    J = pts['J'] + deltaJ
    print('initial J:', pts['J'])
    print("real J:", pts_end['J'])
    print("jacob J", J, "\n")


def test_comparaison_minimize_vs_jacob_indirect(v1, v2, dx, dy):
    """
    test de comparaison des méthodes de cinématique indirecte
    """
    leg_id = kin.FL
    lpl = kin.LEGS[leg_id]['lengths']
    pts = kin.get_leg_points_V1_V2(v1, v2, lpl)
    x0, y0 = pts['J']
    dX = np.array([dx, dy])
    x, y = x0 + dx, y0 + dy
    print("initialement v1, v2 =", v1, v2, "_________position en x, y =", x0, y0, "_______cible x, y =", x, y)

    r = kin.inverse_kinetic(x, y, lpl)
    v1_minimize, v2_minimize = r[0], r[1]
    print("selon minimize, V =", r[0], r[1], "________________________________position atteinte en x, y =",
          kin.get_leg_points_V1_V2(r[0], r[1], lpl)['J'])

    J = gen_jacob_2(pts, v1, v2)
    dV = J @ dX
    print("selon la methode jacobienne, V =", v1 + dV[0], v2 + dV[1], "___________________position atteinte en x, y =",
          kin.get_leg_points_V1_V2(v1 + dV[0], v2 + dV[1], lpl)['J'])

    P = 2 * J.T @ J
    q = J.T @ (np.array([x0, y0]) - np.array([x, y]))
    lb = np.array([(0.450 - v1), (0.450 - v2)])
    ub = np.array([(0.650 - v1), (0.650 - v2)])
    dV = solve_qp(P, q, lb=lb, ub=ub)
    print("en appliquant un solveur sur la jacobienne, V =", v1 + dV[0], v2 + dV[1], "_____position atteinte en x, y =",
          kin.get_leg_points_V1_V2(v1 + dV[0], v2 + dV[1], lpl)['J'])
    print("\n")

def gen_jacob_direct(pts, v1, v2):
  """
  Retourne la Jacobienne correspondant au modèle cinématique indirect dans le plan de la patte
  Prend en argument la position des points de la patte et l'élongation des verrins en m

  >>> gen_jacob_2(kin.get_leg_points_V1_V2(0.495, 0.585, kin.LEGS[kin.FL]['lengths']), 0.495, 0.585) @ np.array([0, 0])
  array([0., 0.])
  """
  x_E, y_E = pts['E']
  x_F, y_F = pts['F']
  x_G, y_G = pts['G']
  x_H, y_H = pts['H']
  x_I, y_I = pts['I']

  A = distance_3_points(pts['D'], pts['A'], pts['C'])
  B = np.array([0, 2*v1])
  M_D = inv(A) @ B

  M_E = (kin.LEGS[kin.FL]['lengths']['ae']/(kin.LEGS[kin.FL]['lengths']['ae'] - kin.LEGS[kin.FL]['lengths']['de'])) * M_D

  A = distance_3_points(pts['F'], pts['E'], pts['B'])
  B = np.array([
    [2*(x_F - x_E), 2*(y_F - y_E)],
    [0, 0]])
  M_F = (inv(A) @ B) @ M_E

  M_G =((kin.LEGS[kin.FL]['lengths']['ef']+kin.LEGS[kin.FL]['lengths']['fg']) / kin.LEGS[kin.FL]['lengths']['ef']) * M_F

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

  return (kin.LEGS[kin.FL]['lengths']['gj']/kin.LEGS[kin.FL]['lengths']['gi']) * M_I

def make_a_penalty(v1, v2, d, eps, leg_id):
  lpl = kin.LEGS[leg_id]['lengths']
  pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
  X0, Z0 = pts['J'][0] * 1000, pts['J'][1] * 1000
  alpha = np.cos(np.pi / 4)
  res = []

  # drawing the shoot trajectoire
  Lx = np.zeros(d//eps + 1)
  Lz = np.zeros(d//eps + 1)
  for k in range(d//eps + 1):
    X = X0 + k * eps
    Z = Z0
    Lx[k], Lz[k] = X, Z
  print(Lx, Lz)
  # shooting that penalty
  for k in range(1, d//eps + 1):
    X, Z = Lx[k], Lz[k]
    print("POSITION ______actual :",X0, Z0,"__________cible :", X, Z)
    print("VERINS_________actual :", v1, v2)
    dX = np.array([X - X0, Z - Z0])
    J = gen_jacob_2(pts, v1 / 1000, v2 / 1000)
    P = 2 * J.T @ J
    q = J.T @ (np.array([X0 / 1000, Z0 / 1000]) - np.array([X / 1000, Z / 1000]))
    lb = np.array([(450.0 - v1), (450.0 - v2)]) / 1000
    ub = np.array([(650.0 - v1), (650.0 - v2)]) / 1000
    dV = J @ dX
    # print(dV)
    # dV = solve_qp(P, q, lb=lb, ub=ub)
    # print(dV)
    v1 += dV[0]
    v2 += dV[1]
    res.append((v1, v2))
    pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
    X0, Z0 = pts['J'][0] * 1000, pts['J'][1] * 1000
  return res

def test_circle_2(v1, v2, r, n, leg_id):
    """
    Trace la trajectoire de la patte du robot réalisant des petits cercles
    """
    # Positions des vérins
    lpl = kin.LEGS[kin.FL]['lengths']
    L = draw_circle_2(v1, v2, r, n, leg_id)
    V1 = [v[0] for v in L]
    V2 = [v[1] for v in L]
    T = np.linspace(0, 1, len(L))


    # Positions du bout de la patte
    Pos = list(map(lambda v1, v2, lpl: kin.get_leg_points_V1_V2(v1, v2, lpl)['J'], [v[0]/1000 for v in L], [v[1]/1000 for v in L], [lpl for v in L]))
    Xp = [p[0]*1000 for p in Pos]
    Yp = [p[1]*1000 for p in Pos]

    # Tracé de la position du bout de la patte
    fig = plt.figure()
    # Cercle théorique
    pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
    X0, Z0 = pts['J'][0] * 1000, pts['J'][1] * 1000
    Lx = np.zeros(n + 1)
    Lz = np.zeros(n + 1)
    for k in range(n + 1):
        X = X0 + r * np.cos(2 * k * np.pi / n) - r
        Z = Z0 + r * np.sin(2 * k * np.pi / n)
        Lx[k], Lz[k] = X, Z
    plt.plot(Lx, Lz)

    plt.scatter(Xp, Yp, label='Positions',
               marker='.',s=3, c='red')  # Tracé des points 3D
    plt.title("Trajectoire patte FL rotation circulaire dand le plan XZ")

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')
    plt.show()

    # Tracé des positions des vérins
    plt.plot(T, V1, label='V1' )
    plt.plot(T, V2, label='V2' )
    plt.title("Elongations des vérins dans le mouvement")
    # plt.plot.set_xlabel('L (en mm)')
    # plt.plot.set_ylabel('T')
    plt.show()

def draw_move_leg(traj, v1, v2, v3, leg_id, upgrade=False):
    """
        Trace la trajectoire de la patte du robot suivant traj avec move_leg
    """
    # Trajectoire
    Xt = [p[0] - 500 for p in traj]
    Yt = [p[1] - 500 for p in traj]
    Zt = [p[2] for p in traj]

    # Elongations des vérins
    Ver = move_leg(traj, v1, v2, v3, leg_id, upgrade=upgrade)
    V1 = [v[0] for v in Ver]
    V2 = [v[1] for v in Ver]
    V3 = [v[2] for v in Ver]
    T = np.linspace(0, 1, len(Ver))

    # Positions du bout de la patte
    Pos = list(map(direct_xyz, [v[0] for v in Ver], [v[1] for v in Ver], [v[2] for v in Ver], [kin.FL for v in Ver]))
    Xp = [p[0] for p in Pos]
    Yp = [p[1] for p in Pos]
    Zp = [p[2] for p in Pos]

    fig = plt.figure()
    ax = fig.gca(projection='3d')  # Affichage en 3D
    ax.plot(Xt, Yt, Zt, label='Théorique')  # Tracé de la courbe théorique
    ax.scatter(Xp, Yp, Zp, label='Positions', marker='.', s=3, c='red')  # Tracé des points réels
    plt.title("Trajectoire du bout de la patte dans l'espace")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xbound(0, 1000)
    ax.set_ybound(0, 1000)
    ax.set_zbound(0, -1000)
    plt.show()

    plt.plot(T, V1, label='V1' )
    plt.plot(T, V2, label='V2' )
    plt.plot(T, V3, label='V3' )
    plt.title("Elongations des vérins dans le mouvement")
    plt.show()

    ErrX ,ErrY, ErrZ = [], [], []
    for k in range(len(Ver)):
        ErrX.append(Xp[k]-Xt[k])
        ErrY.append(Yp[k]-Yt[k])
        ErrZ.append(Zp[k]-Zt[k])
    plt.plot(ErrX, label="x error")
    plt.plot(ErrY, label="y error")
    plt.plot(ErrZ, label="z error")
    plt.legend()
    plt.ylabel('error in the movement (mm)')
    plt.xlabel('steps')
    plt.show()


def test_circle_3(v1, v2, v3, r, n, leg_id):
    """
    Trace la trajectoire de la patte du robot réalisant des petits cercles
    """
    # Positions des vérins
    lpl = kin.LEGS[kin.FL]['lengths']
    L = draw_circle_3(v1, v2, v3, r, n, leg_id)
    V1 = [v[0] for v in L]
    V2 = [v[1] for v in L]
    V3 = [v[2] for v in L]
    T = np.linspace(0, 1, len(L))

    # Positions du bout de la patte
    Pos = list(map(direct_xyz, [v[0] for v in L], [v[1] for v in L], [v[2] for v in L], [kin.FL for v in L]))
    Xp = [p[0] for p in Pos]
    Yp = [p[1] for p in Pos]
    Zp = [p[2] for p in Pos]

    pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
    X0, Z0 = pts['J'][0] * 1000, pts['J'][1] * 1000
    x0, y0, z0 = d2_to_d3(X0, Z0, v3_to_cos_angle(v3))
    res = []

    # Calcul des points du cercle théorique
    Lx = []
    Ly = []
    Lz = []
    for k in range(n + 1):
        Lx.append(x0 + r * np.cos(2 * k * np.pi / n) - r)
        Ly.append(y0 + r * np.sin(2 * k * np.pi / n))
        Lz.append(z0)

    # Tracé de la position du bout de la patte
    fig = plt.figure()
    ax = fig.gca(projection='3d')  # Affichage en 3D
    ax.plot(Lx, Ly, Lz, label='Théorique')  # Tracé de la courbe 3D
    ax.scatter(Xp, Yp, Zp, label='Positions',
               marker='.', s=3, c='red')  # Tracé des points 3D
    plt.title("Trajectoire cercle dans l'espace")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xbound(200, 700)
    ax.set_ybound(400, 900)
    ax.set_zbound(-200, -700)
    plt.tight_layout()
    plt.show()

    # Tracé des positions des vérins
    plt.plot(T, V1, label='V1' )
    plt.plot(T, V2, label='V2' )
    plt.plot(T, V3, label='V3' )
    plt.title("Elongations des vérins dans le mouvement")
    plt.show()

def test_penalty_move_XZ(v1, v2, d, eps, leg_id):
    """
    Trace la trajectoire de la patte du robot réalisant des petits cercles
    """
    # Positions des vérins
    lpl = kin.LEGS[kin.FL]['lengths']
    L = make_a_penalty(v1, v2, d, eps, leg_id)
    V1 = [v[0] for v in L]
    V2 = [v[1] for v in L]
    T = np.linspace(0, 1, len(L))


    # Positions du bout de la patte
    Pos = list(map(lambda v1, v2, lpl: kin.get_leg_points_V1_V2(v1, v2, lpl)['J'], [v[0]/1000 for v in L], [v[1]/1000 for v in L], [lpl for v in L]))
    Xp = [p[0]*1000 for p in Pos]
    Yp = [p[1]*1000 for p in Pos]

    # Tracé de la position du bout de la patte
    fig = plt.figure()
    # drawing the shoot trajectoire
    pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
    X0, Z0 = pts['J'][0] * 1000, pts['J'][1] * 1000
    Lx = np.zeros(d // eps + 1)
    Lz = np.zeros(d // eps + 1)
    for k in range(d // eps + 1):
        X = X0 + k * eps
        Z = Z0
        Lx[k], Lz[k] = X, Z

    plt.plot(Lx, Lz)

    plt.scatter(Xp, Yp, label='Positions',
                marker='.', c = 'red', s = 3)  # Tracé des points 3D
    plt.title("Trajectoire rectiligne de la patte FL dans le plan XZ")

    plt.xlabel('X')

    plt.ylabel('Y')
    plt.axis('equal')
    plt.show()

    # Tracé des positions des vérins
    plt.plot(T, V1, label='V1' )
    plt.plot(T, V2, label='V2' )
    plt.title("Elongations des vérins dans le mouvement")
    # plt.plot.set_xlabel('L (en mm)')
    # plt.plot.set_ylabel('T')
    plt.show()

################################# TESTS ####################################

t_move = 0
# test de normalized_move_xyz
if t_move:
    x0, y0, z0 = direct_xyz(V[0], V[1], V[2], kin.FL)
    print(x0, y0, z0)
    test_normalized_move_xyz(x0, y0, z0, 750, y0, z0, 1, 1, 5, kin.FL)
    test_normalized_move_xyz(x0, y0, z0, x0, 700, z0, 1, 1, 5, kin.FL)
    test_normalized_move_xyz(x0, y0, z0, x0, y0, -550, 1, 1, 5, kin.FL)
    test_normalized_move_xyz(x0, y0, z0, 750, 600, -400, 1, 1, 5, kin.FL)

t_jacob = 0
# test de l'erreur relative de position en appliquant la Jacobienne
if t_jacob:
    test_jacob_2(495, 585, 0, 0)
    test_jacob_2(505, 585, 0, 0)
    test_jacob_2(515, 585, 0, 0)

    test_jacob_2(495, 585, 10, 0)
    test_jacob_2(505, 585, 10, 0)
    test_jacob_2(515, 585, 10, 0)

    # test_jacob_2(495, 585, -10, 0)
    # test_jacob_2(505, 585, -10, 0)
    # test_jacob_2(515, 585, -10, 0)

    test_jacob_2(495, 585, 0, 10)
    test_jacob_2(505, 585, 0, 10)
    test_jacob_2(515, 585, 0, 10)

    # test_jacob_2(495, 585, 0, -10)
    # test_jacob_2(505, 585, 0, -10)
    # test_jacob_2(515, 585, 0, -10)

    # test_jacob_2(495, 585, 10, 10)
    # test_jacob_2(505, 585, 10, 10)
    # test_jacob_2(515, 585, 10, 10)

    # test_jacob_2(495, 585, 10, -10)
    # test_jacob_2(505, 585, 10, -10)
    # test_jacob_2(515, 585, 10, -10)

    # test_jacob_2(495, 585, -10, 10)
    # test_jacob_2(505, 585, -10, 10)
    # test_jacob_2(515, 585, -10, 10)

    # test_jacob_2(495, 585, -10, -10)
    # test_jacob_2(505, 585, -10, -10)
    # test_jacob_2(515, 585, -10, -10)

t_a = 0
if t_a:
    print (test_A(1, -1, 0.1, V[0], V[1], V[2]))

t_jacob_direct = 0
if t_jacob_direct:
    test_jacob_2_direct()
    test_precision_jacobienne_en_direct(0.495, 0.555, 0.015, -0.015, kin.FL)

test_comp_indirect = 0
if test_comp_indirect:
    test_comparaison_minimize_vs_jacob_indirect(0.485, 0.565, 0.001, -0.0015)
    test_comparaison_minimize_vs_jacob_indirect(0.485, 0.565, 0.01, -0.015)
    test_comparaison_minimize_vs_jacob_indirect(0.485, 0.565, 0.1, -0.15)
    test_comparaison_minimize_vs_jacob_indirect(0.485, 0.565, -0.01, +0.015)

t_different_moves = 1
if t_different_moves:
    draw_move_leg(draw_circle(200, 20, 550, 600, 515, kin.FL), 550, 600, 515, kin.FL, upgrade=True)
    draw_move_leg(draw_circle(200, 20, 550, 600, 515, kin.FL), 550, 600, 515, kin.FL, upgrade=False)
    # test_circle_2(450, 500, 200, 200, kin.FL)
    # test_circle_3(550, 600, 515, 200, 200, kin.FL)
    # test_penalty_move_XZ(450, 500, 500, 10, kin.FL)

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
        kin.get_leg_points_V1_V2(495 / 1000, 585 / 1000, kin.LEGS[kin.FR]['lengths'])['J']
    t2 = time.time() - t

    print("direct_v1 prend ", t1 * 100, " us")
    print("L'algo de Julien prend ", t2 * 100, " us")

    print("direct v1 retourne : ", direct_v1(495, 585))
    print("L'algo de Julien retourne :", kin.get_leg_points_V1_V2(495 / 1000, 585 / 1000, kin.LEGS[kin.FR]['lengths'])['J'])

############################################################################

