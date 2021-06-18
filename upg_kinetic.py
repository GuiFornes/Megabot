# L'ensemble des distances sont exprimées en mm : segments de patte et élongations des vérins

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import inv
from qpsolvers import solve_qp
import kinetic as kin

x_A = - kin.LEGS[kin.FL]['lengths']['ao']
y_A = 0
x_B = kin.LEGS[kin.FL]['lengths']['bo']
y_B = 0

AB = kin.LEGS[kin.FL]['lengths']['ao'] + kin.LEGS[kin.FL]['lengths']['bo']
AC = np.sqrt((AB + kin.LEGS[kin.FL]['lengths']['bcx']) ** 2 + kin.LEGS[kin.FL]['lengths']['bcy'] ** 2)
BC = np.sqrt(kin.LEGS[kin.FL]['lengths']['bcx'] ** 2 + kin.LEGS[kin.FL]['lengths']['bcy'] ** 2)
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
LO = np.sqrt(LM ** 2 + MO ** 2)

# Matrices de rotation permettant le passage d'une patte à l'autre
MR = np.array([[[1, 0, 0],
                [0, 1, 0],
                [0, 0, 1]],

               [[0, -1, 0],
                [1, 0, 0],
                [0, 0, 1]],

               [[0, 1, 0],
                [-1, 0, 0],
                [0, 0, 1]],

               [[-1, 0, 0],
                [0, -1, 0],
                [0, 0, 1]]])

# Décalage entre le centre du robot et le point O de la jambe FL
L = np.array([500, 500, 0])

################################# 2D #######################################

def distance_3_points(A, B, C):
    """
    Fonction auxiliaire de gen_jacob_2
    Retourne la matrice des distances de A à B et à C
    """
    return 2 * np.array([[A[0] - B[0], A[1] - B[1]], [A[0] - C[0], A[1] - C[1]]])


def gen_jacob_2(pts, v1, v2):
    """
    Retourne la jacobienne correspondant au modèle cinématique indirect dans le plan de la patte
    Prend en argument la position des points de la patte et l'élongation des verrins en m
    La jacobienne doit être appliqué sur des élongations en m et retourne des position en m

    >>> gen_jacob_2(kin.get_leg_points_V1_V2(0.495, 0.585, kin.LEGS[kin.FL]['lengths']), 0.495, 0.585) @ np.array([0, 0])
    array([0., 0.])

    # >>> gen_jacob_2(kin.get_leg_points_V1_V2(0.495, 0.585, kin.LEGS[kin.FL]['lengths']), 0.495, 0.585) @ np.array([1, 0])
    # array([1, 0])
    #
    # >>> gen_jacob_2(kin.get_leg_points_V1_V2(0.495, 0.585, kin.LEGS[kin.FL]['lengths']), 0.495, 0.585) @ np.array([1, 0])

    """
    x_E, y_E = pts['E']
    x_F, y_F = pts['F']
    x_G, y_G = pts['G']
    x_H, y_H = pts['H']
    x_I, y_I = pts['I']

    A = distance_3_points(pts['D'], pts['A'], pts['C'])
    B = np.array([0, 2 * v1])
    M_D = inv(A) @ B

    M_E = (kin.LEGS[kin.FL]['lengths']['ae'] / (
                kin.LEGS[kin.FL]['lengths']['ae'] - kin.LEGS[kin.FL]['lengths']['de'])) * M_D

    A = distance_3_points(pts['F'], pts['E'], pts['B'])
    B = np.array([
        [2 * (x_F - x_E), 2 * (y_F - y_E)],
        [0, 0]])
    M_F = (inv(A) @ B) @ M_E

    M_G = ((kin.LEGS[kin.FL]['lengths']['ef'] + kin.LEGS[kin.FL]['lengths']['fg']) / kin.LEGS[kin.FL]['lengths'][
        'ef']) * M_F

    M_H = (BH / BF) * M_F

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

    Jacob = inv((kin.LEGS[kin.FL]['lengths']['gj'] / kin.LEGS[kin.FL]['lengths']['gi']) * M_I)

    return Jacob


def d3_to_d2(x, y, z):
    """
    Retourne (X, Z, alpha) les coordonnées cylindriques du bout de la patte
    """
    X = np.sqrt(x ** 2 + y ** 2)
    Z = z
    calpha = y / X
    return X, Z, calpha


def d2_to_d3(X, Z, calpha):
    """
    Retourne (x, y, z) les coordonnées cartésiennes du bout de la patte
    """
    x = X * np.sqrt(1 - calpha ** 2)
    y = X * calpha
    z = Z
    return x, y, z


def solve_indirect_cyl(x, y, z, x0, y0, z0, v1, v2, v3, lpl, pts):
    """
    Retourne les élongations (v1, v2, v3) permettant de placer le bout de la patte en (x, y, z)
    en minimisant les erreurs d'élongation sur v1 et v2 (utilise Jacob_2)
    Toutes les longueurs sont en mm (entrées comme sorties)
    """
    X, Z, calpha = d3_to_d2(x, y, z)
    Xt = np.array([X, Z]) / 1000
    new_v3 = cos_angle_to_v3(calpha, lpl)

    X, Z, calpha = d3_to_d2(x0, y0, z0)
    X0 = np.array([X, Z]) / 1000

    J = gen_jacob_2(pts, v1 / 1000, v2 / 1000)

    P = 2 * J.T @ J
    q = J.T @ (X0 - Xt)
    lb = np.array([(450.0 - v1) / 1000, (450.0 - v2) / 1000])
    ub = np.array([(650.0 - v1) / 1000, (650.0 - v2) / 1000])
    dV = solve_qp(P, q, lb=lb, ub=ub)

    return v1 + dV[0] * 1000, v2 + dV[1] * 1000, new_v3


################################# 3D #######################################

def mat_A(pts, v1, v2, v3, alpha):
    """
    Fonction auxiliaire de gen_jacob_3
    Génère la matrice A conformément à nos équations (cf. indirect.pdf)
    Toutes les longueurs en m
    """
    Jacob = gen_jacob_2(pts, v1, v2)

    A = np.array([
        [Jacob[0][0], Jacob[0][1], 0],
        [Jacob[1][0], Jacob[1][1], 0],
        [0, 0, np.sin(alpha - np.arccos(MO / LO)) * (KO / 1000) * (LO / 1000) / v3]])

    return A


def mat_B(pts, alpha):
    """
    Fonction auxiliaire de gen_jacob_3
    Génère la matrice B conformément à nos équations (cf. indirect.pdf)
    Toutes les longueurs en m
    """
    X = pts['J'][0]

    B = np.array([
        [np.cos(np.pi / 2 - alpha), 0, X * np.sin(np.pi / 2 - alpha)],
        [np.sin(np.pi / 2 - alpha), 0, -X * np.cos(np.pi / 2 - alpha)],
        [0, 1, 0]])

    return B


def gen_jacob_3(v1, v2, v3, alpha, lpl):
    """
    Retourne la jacobienne correspondant au modèle cinématique indirect dans le repère cartésien centré en O
    Prend en argument l'élongation des verrins en m et l'angle alpha en radian
    La jacobienne doit être appliquée sur des élongations en m et retourne des position en m
    """
    pts = kin.get_leg_points_V1_V2(v1, v2, lpl)
    A = mat_A(pts, v1, v2, v3, alpha)
    B = mat_B(pts, alpha)

    return A @ inv(B)


def solve_indirect_cart(x, y, z, x0, y0, z0, v1, v2, v3, lpl, pts):
    """
    Retourne les élongations (v1, v2, v3) permettant de placer le bout de la patte en (x, y, z)
    en minimisant les erreurs d'élongation sur v1, v2 et v3 (utilise Jacob_3)
    Toutes les longueurs sont en mm (entrées comme sorties)
    """
    Xt = np.array([x, y, z]) / 1000
    X0 = np.array([x0, y0, z0]) / 1000

    J = gen_jacob_3(v1 / 1000, v2 / 1000, v3 / 1000, np.arctan(x0 / y0), lpl)

    P = 2 * J.T @ J
    q = J.T @ (X0 - Xt)
    lb = np.array([450.0 - v1, 450.0 - v2, 450.0 - v3]) / 1000
    ub = np.array([650.0 - v1, 650.0 - v2, 650.0 - v3]) / 1000
    dV = solve_qp(P, q, lb=lb, ub=ub)

    return v1 + dV[0] * 1000, v2 + dV[1] * 1000, v3 + dV[2] * 1000


################################ MOVE ######################################

def gen_jacob_12(V):
    J_12 = []
    for i in range(4):
        # Calcul de la jacobienne de la patte
        v1, v2, v3 = V[i*3], V[i*3 + 1], V[i*3 + 2]
        lpl = kin.LEGS[i]['lengths']
        alpha = np.arccos(v3_to_cos_angle(v3, lpl))
        J_3 = gen_jacob_3(v1 / 1000, v2 / 1000, v3 / 1000, alpha, lpl) @ MR[i].T

        # Ajout à J_12
        L = np.zeros((3, 12))
        for j in range(3):
            for k in range(3):
                L[j][i*3 + k] = J_3[j][k]
        J_12.append(L[0])
        J_12.append(L[1])
        J_12.append(L[2])
    return J_12

def direct_12(V):
    """
    Retourne les positions des extrémités des 4 pattes correspondant aux élongations V des vérins
    Les distances sont exprimées en mm et les coordonnées sont exprimées dans le référentiel du robot
    """
    R = []
    for i in range(4):
        R = np.append(R, direct_robot(V[i*3], V[i*3 + 1], V[i*3 + 2], i))
    return R

def move_12(traj, V):
    """
    Retourne le tableau des élongations successives des 12 vérins permettant aux extrémités des 4 pattes de suivre les trajectoires qui leur ont été attribuées par traj
    traj : liste des positions successives des extrémités des 4 pattes sous la forme [[FL_x, FL_y, FL_z, FR_x, FR_y, FR_z, RL_x, RL_y, RL_z, RR_x, RR_y, RR_z], ...]
    V : liste des élongations initiales des 12 vérins (dans l'ordre FL, FR, RL, RR) sous la forme [v1, v2, v3, v4, ..., v12]
    V doit correspondre à la première position de traj
    Toutes les longueurs sont en mm, les coordonnées des trajectoires sont exprimées dans le référentiel du robot
    """
    V0 = V
    R = [V0]
    for i in range(1, len(traj)):
        X0 = direct_12(V)
        dX = traj[i] - X0
        J = gen_jacob_12(V)
        dV = J @ dX
        V0 = V0 + dV
        #for v in V0: assert 450 < v < 650, 'Elongation de vérin invalide'
        R.append(V0)
    return R

def draw_circle_12(n, r, V):
    traj_FL = draw_circle(n, r, V[0], V[1], V[2], 0)
    traj_FR = draw_circle(n, r, V[3], V[4], V[5], 1)
    traj_RL = draw_circle(n, r, V[6], V[7], V[8], 2)
    traj_RR = draw_circle(n, r, V[9], V[10], V[11], 3)
    traj = []
    for i in range(n):
        t = traj_FL[i]
        t = np.append(t, MR[1] @ traj_FR[i])
        t = np.append(t, MR[2] @ traj_RL[i])
        t = np.append(t, MR[3] @ traj_RR[i])
        traj.append(t)
    return traj


def move_4_legs(traj, V, upgrade=False):
    """
    Retourne le tableau des élongations successives des 12 vérins (3 par 3) permettant aux
    extrémités des 4 pattes de suivre les trajectoires qui leur ont été attribuées par traj
    traj : [traj_0, traj_1, traj_2, traj_3], liste des 4 trajectoires suivies par les 4 pattes avec traj_i = [[x0, y0, z0], [x1, y1, z1], ...]
    V : [V_0, V_1, V_2, V_3], liste des élongations initiales des 3 vérins dans chacune des 4 pattes avec V_i = [v1, v2, v3]
    Toutes les longueurs sont en mm, les coordonnées des trajectoires sont exprimées dans le repère du robot
    """
    # Calcul des élongations de chacunes des pattes
    Ver = []
    for i in range(4):
        Ver.append(move_leg(traj[i], V[i][0], V[i][1], V[i][2], i, upgrade=upgrade))
    # Mise sous le format attendu
    R = []
    for k in range(len(Ver[0])):
        r = [Ver[0][k],
             Ver[1][k],
             Ver[2][k],
             Ver[3][k]]
        R.append(r)
    return R


def move_leg(traj, v1, v2, v3, leg_id, display=False, upgrade=False):
    """
    Retourne la liste des élongations des vérins permettant au bout de la patte de suivre traj
    Prend en argument la succession de positions formant traj, les élongations initiales des vérins et l'id de la patte
    Les élongations initiales des vérins doivent placer la patte au premier point de traj
    Toutes les longueurs sont en mm (entrée comme sortie)
    """
    R = [(v1, v2, v3)]
    err = 0
    if upgrade : prev_T = MR[leg_id].T @ traj[0] - L

    # Parcours de traj
    for i in range(1, len(traj)):
        lpl = kin.LEGS[leg_id]['lengths']
        pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
        X, Z = pts['J'][0] * 1000, pts['J'][1] * 1000
        x0, y0, z0 = d2_to_d3(X, Z, v3_to_cos_angle(v3, lpl))

        if display:
            print("POSITIONS ______actual :", x0, y0, z0, "__________target :", traj[i][0], traj[i][1], traj[i][2])
            print("VERINS_________actual :", v1, v2, v3)

        T = MR[leg_id].T @ traj[i] - L
        dX = np.array([T[0] - x0, T[1] - y0, T[2] - z0])
        if upgrade:
            dX += err
            err = np.array([prev_T[0] - x0, prev_T[1] - y0, prev_T[2] - z0])
            prev_T = T
        J = gen_jacob_3(v1 / 1000, v2 / 1000, v3 / 1000, np.arccos(v3_to_cos_angle(v3, lpl)), lpl)
        dV = J @ dX
        v1, v2, v3 = v1 + dV[0], v2 + dV[1], v3 + dV[2]
        R.append((v1, v2, v3))
    return R


def draw_circle(r, n, v1, v2, v3, leg_id):
    lpl = kin.LEGS[leg_id]['lengths']
    pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
    X0, Z0 = pts['J'][0] * 1000, pts['J'][1] * 1000
    x0, y0, z0 = d2_to_d3(X0, Z0, v3_to_cos_angle(v3, lpl))
    R = []
    for k in range(n + 1):
        R.append(np.array([x0+ r * np.cos(2 * k * np.pi / n) - r,
                           y0 + r * np.sin(2 * k * np.pi / n),
                           z0]) + L)
    return R


def normalized_move_xyz(x, y, z, v1, v2, v3, dstep, p, eps, leg_id):
    """
    Retourne la liste des élongations successives des vérins permettant de placer le bout de la patte en (x,y,z)
    Effectue une normalisation du déplacement
    """
    L = []
    c = 0

    pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, kin.LEGS[kin.FL]['lengths'])
    X, Z = pts['J'][0] * 1000, pts['J'][1] * 1000
    x0, y0, z0 = d2_to_d3(X, Z, v3_to_cos_angle(v3))
    dist = distance(x - x0, y - y0, z - z0)
    while dist > eps and c < 300:
        c += 1
        U = np.array([(x - x0), (y - y0), (z - z0)])
        U = dstep / 100 * U  # / np.linalg.norm(U)**p
        dx, dy, dz = U[0], U[1], U[2]

        v1, v2, v3 = solve_indirect_cart(x0 + dx, y0 + dy, z0 + dz, x0, y0, z0, v1, v2, v3, leg_id, pts)
        L.append((v1, v2, v3))

        pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, kin.LEGS[kin.FL]['lengths'])
        X, Z = pts['J'][0] * 1000, pts['J'][1] * 1000
        x1, y1, z1 = x0 + dx, y0 + dy, z0 + dz
        x0, y0, z0 = d2_to_d3(X, Z, v3_to_cos_angle(v3))

        dist = distance(x - x0, y - y0, z - z0)
    return L


def draw_circle_2(v1, v2, r, n, leg_id, solved=False):
    lpl = kin.LEGS[leg_id]['lengths']
    pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
    X0, Z0 = pts['J'][0] * 1000, pts['J'][1] * 1000
    alpha = np.cos(np.pi / 4)
    res = []

    # Calcul des points du cercle
    Lx = []
    Lz = []
    for k in range(n + 1):
        Lx.append(X0 + r * np.cos(2 * k * np.pi / n) - r)
        Lz.append(Z0 + r * np.sin(2 * k * np.pi / n))

    # Parcours du cercle
    for k in range(1, n + 1):
        # print("POSITIONS ______actual :", X0, Z0,"__________cible :", Lx[k], Lz[k])
        # print("VERINS_________actual :", v1, v2)
        dX = np.array([Lx[k] - X0, Lz[k] - Z0])
        J = gen_jacob_2(pts, v1 / 1000, v2 / 1000)

        if solved:  # Utilisation du solveur
            P = 2 * J.T @ J
            q = J.T @ (np.array([X0 / 1000, Z0 / 1000]) - np.array([Lx[k] / 1000, Lz[k] / 1000]))
            lb = np.array([(450.0 - v1), (450.0 - v2)]) / 1000
            ub = np.array([(650.0 - v1), (650.0 - v2)]) / 1000
            dV = solve_qp(P, q, lb=lb, ub=ub) * 1000
        else:  # Utilisation de la jacobienne sans solveur
            dV = J @ dX

        v1 += dV[0]
        v2 += dV[1]
        res.append((v1, v2))
        pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
        X0, Z0 = pts['J'][0] * 1000, pts['J'][1] * 1000
    return res


def draw_circle_3(v1, v2, v3, r, n, leg_id, solved=False):
    lpl = kin.LEGS[leg_id]['lengths']
    pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
    X0, Z0 = pts['J'][0] * 1000, pts['J'][1] * 1000
    x0, y0, z0 = d2_to_d3(X0, Z0, v3_to_cos_angle(v3))
    res = []

    # Calcul des points du cercle
    Lx = []
    Ly = []
    Lz = []
    for k in range(n + 1):
        Lx.append(x0 + r * np.cos(2 * k * np.pi / n) - r)
        Ly.append(y0 + r * np.sin(2 * k * np.pi / n))
        Lz.append(z0)

    # Parcours du cercle
    for k in range(1, n + 1):
        # print("POSITIONS ______actual :",x0, y0, z0,"__________cible :", Lx[k], Ly[k], Lz[k])
        # print("VERINS_________actual :", v1, v2, v3)
        dX = np.array([Lx[k] - x0, Ly[k] - y0, Lz[k] - z0])
        alpha0 = np.arccos(v3_to_cos_angle(v3))
        J = gen_jacob_3(v1 / 1000, v2 / 1000, v3 / 1000, alpha0, kin.FL)

        if solved:  # Utilisation du solveur
            dV = np.array([0, 0, 0])
            # P = 2 * J.T @ J
            # q = J.T @ (np.array([X0/1000, Z0/1000]) - np.array([X/1000, Z/1000]))
            # lb = np.array([(450.0 - v1), (450.0 - v2)])/1000
            # ub = np.array([(650.0 - v1), (650.0 - v2)])/1000
            # dV = solve_qp(P, q, lb=lb, ub=ub) * 1000
        else:  # Utilisation de la jacobienne sans solveur
            dV = J @ dX

        v1 += dV[0]
        v2 += dV[1]
        v3 += dV[2]
        res.append((v1, v2, v3))
        pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
        X0, Z0 = pts['J'][0] * 1000, pts['J'][1] * 1000
        x0, y0, z0 = d2_to_d3(X0, Z0, v3_to_cos_angle(v3))
    return res


def direct_leg(v1, v2, v3):
    """
    Retourne les positions x, y, z du bout de la patte en fonctions de v1, v2, v3 dans le référentiel de la patte FL
    Se base sur le modèle direct de Julien
    """
    lpl = kin.LEGS[kin.FL]['lengths']
    X, Z = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)['J']
    calpha = v3_to_cos_angle(v3, lpl)
    x = X * np.sqrt(1 - calpha ** 2) * 1000
    y = X * calpha * 1000
    z = Z * 1000
    return x, y, z


def direct_robot(v1, v2, v3, leg_id):
    """
    Retourne les positions x, y, z du bout de la patte en fonctions de v1, v2, v3 dans le référentiel du robot
    Se base sur le modèle direct de Julien
    """
    lpl = kin.LEGS[leg_id]['lengths']
    X, Z = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)['J']
    calpha = v3_to_cos_angle(v3, lpl)
    Pos = np.array([X * np.sqrt(1 - calpha ** 2) * 1000, X * calpha * 1000, Z * 1000])
    return MR[leg_id] @ (Pos + L)

################################# OUTILS ###################################

def cos_angle_to_v3(cangle, lpl):
    """
    Fonction auxiliaire de normalized_move_xyz
    Retourne l'élongation de v3 en mm en fonction du cosinus de l'angle de la patte au chassis

    >>> cos_angle_to_v3(np.cos(np.pi/4), kin.LEGS[0]['lengths']) - al_kashi_longueur(KO, LO, np.pi/4 - np.arccos(MO/LO))
    0.0
    >>> v3_to_cos_angle(cos_angle_to_v3(np.cos(np.pi/4), kin.LEGS[0]['lengths']), kin.LEGS[0]['lengths']) - np.cos(np.pi/4) < 0.0000001
    True
    """
    KO = lpl['yaw_c']
    LM = lpl['yaw_b']
    MO = lpl['yaw_a']
    LO = np.sqrt(LM ** 2 + MO ** 2)
    return np.sqrt(
        LO ** 2 + KO ** 2 - 2 * LO * KO * (cangle * MO / LO + np.sqrt(1 - cangle ** 2) * np.sqrt(1 - (MO / LO) ** 2)))


def v3_to_cos_angle(v3, lpl):
    """
    Fonction auxiliaire de normalized_move_xyz
    Retourne le cosinus de l'angle de la patte au chassis en fonction de l'élongation de v3 en mm

    >>> v3_to_cos_angle(500, kin.LEGS[0]['lengths']) - np.cos(al_kashi_angle(LO, KO, 500) + np.arccos(MO/LO))
    0.0
    >>> cos_angle_to_v3(v3_to_cos_angle(500, kin.LEGS[0]['lengths']), kin.LEGS[0]['lengths']) - 500 < 0.0000001
    True
    """
    KO = lpl['yaw_c']
    LM = lpl['yaw_b']
    MO = lpl['yaw_a']
    LO = np.sqrt(LM ** 2 + MO ** 2)
    return (KO ** 2 + LO ** 2 - v3 ** 2) / (2 * KO * LO) * MO / LO - np.sqrt(
        1 - ((KO ** 2 + LO ** 2 - v3 ** 2) / (2 * KO * LO)) ** 2) * np.sqrt(1 - (MO / LO) ** 2)


def distance(x, y, z=0):
    """
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
    """
    return np.sqrt(x ** 2 + y ** 2 + z ** 2)


def al_kashi_longueur(a, b, alpha):
    return np.sqrt(a ** 2 + b ** 2 - 2 * a * b * np.cos(alpha))


def al_kashi_angle(a, b, c):
    return np.arccos((a ** 2 + b ** 2 - c ** 2) / (2 * a * b))


############################################################################

if __name__ == "__main__":
    import doctest

    doctest.testmod()
