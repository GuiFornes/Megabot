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

ori = np.array([[1, 1, -1, -1], [1, -1, -1, 1]])  # [[oritentation selon x][orientation selon y]]

V = 460, 565, 500


################################# 2D #######################################

def distance_3_points(A, B, C):
    '''
    Fonction auxiliaire de gen_jacob_2
    Retourne la matrice des distances de A à B et à C
    '''
    return 2 * np.array([[A[0] - B[0], A[1] - B[1]], [A[0] - C[0], A[1] - C[1]]])


def gen_jacob_2(pts, v1, v2):
    '''
    Retourne la jacobienne correspondant au modèle cinématique indirect dans le plan de la patte
    Prend en argument la position des points de la patte et l'élongation des verrins en m
    La jacobienne doit être appliqué sur des élongations en m et retourne des position en m

    >>> gen_jacob_2(kin.get_leg_points_V1_V2(0.495, 0.585, kin.LEGS[kin.FL]['lengths']), 0.495, 0.585) @ np.array([0, 0])
    array([0., 0.])

    # >>> gen_jacob_2(kin.get_leg_points_V1_V2(0.495, 0.585, kin.LEGS[kin.FL]['lengths']), 0.495, 0.585) @ np.array([1, 0])
    # array([1, 0])
    #
    # >>> gen_jacob_2(kin.get_leg_points_V1_V2(0.495, 0.585, kin.LEGS[kin.FL]['lengths']), 0.495, 0.585) @ np.array([1, 0])

    '''
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
    '''
    Retourne (X, Z, alpha) les coordonnées cylindriques du bout de la patte
    '''
    X = np.sqrt(x ** 2 + y ** 2)
    Z = z
    calpha = y / X
    return X, Z, calpha


def d2_to_d3(X, Z, calpha):
    '''
    Retourne (x, y, z) les coordonnées cartésiennes du bout de la patte
    '''
    x = X * np.sqrt(1 - calpha ** 2)
    y = X * calpha
    z = Z
    return x, y, z


def solve_indirect_cyl(x, y, z, x0, y0, z0, v1, v2, v3, leg_id, pts):
    '''
    Retourne les élongations (v1, v2, v3) permettant de placer le bout de la patte en (x, y, z)
    en minimisant les erreurs d'élongation sur v1 et v2 (utilise Jacob_2)
    Toutes les longueurs sont en mm (entrées comme sorties)
    '''
    X, Z, calpha = d3_to_d2(x, y, z)
    Xt = np.array([X, Z]) / 1000
    new_v3 = cos_angle_to_v3(calpha)

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
    '''
    Fonction auxiliaire de gen_jacob_3
    Génère la matrice A conformément à nos équations (cf. indirect.pdf)
    Toutes les longueurs en m
    '''
    Jacob = gen_jacob_2(pts, v1, v2)

    A = np.array([
        [Jacob[0][0], Jacob[0][1], 0],
        [Jacob[1][0], Jacob[1][1], 0],
        [0, 0, np.sin(alpha - np.arccos(MO / LO)) * (KO / 1000) * (LO / 1000) / v3]])

    return A


def mat_B(pts, alpha, leg_id):
    '''
    Fonction auxiliaire de gen_jacob_3
    Génère la matrice B conformément à nos équations (cf. indirect.pdf)
    Toutes les longueurs en m
    '''
    X = pts['J'][0]

    B = np.array([
        [np.cos(np.pi / 2 - alpha), 0, X * np.sin(np.pi / 2 - alpha)],
        [np.sin(np.pi / 2 - alpha), 0, -X * np.cos(np.pi / 2 - alpha)],
        [0, 1, 0]])

    return B


def gen_jacob_3(v1, v2, v3, alpha, leg_id):
    '''
    Retourne la jacobienne correspondant au modèle cinématique indirect dans le repère cartésien centré en O
    Prend en argument l'élongation des verrins en m et l'angle alpha en radian
    La jacobienne doit être appliquée sur des élongations en m et retourne des position en m
    '''
    pts = kin.get_leg_points_V1_V2(v1, v2, kin.LEGS[kin.FL]['lengths'])
    A = mat_A(pts, v1, v2, v3, alpha)
    B = mat_B(pts, alpha, leg_id)

    return A @ inv(B)


def solve_indirect_cart(x, y, z, x0, y0, z0, v1, v2, v3, leg_id, pts):
    '''
    Retourne les élongations (v1, v2, v3) permettant de placer le bout de la patte en (x, y, z)
    en minimisant les erreurs d'élongation sur v1, v2 et v3 (utilise Jacob_3)
    Toutes les longueurs sont en mm (entrées comme sorties)
    '''
    Xt = np.array([x, y, z]) / 1000
    X0 = np.array([x0, y0, z0]) / 1000

    J = gen_jacob_3(v1 / 1000, v2 / 1000, v3 / 1000, np.arctan(x0 / y0), kin.FL)

    P = 2 * J.T @ J
    q = J.T @ (X0 - Xt)
    lb = np.array([450.0 - v1, 450.0 - v2, 450.0 - v3]) / 1000
    ub = np.array([650.0 - v1, 650.0 - v2, 650.0 - v3]) / 1000
    dV = solve_qp(P, q, lb=lb, ub=ub)

    return v1 + dV[0] * 1000, v2 + dV[1] * 1000, v3 + dV[2] * 1000


################################ MOVE ######################################

def normalized_move_xyz(x, y, z, v1, v2, v3, dstep, p, eps, leg_id):
    '''
    Retourne la liste des élongations successives des verrins permettant de placer le bout de la patte en (x,y,z)
    Effectue une normalisation du déplacement
    '''
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


def direct_xyz(v1, v2, v3, leg_id):
    '''
    Retourne les positions x, y, z du bout de la patte en fonctions de v1, v2, v3
    Se base sur le modèle direct de Julien
    '''
    X, Z = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, kin.LEGS[kin.FL]['lengths'])['J']
    calpha = v3_to_cos_angle(v3)
    x = ori[leg_id][0] * X * np.sqrt(1 - calpha ** 2) * 1000
    y = ori[leg_id][1] * X * calpha * 1000
    z = Z * 1000
    return x, y, z


################################# OUTILS ###################################

def cos_angle_to_v3(cangle):
    '''
    Fonction auxiliaire de normalized_move_xyz
    Retourne l'élongation de v3 en mm en fonction du cosinus de l'angle de la patte au chassis

    >>> cos_angle_to_v3(np.cos(np.pi/4)) - al_kashi_longueur(KO, LO, np.pi/4 - np.arccos(MO/LO))
    0.0
    >>> v3_to_cos_angle(cos_angle_to_v3(np.cos(np.pi/4))) - np.cos(np.pi/4) < 0.0000001
    True
    '''
    return np.sqrt(
        LO ** 2 + KO ** 2 - 2 * LO * KO * (cangle * MO / LO + np.sqrt(1 - cangle ** 2) * np.sqrt(1 - (MO / LO) ** 2)))


def v3_to_cos_angle(v3):
    '''
    Fonction auxiliaire de normalized_move_xyz
    Retourne le cosinus de l'angle de la patte au chassis en fonction de l'élongation de v3 en mm

    >>> v3_to_cos_angle(500) - np.cos(al_kashi_angle(LO, KO, 500) + np.arccos(MO/LO))
    0.0
    >>> cos_angle_to_v3(v3_to_cos_angle(500)) - 500 < 0.0000001
    True
    '''
    return (KO ** 2 + LO ** 2 - v3 ** 2) / (2 * KO * LO) * MO / LO - np.sqrt(
        1 - ((KO ** 2 + LO ** 2 - v3 ** 2) / (2 * KO * LO)) ** 2) * np.sqrt(1 - (MO / LO) ** 2)


def distance(x, y, z=0):
    '''
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
    '''
    return np.sqrt(x ** 2 + y ** 2 + z ** 2)


def al_kashi_longueur(a, b, alpha):
    return np.sqrt(a ** 2 + b ** 2 - 2 * a * b * np.cos(alpha))


def al_kashi_angle(a, b, c):
    return np.arccos((a ** 2 + b ** 2 - c ** 2) / (2 * a * b))


############################################################################

if __name__ == "__main__":
    import doctest

    doctest.testmod()
