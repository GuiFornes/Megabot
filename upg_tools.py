# L'ensemble des distances sont exprimées en mm : segments de patte et élongations des vérins

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