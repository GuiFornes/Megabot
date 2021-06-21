# L'ensemble des distances sont exprimées en mm : segments de patte et élongations des vérins

import numpy as np
from numpy.linalg import inv
from qpsolvers import solve_qp
import kinetic as kin
from upg_tools import *

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