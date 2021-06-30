# L'ensemble des distances sont exprimées en mm : segments de patte et élongations des vérins

from upg_kinetic import *
import numpy as np
from numpy.linalg import inv
from qpsolvers import solve_qp
from upg_tools import *

def distance_3_points(A, B, C):
    """
    Fonction auxiliaire de gen_jacob_2
    Retourne la matrice des distances de A à B et à C
    """
    return 2 * np.array([[A[0] - B[0], A[1] - B[1]], [A[0] - C[0], A[1] - C[1]]])

def gen_jacob_2(pts, lpl, v1, v2):
    """
    Retourne la jacobienne correspondant au modèle cinématique indirect dans le plan de la patte
    Prend en argument la position des points de la patte et l'élongation des verrins en m
    La jacobienne doit être appliqué sur des élongations en m et retourne des position en m

    >>> gen_jacob_2(get_leg_points_V1_V2(0.495, 0.585, LEGS[FL]['lengths']), 0.495, 0.585) @ np.array([0, 0])
    array([0., 0.])

    # >>> gen_jacob_2(get_leg_points_V1_V2(0.495, 0.585, LEGS[FL]['lengths']), 0.495, 0.585) @ np.array([1, 0])
    # array([1, 0])
    #
    # >>> gen_jacob_2(get_leg_points_V1_V2(0.495, 0.585, LEGS[FL]['lengths']), 0.495, 0.585) @ np.array([1, 0])

    """
    x_E, y_E = pts['E']
    x_F, y_F = pts['F']
    x_G, y_G = pts['G']
    x_H, y_H = pts['H']
    x_I, y_I = pts['I']

    A = distance_3_points(pts['D'], pts['A'], pts['C'])
    B = np.array([0, 2 * v1])
    M_D = inv(A) @ B

    M_E = (lpl['ae'] / (lpl['ae'] - lpl['de'])) * M_D

    A = distance_3_points(pts['F'], pts['E'], pts['B'])
    B = np.array([
        [2 * (x_F - x_E), 2 * (y_F - y_E)],
        [0, 0]])
    M_F = (inv(A) @ B) @ M_E

    M_G = ((lpl['ef'] + lpl['fg']) / lpl['ef']) * M_F

    M_H = ((lpl['bf'] - lpl['fh']) / lpl['bf']) * M_F

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

    Jacob = inv((lpl['gj'] / lpl['gi']) * M_I)

    return Jacob

def mat_A(pts, lpl, v1, v2, v3, alpha):
    """
    Fonction auxiliaire de gen_jacob_3
    Génère la matrice A conformément à nos équations (cf. indirect.pdf)
    Toutes les longueurs en m
    """
    Jacob = gen_jacob_2(pts, lpl, v1, v2)

    KO = lpl['yaw_c']
    LM = lpl['yaw_b']
    MO = lpl['yaw_a']
    LO = np.sqrt(LM ** 2 + MO ** 2)

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
    pts = get_leg_points_V1_V2(v1, v2, lpl)
    A = mat_A(pts, lpl, v1, v2, v3, alpha)
    B = mat_B(pts, alpha)

    return A @ inv(B)

def gen_jacob_12(V):
    J_12 = []
    for i in range(4):
        # Calcul de la jacobienne de la patte
        v1, v2, v3 = V[i*3], V[i*3 + 1], V[i*3 + 2]
        lpl = LEGS[i]['lengths']
        alpha = np.arccos(v3_to_cos_angle(v3, lpl))
        J_3 = gen_jacob_3(v1 / 1000, v2 / 1000, v3 / 1000, alpha, lpl) @ LEGS[i]['matrix'].T

        # Ajout à J_12
        L = np.zeros((3, 12))
        for j in range(3):
            for k in range(3):
                L[j][i*3 + k] = J_3[j][k]
        J_12.append(L[0])
        J_12.append(L[1])
        J_12.append(L[2])
    return J_12

def gen_matrix_rota(angle, axis):
    if axis == 0: # Rotation selon x
        return np.array([[1, 0, 0],
                         [0, np.cos(angle), -np.sin(angle)],
                         [0, np.sin(angle), np.cos(angle)]])
    elif axis == 1: # Rotation selon y
        return np.array([[np.cos(angle), 0, -np.sin(angle)],
                         [0, 1, 0],
                         [np.sin(angle), 0, np.cos(angle)]])
    else: # Rotation selon z
        return np.array([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle), np.cos(angle), 0],
                         [0, 0, 1]])

def gen_deriv_matrix_rota(angle, axis):
    if axis == 0: # Rotation selon x
        return np.array([[0, 0, 0],
                         [0, -np.sin(angle), -np.cos(angle)],
                         [0, np.cos(angle), -np.sin(angle)]])
    elif axis == 1: # Rotation selon y
        return np.array([[-np.sin(angle), 0, -np.cos(angle)],
                         [0, 0, 0],
                         [np.cos(angle), 0, -np.sin(angle)]])
    else: # Rotation selon z
        return np.array([[-np.sin(angle), -np.cos(angle), 0],
                         [np.cos(angle), -np.sin(angle), 0],
                         [0, 0, 0]])

def gen_R(l, m, n):
    return gen_matrix_rota(l, 2) @ gen_matrix_rota(m, 0) @ gen_matrix_rota(n, 2)

def gen_dRdl(l, m, n):
    return gen_deriv_matrix_rota(l, 2) @ gen_matrix_rota(m, 0) @ gen_matrix_rota(n, 2)

def gen_dRdm(l, m, n):
    return gen_matrix_rota(l, 2) @ gen_deriv_matrix_rota(m, 0) @ gen_matrix_rota(n, 2)

def gen_dRdn(l, m, n):
    return gen_matrix_rota(l, 2) @ gen_matrix_rota(m, 0) @ gen_deriv_matrix_rota(n, 2)

def gen_jacob_12x18(V, R, dRdl, dRdm, dRdn, X):
    J_l = np.block([[dRdl, np.zeros((3, 9))],
                    [np.zeros((3, 3)), dRdl, np.zeros((3, 6))],
                    [np.zeros((3, 6)), dRdl, np.zeros((3, 3))],
                    [np.zeros((3, 9)), dRdl]]) @ X
    J_m = np.block([[dRdm, np.zeros((3, 9))],
                    [np.zeros((3, 3)), dRdm, np.zeros((3, 6))],
                    [np.zeros((3, 6)), dRdm, np.zeros((3, 3))],
                    [np.zeros((3, 9)), dRdm]]) @ X
    J_n = np.block([[dRdn, np.zeros((3, 9))],
                    [np.zeros((3, 3)), dRdn, np.zeros((3, 6))],
                    [np.zeros((3, 6)), dRdn, np.zeros((3, 3))],
                    [np.zeros((3, 9)), dRdn]]) @ X
    J_Omega = np.concatenate((-J_l.reshape((12,1)), -J_m.reshape((12,1)), -J_n.reshape((12,1))), axis=1)

    J_O = np.concatenate((-np.eye(3), -np.eye(3), -np.eye(3), -np.eye(3)))

    M = np.concatenate((np.eye(12), J_Omega, J_O), axis=1)

    J_12 = gen_jacob_12(V)

    Big_inv_R = np.block([[R.T, np.zeros((3, 9))],
                          [np.zeros((3, 3)), R.T, np.zeros((3, 6))],
                          [np.zeros((3, 6)), R.T, np.zeros((3, 3))],
                          [np.zeros((3, 9)), R.T]])

    return J_12 @ Big_inv_R @ M

def pseudo_inv(M):
    return inv(M.T @ M) @ M.T

######################### 1 ########################

def legs_constraints():
    C = np.zeros((12, 18))
    for i in range(4):
        if get_og(i):
            for j in range(3):
                C[i * 3 + j][i * 3 + j] = 1
    return C

def gen_jacob_24x18(V, X, l, m, n):
    """
    A voir
    """
    R = gen_R(l, m, n)
    dRdl = gen_dRdl(l, m, n)
    dRdm = gen_dRdm(l, m, n)
    dRdn = gen_dRdn(l, m, n)
    J_12x18 = gen_jacob_12x18(V, R, dRdl, dRdm, dRdn, X)
    Legs_constraints = legs_constraints()
    return np.concatenate((J_12x18, Legs_constraints))

def gen_jacob(V, X, l, m, n):
    J_24x18 = gen_jacob_24x18(V, X, l, m, n)
    J_24x15 = J_24x18[0:24, 0:15]
    J_O = J_24x18[0:24, 15:18]
    return pseudo_inv(J_O) @ np.concatenate((-J_24x15, np.eye(24)), axis=1)

######################### 2 ########################

def legs_constraints_2():
    """
    3 pattes au sol
    """
    C = np.zeros((3, 18))
    l = 0
    for i in range(4):
        if l == 3: break
        if get_og(i):
            C[l][i * 3 + 2] = 1
            l += 1
    return C

def gen_jacob_15x18(V, X, l, m, n):
    """
    A voir
    """
    R = gen_R(l, m, n)
    dRdl = gen_dRdl(l, m, n)
    dRdm = gen_dRdm(l, m, n)
    dRdn = gen_dRdn(l, m, n)
    J_12x18 = gen_jacob_12x18(V, R, dRdl, dRdm, dRdn, X)
    Legs_constraints = legs_constraints_2()
    return np.concatenate((J_12x18, Legs_constraints))

def gen_jacob_alt(V, X, l, m, n):
    J_15x18 = gen_jacob_15x18(V, X, l, m, n)
    J_15x15 = J_15x18[0:15, 0:15]
    J_O = J_15x18[0:15, 15:18]
    return pseudo_inv(J_O) @ np.concatenate((-J_15x15, np.eye(15)), axis=1)

#####################################################################################

def gen_VAR(V, Omega, X_rel):
    dRdl = gen_dRdl(Omega[0], Omega[1], Omega[2])
    dRdm = gen_dRdm(Omega[0], Omega[1], Omega[2])
    dRdn = gen_dRdn(Omega[0], Omega[1], Omega[2])
    J_l = np.block([[dRdl, np.zeros((3, 9))],
                    [np.zeros((3, 3)), dRdl, np.zeros((3, 6))],
                    [np.zeros((3, 6)), dRdl, np.zeros((3, 3))],
                    [np.zeros((3, 9)), dRdl]]) @ X_rel
    J_m = np.block([[dRdm, np.zeros((3, 9))],
                    [np.zeros((3, 3)), dRdm, np.zeros((3, 6))],
                    [np.zeros((3, 6)), dRdm, np.zeros((3, 3))],
                    [np.zeros((3, 9)), dRdm]]) @ X_rel
    J_n = np.block([[dRdn, np.zeros((3, 9))],
                    [np.zeros((3, 3)), dRdn, np.zeros((3, 6))],
                    [np.zeros((3, 6)), dRdn, np.zeros((3, 3))],
                    [np.zeros((3, 9)), dRdn]]) @ X_rel

    J_12 = gen_jacob_12(V)
    R = gen_R(Omega[0], Omega[1], Omega[2])
    Big_R = np.block([[R, np.zeros((3, 9))],
                          [np.zeros((3, 3)), R, np.zeros((3, 6))],
                          [np.zeros((3, 6)), R, np.zeros((3, 3))],
                          [np.zeros((3, 9)), R]])
    return np.concatenate((np.concatenate((Big_R @ inv(J_12), J_l.reshape((12, 1)), J_m.reshape((12, 1)), J_n.reshape((12, 1))), axis=1), np.zeros((3,15))))

def legs_constraints_3():
    """
    3 pattes au sol
    """
    C = np.zeros((3, 15))
    l = 0
    for i in range(4):
        if l == 3: break
        if get_og(i):
            C[l][i * 3 + 2] = 1
            l += 1
    return C

def legs_constraints_4():
    """
    1 patte au sol
    """
    C = np.zeros((3, 15))
    for i in range(4):
        if get_og(i):
            C[0][i * 3] = 1
            C[1][i * 3 + 1] = 1
            C[2][i * 3 + 2] = 1
            break
    return C

def gen_OBJ():
    return np.concatenate((np.concatenate((np.eye(12), - np.concatenate((np.eye(3), np.eye(3), np.eye(3), np.eye(3)))), axis=1), legs_constraints_4()))

def gen_new_jacob(V, Omega, X_rel):
    return inv(gen_OBJ()) @ gen_VAR(V, Omega, X_rel)


# V = [550, 550, 550, 550, 550, 550, 550, 550, 550, 550, 550, 550]
# print(gen_new_jacob(V, np.array([0, 0, 0])))