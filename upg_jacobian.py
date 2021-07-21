# L'ensemble des distances sont exprimées en mm : segments de patte et élongations des vérins

from upg_tools import *


################# CINEMATIQUE 2D DANS LE REF DE LA PATTE #################

def distance_3_points(A, B, C):
    """
    Fonction auxiliaire de gen_MJ
    Retourne la matrice des distances de A à B et à C
    """
    return 2 * np.array([[A[0] - B[0], A[1] - B[1]], [A[0] - C[0], A[1] - C[1]]])


def gen_MD(pts, v1):
    A = distance_3_points(pts['D'], pts['A'], pts['C'])
    B = np.array([0, 2 * v1])
    return inv(A) @ B


def gen_ME(pts, lpl, v1):
    return (lpl['ae'] / (lpl['ae'] - lpl['de'])) * gen_MD(pts, v1)


def gen_MF(pts, lpl, v1):
    x_E, y_E = pts['E']
    x_F, y_F = pts['F']

    A = distance_3_points(pts['F'], pts['E'], pts['B'])
    B = np.array([[2 * (x_F - x_E), 2 * (y_F - y_E)],
                  [0, 0]])
    return (inv(A) @ B) @ gen_ME(pts, lpl, v1)


def gen_MG(pts, lpl, v1):
    return ((lpl['ef'] + lpl['fg']) / lpl['ef']) * gen_MF(pts, lpl, v1)


def gen_MH(pts, lpl, v1):
    return ((lpl['bf'] - lpl['fh']) / lpl['bf']) * gen_MF(pts, lpl, v1)


def gen_MI(pts, lpl, v1, v2):
    x_G, y_G = pts['G']
    x_H, y_H = pts['H']
    x_I, y_I = pts['I']

    A = distance_3_points(pts['I'], pts['G'], pts['H'])
    B = np.array([[2 * (x_I - x_G), 2 * (y_I - y_G)],
                  [0, 0]])
    C = np.array([[0, 0],
                  [2 * (x_I - x_H), 2 * (y_I - y_H)]])
    D = np.array([0, 2 * v2])

    V1 = inv(A) @ (B @ gen_MG(pts, lpl, v1) + C @ gen_MH(pts, lpl, v1))
    V2 = inv(A) @ D
    return np.array([[V1[0], V2[0]],
                     [V1[1], V2[1]]])


def gen_MJ(pts, lpl, v1, v2):
    """
    Retourne la jacobienne correspondant au modèle cinématique inverse dans le plan de la patte
    Prend en argument la position des points de la patte et l'élongation des verrins en m
    La jacobienne doit être appliqué sur des élongations en m et retourne des position en m

    >>> gen_MJ(get_leg_points_V1_V2(0.495, 0.585, ROBOT['legs'][FL]['lengths']), ROBOT['legs'][FL]['lengths'], 0.495, 0.585) @ np.array([0, 0])
    array([0., 0.])

    # >>> gen_MJ(get_leg_points_V1_V2(0.495, 0.585, ROBOT['legs'][FL]['lengths']), ROBOT['legs'][FL]['lengths'], 0.495, 0.585) @ np.array([1, 0])
    # array([1, 0])
    #
    # >>> gen_MJ(get_leg_points_V1_V2(0.495, 0.585, ROBOT['legs'][FL]['lengths']), ROBOT['legs'][FL]['lengths'], 0.495, 0.585) @ np.array([1, 0])

    """
    return (lpl['gj'] / lpl['gi']) * gen_MI(pts, lpl, v1, v2)


################# CINEMATIQUE 3D DANS LE REF DE LA PATTE #################

def mat_A(J, lpl, v3, alpha):
    """
    Fonction auxiliaire de gen_jacob_3
    Génère la matrice A conformément à nos équations (cf. indirect.pdf)
    Toutes les longueurs en m
    """

    KO = lpl['yaw_c']
    LM = lpl['yaw_b']
    MO = lpl['yaw_a']
    LO = np.sqrt(LM ** 2 + MO ** 2)

    A = np.array([
        [J[0][0], J[0][1], 0],
        [J[1][0], J[1][1], 0],
        [0, 0, v3 / (np.sin(alpha - np.arccos(MO / LO)) * (KO / 1000) * (LO / 1000))]])

    return A


def mat_B(X, alpha):
    """
    Fonction auxiliaire de gen_jacob_3
    Génère la matrice B conformément à nos équations (cf. indirect.pdf)
    Toutes les longueurs en m
    """
    B = np.array([
        [np.cos(np.pi / 2 - alpha), 0, X * np.sin(np.pi / 2 - alpha)],
        [np.sin(np.pi / 2 - alpha), 0, -X * np.cos(np.pi / 2 - alpha)],
        [0, 1, 0]])

    return B


def gen_jacob_3(v1, v2, v3, alpha, lpl):
    """
    Retourne la jacobienne correspondant au modèle cinématique inverse dans le repère cartésien centré en O
    Prend en argument l'élongation des verrins en m et l'angle alpha en radian
    La jacobienne doit être appliquée sur des élongations en m et retourne des position en m
    """
    pts = get_leg_points_V1_V2(v1, v2, lpl)
    J = gen_MJ(pts, lpl, v1, v2)
    A = mat_A(J, lpl, v3, alpha)
    B = mat_B(pts['J'][0], alpha)

    return B @ A


#################### CINEMATIQUE DANS LE REF DU ROBOT ####################

def gen_jacob_12(V):
    J = []
    for i in range(4):
        # Calcul de la jacobienne de la patte
        v1, v2, v3 = V[i * 3], V[i * 3 + 1], V[i * 3 + 2]
        lpl = ROBOT['legs'][i]['lengths']
        alpha = np.arccos(v3_to_cos_angle(v3, lpl))
        J.append(ROBOT['legs'][i]['matrix'].T @ gen_jacob_3(v1 / 1000, v2 / 1000, v3 / 1000, alpha, lpl))

    return np.concatenate((np.concatenate((J[0], np.zeros((3, 9))), axis=1),
                           np.concatenate((np.zeros((3, 3)), J[1], np.zeros((3, 6))), axis=1),
                           np.concatenate((np.zeros((3, 6)), J[2], np.zeros((3, 3))), axis=1),
                           np.concatenate((np.zeros((3, 9)), J[3]), axis=1)))


def jacob_dX_to_dV_dO_dOmega(V, Omega, X_rel):
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
    return np.concatenate((Big_R @ J_12,
                           np.concatenate((np.eye(3), np.eye(3), np.eye(3), np.eye(3))),
                           np.reshape(J_l, (12, 1)),
                           np.reshape(J_m, (12, 1)),
                           np.reshape(J_n, (12, 1))), axis=1)


def jacob_dPf_dO_dOmega_to_dV_dPg_dO_dOmega(leg_id, V, Omega, X_rel):
    v = []
    Pf = np.zeros((3, 1))
    for j in range(3):
        v.append(V[3 * leg_id + j])
        Pf[j][0] = X_rel[3 * leg_id + j]
    lpl = ROBOT['legs'][leg_id]['lengths']
    Jf = ROBOT['legs'][leg_id]['matrix'].T @ \
         gen_jacob_3(v[0] / 1000, v[1] / 1000, v[2] / 1000, np.arccos(v3_to_cos_angle(v[2], lpl)), lpl)
    R = gen_R(Omega[0], Omega[1], Omega[2])
    dRdl = gen_dRdl(Omega[0], Omega[1], Omega[2])
    dRdm = gen_dRdm(Omega[0], Omega[1], Omega[2])
    dRdn = gen_dRdn(Omega[0], Omega[1], Omega[2])
    if leg_id == 0:
        JPf = np.concatenate((R @ Jf, np.zeros((3, 18)),
                              np.eye(3), dRdl @ Pf, dRdm @ Pf, dRdn @ Pf), axis=1)
    elif leg_id == 3:
        JPf = np.concatenate((np.zeros((3, 9)), R @ Jf, np.zeros((3, 9)),
                              np.eye(3), dRdl @ Pf, dRdm @ Pf, dRdn @ Pf), axis=1)
    else:
        JPf = np.concatenate((np.zeros((3, 3 * leg_id)), R @ Jf, np.zeros((3, 3 * (6 - leg_id))),
                              np.eye(3), dRdl @ Pf, dRdm @ Pf, dRdn @ Pf), axis=1)
    JO = np.concatenate((np.zeros((3, 21)), np.eye(3), np.zeros((3, 3))), axis=1)
    JOmega = np.concatenate((np.zeros((3, 24)), np.eye(3)), axis=1)
    return np.concatenate((JPf, JO, JOmega))

# def gen_jacob_12x18(V, R, dRdl, dRdm, dRdn, X):
#     J_l = np.block([[dRdl, np.zeros((3, 9))],
#                     [np.zeros((3, 3)), dRdl, np.zeros((3, 6))],
#                     [np.zeros((3, 6)), dRdl, np.zeros((3, 3))],
#                     [np.zeros((3, 9)), dRdl]]) @ X
#     J_m = np.block([[dRdm, np.zeros((3, 9))],
#                     [np.zeros((3, 3)), dRdm, np.zeros((3, 6))],
#                     [np.zeros((3, 6)), dRdm, np.zeros((3, 3))],
#                     [np.zeros((3, 9)), dRdm]]) @ X
#     J_n = np.block([[dRdn, np.zeros((3, 9))],
#                     [np.zeros((3, 3)), dRdn, np.zeros((3, 6))],
#                     [np.zeros((3, 6)), dRdn, np.zeros((3, 3))],
#                     [np.zeros((3, 9)), dRdn]]) @ X
#     J_Omega = np.concatenate((-J_l.reshape((12, 1)), -J_m.reshape((12, 1)), -J_n.reshape((12, 1))), axis=1)
#
#     J_O = np.concatenate((-np.eye(3), -np.eye(3), -np.eye(3), -np.eye(3)))
#
#     M = np.concatenate((np.eye(12), J_Omega, J_O), axis=1)
#
#     J_12 = gen_jacob_12(V)
#
#     Big_inv_R = np.block([[R.T, np.zeros((3, 9))],
#                           [np.zeros((3, 3)), R.T, np.zeros((3, 6))],
#                           [np.zeros((3, 6)), R.T, np.zeros((3, 3))],
#                           [np.zeros((3, 9)), R.T]])
#
#     return J_12 @ Big_inv_R @ M
#
#
# ######################### 1 ########################
#
# def legs_constraints():
#     C = np.zeros((12, 18))
#     for i in range(4):
#         if get_og(i):
#             for j in range(3):
#                 C[i * 3 + j][i * 3 + j] = 1
#     return C
#
# def gen_jacob_24x18(V, X, l, m, n):
#     """
#     A voir
#     """
#     R = gen_R(l, m, n)
#     dRdl = gen_dRdl(l, m, n)
#     dRdm = gen_dRdm(l, m, n)
#     dRdn = gen_dRdn(l, m, n)
#     J_12x18 = gen_jacob_12x18(V, R, dRdl, dRdm, dRdn, X)
#     Legs_constraints = legs_constraints()
#     return np.concatenate((J_12x18, Legs_constraints))
#
# def gen_jacob(V, X, l, m, n):
#     J_24x18 = gen_jacob_24x18(V, X, l, m, n)
#     J_24x15 = J_24x18[0:24, 0:15]
#     J_O = J_24x18[0:24, 15:18]
#     return pseudo_inv(J_O) @ np.concatenate((-J_24x15, np.eye(24)), axis=1)
#
# ######################### 2 ########################
#
# def legs_constraints_2():
#     """
#     3 pattes au sol
#     """
#     C = np.zeros((3, 18))
#     l = 0
#     for i in range(4):
#         if l == 3: break
#         if get_og(i):
#             C[l][i * 3 + 2] = 1
#             l += 1
#     return C
#
# def gen_jacob_15x18(V, X, l, m, n):
#     """
#     A voir
#     """
#     R = gen_R(l, m, n)
#     dRdl = gen_dRdl(l, m, n)
#     dRdm = gen_dRdm(l, m, n)
#     dRdn = gen_dRdn(l, m, n)
#     J_12x18 = gen_jacob_12x18(V, R, dRdl, dRdm, dRdn, X)
#     Legs_constraints = legs_constraints_2()
#     return np.concatenate((J_12x18, Legs_constraints))
#
# def gen_jacob_alt(V, X, l, m, n):
#     J_15x18 = gen_jacob_15x18(V, X, l, m, n)
#     J_15x15 = J_15x18[0:15, 0:15]
#     J_O = J_15x18[0:15, 15:18]
#     return pseudo_inv(J_O) @ np.concatenate((-J_15x15, np.eye(15)), axis=1)
#
# #####################################################################################
#
# def gen_VAR(V, Omega, X_rel):
#     dRdl = gen_dRdl(Omega[0], Omega[1], Omega[2])
#     dRdm = gen_dRdm(Omega[0], Omega[1], Omega[2])
#     dRdn = gen_dRdn(Omega[0], Omega[1], Omega[2])
#     J_l = np.block([[dRdl, np.zeros((3, 9))],
#                     [np.zeros((3, 3)), dRdl, np.zeros((3, 6))],
#                     [np.zeros((3, 6)), dRdl, np.zeros((3, 3))],
#                     [np.zeros((3, 9)), dRdl]]) @ X_rel
#     J_m = np.block([[dRdm, np.zeros((3, 9))],
#                     [np.zeros((3, 3)), dRdm, np.zeros((3, 6))],
#                     [np.zeros((3, 6)), dRdm, np.zeros((3, 3))],
#                     [np.zeros((3, 9)), dRdm]]) @ X_rel
#     J_n = np.block([[dRdn, np.zeros((3, 9))],
#                     [np.zeros((3, 3)), dRdn, np.zeros((3, 6))],
#                     [np.zeros((3, 6)), dRdn, np.zeros((3, 3))],
#                     [np.zeros((3, 9)), dRdn]]) @ X_rel
#
#     J_12 = gen_jacob_12(V)
#     R = gen_R(Omega[0], Omega[1], Omega[2])
#     Big_R = np.block([[R, np.zeros((3, 9))],
#                           [np.zeros((3, 3)), R, np.zeros((3, 6))],
#                           [np.zeros((3, 6)), R, np.zeros((3, 3))],
#                           [np.zeros((3, 9)), R]])
#     old_VAR = np.concatenate((np.concatenate((Big_R @ inv(J_12), J_l.reshape((12, 1)), J_m.reshape((12, 1)), J_n.reshape((12, 1))), axis=1), np.zeros((3,15))))
#     return np.concatenate((old_VAR, np.concatenate((np.zeros((12, 3)), np.eye(3)))), axis=1)
#
# def legs_constraints_3():
#     """
#     3 pattes au sol
#     """
#     C = np.zeros((3, 15))
#     l = 0
#     for i in range(4):
#         if l == 3: break
#         if get_og(i):
#             C[l][i * 3 + 2] = 1
#             l += 1
#     return C
#
# def legs_constraints_4():
#     """
#     1 patte au sol
#     """
#     C = np.zeros((3, 15))
#     for i in range(4):
#         if get_og(i):
#             C[0][i * 3] = 1
#             C[1][i * 3 + 1] = 1
#             C[2][i * 3 + 2] = 1
#             break
#     return C
#
# def legs_constraints_5():
#     """
#     3 pattes fixes
#     """
#     C = np.zeros((9, 15))
#     c = 0
#     for i in range(4):
#         if get_og(i):
#             c += 1
#             C[3 * i][i * 3] = 1
#             C[3 * i + 1][i * 3 + 1] = 1
#             C[3 * i + 2][i * 3 + 2] = 1
#         if c == 3: break
#     return C
#
# def gen_OBJ():
#     return np.concatenate((np.concatenate((np.eye(12), - np.concatenate((np.eye(3), np.eye(3), np.eye(3), np.eye(3)))), axis=1), legs_constraints_4()))
#
# def gen_weights_regu(V_weights=0.1):
#     W = np.zeros((18, 18))
#     for i in range(12):
#         W[i][i] = V_weights
#     return W
#
# def gen_new_jacob(V, Omega, X_rel):
#     return pseudo_inv(gen_OBJ()) @ gen_VAR(V, Omega, X_rel)
#

#####################################################
