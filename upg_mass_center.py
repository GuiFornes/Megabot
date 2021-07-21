import numpy as np

import kinetic as kin
from upg_tools import *
from upg_jacobian import *

body_weight = 140.0

leg_weights = {}
for p in 'ABCDEFGHIJ':
    leg_weights[p] = {}
leg_weights['A']['B'] = 2.0
leg_weights['A']['C'] = 1.6
leg_weights['A']['E'] = 3.6
leg_weights['E']['G'] = 4.0
leg_weights['B']['F'] = 3.45
leg_weights['G']['J'] = 3.65
leg_weights['D']['C'] = 5.3
leg_weights['H']['I'] = 5.3

leg_tot_weight = 0
for a in leg_weights:
    for b in leg_weights[a]:
        leg_tot_weight += leg_weights[a][b]

for a in leg_weights:
    for b in leg_weights[a]:
        leg_weights[b][a] = leg_weights[a][b]


################## CALCUL DU PROJETE DU CENTRE DE MASSE ##################

def leg_center_of_mass(v1, v2, leg_id):
    """
    Retourne le centre de masse de la patte leg_id dans son référentiel et sa masse équivalente

    :param v1: élongation du vérin 1 de la patte
    :param v2: élongation du vérin 2 de la patte
    :param leg_id: ID de la patte
    """
    points = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, ROBOT['legs'][leg_id]['lengths'])
    lweights = []
    total_weight = 0
    application_points = []
    G = 0
    for i in range(len(points)):
        pt1 = list(points)[i]
        for j in range(i + 1, len(points)):
            pt2 = list(points)[j]
            if pt1 in leg_weights and pt2 in leg_weights[pt1]:
                lweights.append((leg_weights[pt1][pt2], pt1, pt2))
    for w, pt1, pt2 in lweights:
        m = np.mean([points[pt1], points[pt2]], axis=0)
        application_points.append((m, w))
    for m, w in application_points:
        G += m * w
        total_weight += w
    G /= total_weight
    return G * 1000, total_weight


def body_center_of_mass(passenger=True, passenger_weight=80.0):
    """
    Retourne le centre de masse du châssis dans le référentiel du robot et sa masse équivalente
    Peut prendre un passager en argument, et le cas échéant fixer sa masse

    :param passenger: présence d'un passager
    :param passenger_weight: si passager il y a, sa masse
    :return: centre de masse et masse équivalente du châssis du robot
    """
    if passenger:
        total_weight = body_weight + passenger_weight
        return np.array([0.0, 0.0, 750.0 * passenger_weight]) / total_weight, total_weight
    return np.array([0.0, 0.0, 0.0]), body_weight


def center_of_mass(V, passenger=True, passenger_weight=80.0):
    """
    Retourne les coordonnées du centre de masse du robot dans le référentiel du robot

    :param V: élongations des vérins (vecteur de taille 12)
    :param passenger: présence d'un passager
    :param passenger_weight: si passager il y a, sa masse
    :return: centre de masse du robot dans son référentiel
    """
    m, w = body_center_of_mass(passenger=passenger, passenger_weight=passenger_weight)
    G = m * w
    total_weight = w

    for i in range(4):
        m, w = leg_center_of_mass(V[3 * i], V[3 * i + 1], i)
        G = G + leg_ref_to_robot(d2_to_d3(m[0], m[1], v3_to_cos_angle(V[3 * i + 2], ROBOT['legs'][i]['lengths'])),
                                 i) * w
        total_weight += w

    return G / total_weight


def plan_from_points(A, B, C):
    """
    Détermine une équation du plan contenant les 3 points non colinéaires A, B et C

    :param A: 1er point
    :param B: 2nd point
    :param C: 3ème point
    :return: (a, b, c, d) tel que le plan contenant A, B et C ait pour équation ax + by + cz + d = 0

    >>> eq = plan_from_points(np.array([1, 0, 0]), np.array([0, 0, 0]), np.array([0, 1, 0]))
    >>> eq[0] * 1 + eq[1] * 1 + eq[2] * 0 + eq[3]
    0.0
    >>> eq = plan_from_points(np.array([1, 0, 1]), np.array([0, 0, 0]), np.array([0, 1, 1]))
    >>> eq[0] * 1 + eq[1] * 1 + eq[2] * 2 + eq[3]
    0.0
    >>> eq = plan_from_points(np.array([0, 1, 1]), np.array([0, 0, 1]), np.array([1, 0, 0]))
    >>> eq[0] * 1 + eq[1] * 1 + eq[2] * 0 + eq[3]
    0.0
    """
    AB = B - A
    AC = C - A
    if AB[1] * AC[2] - AC[1] * AB[2] > 0.001 and AB[1] > 0.001 and AC[2] > 0.001:  # a = 1
        b = (AC[0] * AB[2] - AB[0] * AC[2]) / (AB[1] * AC[2] - AC[1] * AB[2])
        c = - (b * AC[1] + AC[0]) / AC[2]
        d = - (A[0] + b * A[1] + c * A[2])
        return 1., b, c, d
    elif AB[0] * AC[2] - AC[0] * AB[2] > 0.001 and AB[0] > 0.001 and AC[2] > 0.001:  # b = 1
        a = (AC[1] * AB[2] - AB[1] * AC[2]) / (AB[0] * AC[2] - AC[0] * AB[2])
        c = - (a * AC[0] + AC[1]) / AC[2]
        d = - (a * A[0] + A[1] + c * A[2])
        return a, 1., c, d
    b = (AC[2] * AB[0] - AB[2] * AC[0]) / (AB[1] * AC[0] - AC[1] * AB[0])  # c = 1
    a = - (b * AC[1] + AC[2]) / AC[0]
    d = - (a * A[0] + b * A[1] + A[2])
    return a, b, 1., d


def project(point):
    """
    Calcule le projeté du point sur le sol dans le référentiel du robot

    :param point: coordonnées du point à projeter sur le sol dans le référentiel du robot
    :return: projeté du point sur le sol dans le référnetiel du robot
    """
    V = get_verins_12()
    legs_on_ground = []
    for i in range(4):
        if ROBOT['legs'][i]['og']:
            lpl = ROBOT['legs'][i]['lengths']
            Pos2D = get_leg_points_V1_V2(V[3 * i] / 1000, V[3 * i + 1] / 1000, lpl)['J']
            legs_on_ground.append(
                leg_ref_to_robot(d2_to_d3(Pos2D[0] * 1000, Pos2D[1] * 1000, v3_to_cos_angle(V[3 * i + 2], lpl)), i))
        if len(legs_on_ground) == 3:
            break
    plan = plan_from_points(legs_on_ground[0], legs_on_ground[1], legs_on_ground[2])

    k = - (plan[0] * point[0] + plan[1] * point[1] + plan[2] * point[2] + plan[3]) / (
            plan[0] ** 2 + plan[1] ** 2 + plan[2] ** 2)

    return plan[0] * k + point[0], plan[1] * k + point[1], plan[2] * k + point[2]


def proj_center_of_mass_robot(passenger=True, passenger_weight=80.0):
    """
    Calcule le projeté au sol du centre de masse du robot dans son référentiel

    :param passenger: présence d'un passager
    :param passenger_weight: si passager il y a, sa masse
    :return: projeté au sol du centre de masse dans le référentiel du robot
    """
    return project(center_of_mass(get_verins_12(), passenger=passenger, passenger_weight=passenger_weight))


############### CALCUL DU CENTRE DE MASSE PAR JACOBIENNES ################

def gen_list_M_leg(v1, v2, v3, leg_id, pts):
    """
    Génère la liste des matrices de passage permettant d'obtenir les deltas de position des points de la patte
    leg_id (dx, dy, dz) dans le référentiel du robot à partir des deltas d'élongation des vérins (dv1, dv2, dv3)
    de la patte. Chaque matrice est accessible via son nom ('M_A' pour A, 'M_B' pour B, etc).
    Les élongations des vérins sont en m

    :param v1: élongation de v1 en m
    :param v2: élongation de v2 en m
    :param v3: élongation de v3 en m
    :return: ['M_A': matrice 3x3, 'M_B': matrice 3x3,  ...]
    """
    lpl = ROBOT['legs'][leg_id]['lengths']
    alpha = np.arccos(v3_to_cos_angle(v3*1000, lpl))

    M_A = np.zeros((3, 3))
    M_B = mat_B(pts['B'][0], alpha) @ mat_A(np.zeros((2, 2)), lpl, v3, alpha)
    M_C = mat_B(pts['C'][0], alpha) @ mat_A(np.zeros((2, 2)), lpl, v3, alpha)
    M_D = mat_B(pts['D'][0], alpha) @ \
          mat_A(np.concatenate((np.reshape(gen_MD(pts, v1), (2,1)), np.zeros((2, 1))), axis=1), lpl, v3, alpha)
    M_E = mat_B(pts['E'][0], alpha) @ \
          mat_A(np.concatenate((np.reshape(gen_ME(pts, lpl, v1), (2,1)), np.zeros((2, 1))), axis=1), lpl, v3, alpha)
    M_F = mat_B(pts['F'][0], alpha) @ \
          mat_A(np.concatenate((np.reshape(gen_MF(pts, lpl, v1), (2,1)), np.zeros((2, 1))), axis=1), lpl, v3, alpha)
    M_G = mat_B(pts['G'][0], alpha) @ \
          mat_A(np.concatenate((np.reshape(gen_MG(pts, lpl, v1), (2,1)), np.zeros((2, 1))), axis=1), lpl, v3, alpha)
    M_H = mat_B(pts['H'][0], alpha) @ \
          mat_A(np.concatenate((np.reshape(gen_MH(pts, lpl, v1), (2,1)), np.zeros((2, 1))), axis=1), lpl, v3, alpha)
    M_I = mat_B(pts['I'][0], alpha) @ mat_A(gen_MI(pts, lpl, v1, v2), lpl, v3, alpha)
    M_J = mat_B(pts['J'][0], alpha) @ mat_A(gen_MJ(pts, lpl, v1, v2), lpl, v3, alpha)

    return {'M_A': ROBOT['legs'][leg_id]['matrix'].T @ M_A,
            'M_B': ROBOT['legs'][leg_id]['matrix'].T @ M_B,
            'M_C': ROBOT['legs'][leg_id]['matrix'].T @ M_C,
            'M_D': ROBOT['legs'][leg_id]['matrix'].T @ M_D,
            'M_E': ROBOT['legs'][leg_id]['matrix'].T @ M_E,
            'M_F': ROBOT['legs'][leg_id]['matrix'].T @ M_F,
            'M_G': ROBOT['legs'][leg_id]['matrix'].T @ M_G,
            'M_H': ROBOT['legs'][leg_id]['matrix'].T @ M_H,
            'M_I': ROBOT['legs'][leg_id]['matrix'].T @ M_I,
            'M_J': ROBOT['legs'][leg_id]['matrix'].T @ M_J}


def gen_J_com_leg(v1, v2, v3, leg_id):
    pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, ROBOT['legs'][leg_id]['lengths'])
    lweights = []
    M = gen_list_M_leg(v1/1000, v2/1000, v3/1000, leg_id, pts)
    J = np.zeros((3, 3))
    for i in range(len(pts)):
        pt1 = list(pts)[i]
        for j in range(i + 1, len(pts)):
            pt2 = list(pts)[j]
            if pt1 in leg_weights and pt2 in leg_weights[pt1]:
                lweights.append((leg_weights[pt1][pt2], pt1, pt2))
    for w, pt1, pt2 in lweights:
        J = J + (w / (2 * leg_tot_weight)) * M['M_' + pt1] + (w / (2 * leg_tot_weight)) * M['M_' + pt2]
    return J


def gen_J_com_rel(V, passenger_weight=80):
    tot_weight = 4 * leg_tot_weight + body_weight + passenger_weight
    return np.concatenate((leg_tot_weight / tot_weight * gen_J_com_leg(V[0], V[1], V[2], 0),
                           leg_tot_weight / tot_weight * gen_J_com_leg(V[3], V[4], V[5], 1),
                           leg_tot_weight / tot_weight * gen_J_com_leg(V[6], V[7], V[8], 2),
                           leg_tot_weight / tot_weight * gen_J_com_leg(V[9], V[10], V[11], 3)), axis=1)


def gen_J_com_abs(V, Omega, com, passenger_weight=80):
    J_rel = gen_J_com_rel(V, passenger_weight=passenger_weight)
    return np.concatenate((gen_R(Omega[0], Omega[1], Omega[2]) @ J_rel,
                           np.eye(3),
                           np.reshape(gen_dRdl(Omega[0], Omega[1], Omega[2]) @ com, (3, 1)),
                           np.reshape(gen_dRdm(Omega[0], Omega[1], Omega[2]) @ com, (3, 1)),
                           np.reshape(gen_dRdn(Omega[0], Omega[1], Omega[2]) @ com, (3, 1))), axis=1)


############################################################################

if __name__ == "__main__":
    import doctest

    doctest.testmod()
