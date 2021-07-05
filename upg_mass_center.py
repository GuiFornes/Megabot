import kinetic as kin
from upg_tools import *
import numpy as np

body_weight = 140.0

leg_weights = {}
for p in 'ABCDEFGHIJO':
    leg_weights[p]={}
leg_weights['A']['B']=2.0
leg_weights['A']['C']=1.6
leg_weights['A']['E']=3.6
leg_weights['E']['G']=4.0
leg_weights['B']['F']=3.45
leg_weights['G']['J']=3.65
leg_weights['D']['C']=5.3
leg_weights['H']['I']=5.3
tweights = {}
for a in leg_weights:
    for b in leg_weights[a]:
        leg_weights[b][a]=leg_weights[a][b]

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
    G=0
    for i in range(len(points)):
        pt1 = list(points)[i]
        for j in range(i+1, len(points)):
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
        G += leg_ref_to_robot(d2_to_d3(m[0], m[1], v3_to_cos_angle(V[3 * i + 2], ROBOT['legs'][i]['lengths'])), i) * w
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
    if AB[1] * AC[2] - AC[1] * AB[2] > 0.001 and AB[1] > 0.001 and AC[2] > 0.001 : # a = 1
        b = (AC[0] * AB[2] - AB[0] * AC[2]) / (AB[1] * AC[2] - AC[1] * AB[2])
        c = - (b * AC[1] + AC[0]) / AC[2]
        d = - (A[0] + b * A[1] + c * A[2])
        return 1., b, c, d
    elif AB[0] * AC[2] - AC[0] * AB[2] > 0.001 and AB[0] > 0.001 and AC[2] > 0.001 : # b = 1
        a = (AC[1] * AB[2] - AB[1] * AC[2]) / (AB[0] * AC[2] - AC[0] * AB[2])
        c = - (a * AC[0] + AC[1]) / AC[2]
        d = - (a * A[0] + A[1] + c * A[2])
        return a, 1., c, d
    b = (AC[2] * AB[0] - AB[2] * AC[0]) / (AB[1] * AC[0] - AC[1] * AB[0]) # c = 1
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
            legs_on_ground.append(leg_ref_to_robot(d2_to_d3(Pos2D[0] * 1000, Pos2D[1] * 1000, v3_to_cos_angle(V[3 * i + 2], lpl)), i))
        if len(legs_on_ground) == 3:
            break
    plan = plan_from_points(legs_on_ground[0], legs_on_ground[1], legs_on_ground[2])

    k = - (plan[0] * point[0] + plan[1] * point[1] + plan[2] * point[2] + plan[3]) / (plan[0]**2 + plan[1]**2 + plan[2]**2)

    return plan[0] * k + point[0], plan[1] * k + point[1], plan[2] * k + point[2]

def proj_center_of_mass_robot(passenger=True, passenger_weight=80.0):
    """
    Calcule le projeté au sol du centre de masse du robot dans son référentiel

    :param passenger: présence d'un passager
    :param passenger_weight: si passager il y a, sa masse
    :return: projeté au sol du centre de masse dans le référentiel du robot
    """
    return project(center_of_mass(get_verins_12(), passenger=passenger, passenger_weight=passenger_weight))

print(proj_center_of_mass_robot())

############################################################################
if __name__ == "__main__":
    import doctest

    doctest.testmod()