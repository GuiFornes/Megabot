# L'ensemble des distances sont exprimées en mm : segments de patte et élongations des vérins

import numpy as np
from numpy.linalg import inv
from qpsolvers import solve_qp
from kinetic import get_leg_points_V1_V2

FL = 0  # front left leg
FR = 1  # front right leg
RL = 2  # rear left leg
RR = 3  # rear right leg

ROBOT = {'legs': {FL: {'origin': (-0.5, 0.5, 0),
                       'lengths': {'ao': 135.0, 'bo': 120.0, 'bcx': 290.0, 'bcy': 60.0,
                                   'ae': 500.0, 'de': 100.0, 'ef': 450.0, 'fg': 300.0,
                                   'fh': 200.0, 'gi': 520.0, 'bf': 600.0, 'gj': 1055.0,
                                   'yaw_a': 700.0, 'yaw_b': 50.0, 'yaw_c': 280.0},
                       'verins': [535, 615, 520],
                       'matrix': np.array([[1, 0, 0],
                                           [0, 1, 0],
                                           [0, 0, 1]]),
                       'og': 1,
                       'pos_abs': np.array([0.0, 0.0, 0.0])},
                  FR: {'origin': (0.5, 0.5, 0),
                       'lengths': {'ao': 130.0, 'bo': 120.0, 'bcx': 300.0, 'bcy': 60.0,
                                   'ae': 500.0, 'de': 100.0, 'ef': 445.0, 'fg': 285.0,
                                   'fh': 200.0, 'gi': 500.0, 'bf': 603.0, 'gj': 1035.0,
                                   'yaw_a': 700.0, 'yaw_b': 55.0, 'yaw_c': 280.0},
                       'verins': [535, 615, 520],
                       'matrix': np.array([[0, -1, 0],
                                           [1, 0, 0],
                                           [0, 0, 1]]),
                       'og': 1,
                       'pos_abs': np.array([0.0, 0.0, 0.0])},
                  RL: {'origin': (-0.5, -0.5, 0),
                       'lengths': {'ao': 130.0, 'bo': 120.0, 'bcx': 295.0, 'bcy': 60.0,
                                   'ae': 495.0, 'de': 100.0, 'ef': 450.0, 'fg': 300.0,
                                   'fh': 200.0, 'gi': 515.0, 'bf': 600.0, 'gj': 1055.0,
                                   'yaw_a': 700.0, 'yaw_b': 60.0, 'yaw_c': 280.0},
                       'verins': [535, 615, 520],
                       'matrix': np.array([[0, 1, 0],
                                           [-1, 0, 0],
                                           [0, 0, 1]]),
                       'og': 1,
                       'pos_abs': np.array([0.0, 0.0, 0.0])},
                  RR: {'origin': (0.5, -0.5, 0),
                       'lengths': {'ao': 130.0, 'bo': 120.0, 'bcx': 290.0, 'bcy': 60.0,
                                   'ae': 495.0, 'de': 100.0, 'ef': 445.0, 'fg': 300.0,
                                   'fh': 200.0, 'gi': 500.0, 'bf': 600.0, 'gj': 1045.0,
                                   'yaw_a': 700.0, 'yaw_b': 55.0, 'yaw_c': 280.0},
                       'verins': [535, 615, 520],
                       'matrix': np.array([[-1, 0, 0],
                                           [0, -1, 0],
                                           [0, 0, 1]]),
                       'og': 1,
                       'pos_abs': np.array([0.0, 0.0, 0.0])}
                  },

         'body': {'offset': [500, 500, 0],
                  'center': [0, 0, 0],
                  'omega': {'l': 0, 'm': 0, 'n': 0},
                  'com': [0, 0, 0]},

         'idle_pos': {'verins': [535, 615, 520]}}


######################### ACCES AU ROBOT ##########################

def get_verins_3(leg_id):
    """
    Retourne les valeurs des 3 vérins pour une jambe

    :param leg_id: ID de la patte
    """
    return ROBOT['legs'][leg_id]['verins']


def set_verins_3(v1, v2, v3, leg_id):
    """
    actualise les valeurs de vérins courantes stockées dans la structure ROBOT avec celles entrées en paramètre.

    :param v1: élongation de v1
    :param v2: élongation de v2
    :param v3: élongation de v3
    :param leg_id: ID de la patte
    """
    global ROBOT
    ROBOT['legs'][leg_id]['verins'][0] = v1
    ROBOT['legs'][leg_id]['verins'][1] = v2
    ROBOT['legs'][leg_id]['verins'][2] = v3
    return


def get_verins_12():
    """retourne les valeurs des 12 vérins"""
    res = []
    for i in range(0, 12, 3):
        res.append(ROBOT['legs'][i // 3]['verins'][0])
        res.append(ROBOT['legs'][i // 3]['verins'][1])
        res.append(ROBOT['legs'][i // 3]['verins'][2])
    return res


def set_verins_12(V):
    """
    actualise les 12 vérins

    :param V: Liste des élongations des 12 vérins (en mm)
    """
    for i in range(0, len(V), 3):
        set_verins_3(V[i], V[i + 1], V[i + 2], i / 3)


def get_O():
    """retourne les coordonnées du centre du robot"""
    return ROBOT['body']['center']


def set_O(O):
    """
    actualise les coordonnées du centre du robot

    :param O: coordonnées de O
    """
    global ROBOT
    ROBOT['body']['center'] = O


def get_omega():
    """Retourne les angles que forme le châssis par rapport au repère fixe"""
    return np.array([ROBOT['body']['omega']['l'], ROBOT['body']['omega']['m'], ROBOT['body']['omega']['n']])


def set_omega(Omega):
    """
    Actualise les valeurs des angles du robot par rapport au repère absolu

    :param Omega: Valeur des angles
    """
    global ROBOT
    for i in range(3):
        ROBOT['body']['omega'][i] = Omega[i]


def get_leg_pos(leg_id):
    """
    Renvoie la position absolue d'une patte

    :param leg_id: ID de la patte
    """
    return ROBOT['legs'][leg_id]['pos_abs']


def set_leg_pos(pos, leg_id):
    """
    Actualise la position absolue d'une patte

    :param pos: coordonnées du bout de la patte
    :param leg_id: ID de la patte
    """
    global ROBOT
    ROBOT['legs'][leg_id]['pos_abs'] = np.array(pos)


def get_X():
    """Renvoie la position de la position"""
    X = []
    for i in range(4):
        X = np.append(X, get_leg_pos(i))
    return X


def set_X(X):
    """
    Actualise la valeur de la position

    :param X: coordonnées des 12 pattes
    """
    global ROBOT
    for i in range(4):
        set_leg_pos(X[3 * i:3 * i + 3], i)


def get_og(leg_id):
    """
    Renvoie 1 si la patte est au sol, 0 sinon

    :param leg_id: ID de la patte
    """
    return ROBOT['legs'][leg_id]['og']


def set_og(bool, leg_id):
    """
    Actualise la valeur de OnGround pour la patte leg_id dans ROBOT

    :param bool: booléen
    :param leg_id: ID de la patte
    """
    global ROBOT
    ROBOT['legs'][leg_id]['og'] = bool


def get_com():
    """Retourne de centre de masse du robot"""
    return ROBOT['body']['com']


def set_com(com):
    """
    actualise la position du centre de masse du robot

    :param com: coordonnées du COM
    """
    global ROBOT
    ROBOT['body']['com'] = com


########################### ROTATIONS #############################

def gen_matrix_rota(angle, axis):
    if axis == 0:  # Rotation selon x
        return np.array([[1, 0, 0],
                         [0, np.cos(angle), -np.sin(angle)],
                         [0, np.sin(angle), np.cos(angle)]])
    elif axis == 1:  # Rotation selon y
        return np.array([[np.cos(angle), 0, -np.sin(angle)],
                         [0, 1, 0],
                         [np.sin(angle), 0, np.cos(angle)]])
    else:  # Rotation selon z
        return np.array([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle), np.cos(angle), 0],
                         [0, 0, 1]])


def gen_deriv_matrix_rota(angle, axis):
    if axis == 0:  # Rotation selon x
        return np.array([[0, 0, 0],
                         [0, -np.sin(angle), -np.cos(angle)],
                         [0, np.cos(angle), -np.sin(angle)]])
    elif axis == 1:  # Rotation selon y
        return np.array([[-np.sin(angle), 0, -np.cos(angle)],
                         [0, 0, 0],
                         [np.cos(angle), 0, -np.sin(angle)]])
    else:  # Rotation selon z
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


###################### CHGT DE REFERENTIEL ########################

def robot_ref_to_leg(point, leg_id):
    """
    Calcul les coordonnées de point en effectuant un changement de référentiel, du robot à celui de la patte

    :param point: coordonées du point
    :param leg_id: ID de la patte
    :return: coords dans le ref de la patte
    """
    point = ROBOT['legs'][leg_id]['matrix'] @ point
    point = point - ROBOT['body']['offset']
    return point


def leg_ref_to_robot(point, leg_id):
    """
    Passer d'un point dans le référentiel de la jambe à celui du robot

    :param point: coordonnées du point
    :param leg_id: ID de la patte
    """
    new_point = [point[0] + ROBOT['body']['offset'][0],
                 point[1] + ROBOT['body']['offset'][1],
                 point[2] + ROBOT['body']['offset'][2]]
    return ROBOT['legs'][leg_id]['matrix'].T @ new_point


def robot_ref_to_abs(point, O, Omega):
    return O + gen_R(Omega[0], Omega[1], Omega[2]) @ point


def d3_to_d2(x, y, z):
    """
    Retourne (X, Z, alpha) les coordonnées cylindriques du bout de la patte

    :param x: coordonnée en x
    :param y: coordonnée en y
    :param z: coordonnée en z
    """
    X = np.sqrt(x ** 2 + y ** 2)
    Z = z
    calpha = y / X
    return X, Z, calpha


def d2_to_d3(X, Z, calpha):
    """
    Retourne (x, y, z) les coordonnées cartésiennes du bout de la patte

    :param X: coordonnée en X
    :param Z: coordonnée en Z
    :param calpha: valeur de l'angle alpha
    """
    x = X * np.sqrt(1 - calpha ** 2)
    y = X * calpha
    z = Z
    return x, y, z


############################ GENERAL ##############################

def cos_angle_to_v3(cangle, lpl):
    """
    Fonction auxiliaire de normalized_move_xyz
    Retourne l'élongation de v3 en mm en fonction du cosinus de l'angle de la patte au chassis

    # >>> cos_angle_to_v3(np.cos(np.pi/4), ROBOT['legs'][0]['lengths']) - al_kashi_longueur(KO, LO, np.pi/4 - np.arccos(MO/LO))
    # 0.0
    # >>> v3_to_cos_angle(cos_angle_to_v3(np.cos(np.pi/4), ROBOT['legs'][0]['lengths']), ROBOT['legs'][0]['lengths']) - np.cos(np.pi/4) < 0.0000001
    # True
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

    # >>> v3_to_cos_angle(500, ROBOT['legs'][0]['lengths']) - np.cos(al_kashi_angle(LO, KO, 500) + np.arccos(MO/LO))
    # 0.0
    # >>> cos_angle_to_v3(v3_to_cos_angle(500, ROBOT['legs'][0]['lengths']), ROBOT['legs'][0]['lengths']) - 500 < 0.0000001
    # True
    """
    KO = lpl['yaw_c']
    LM = lpl['yaw_b']
    MO = lpl['yaw_a']
    LO = np.sqrt(LM ** 2 + MO ** 2)
    return (KO ** 2 + LO ** 2 - v3 ** 2) / (2 * KO * LO) * MO / LO - np.sqrt(
        1 - ((KO ** 2 + LO ** 2 - v3 ** 2) / (2 * KO * LO)) ** 2) * np.sqrt(1 - (MO / LO) ** 2)


def distance(pt1, pt2=(0, 0, 0)):
    """
    Calcule la distance à l'origine ou entre 2 points dans l'espace en 2 ou 3 dimensions
    (en fonction de la dimension des points donnés en paramètre et de leur nombre)

    :param pt1: 1er point
    :param pt2: 2nd point (optionnel)

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
    if len(pt1) == 2:
        return np.sqrt((pt1[0] - pt2[0]) ** 2 + (pt1[1] - pt2[1]) ** 2)
    return np.sqrt((pt1[0] - pt2[0]) ** 2 + (pt1[1] - pt2[1]) ** 2 + (pt1[2] - pt2[2]) ** 2)


def al_kashi_longueur(a, b, alpha):
    return np.sqrt(a ** 2 + b ** 2 - 2 * a * b * np.cos(alpha))


def al_kashi_angle(a, b, c):
    """ Retourne l'angle calculé par la formule d'al kashi dans le triangle de cotés de longueurs a, b et c """
    return np.arccos((a ** 2 + b ** 2 - c ** 2) / (2 * a * b))


def normal_vector(v):
    """
    compute the normal unitary vector

    :param v: vector
    :return: normal unitary vector
    >>> normal_vector((1, 0))
    (0, 1)
    >>> normal_vector((-1, 0))
    (0, -1)
    >>> normal_vector((0, 1))
    (-1, 0)
    >>> normal_vector((0, -1))
    (1, 0)
    >>> normal_vector((0.3, -0.6))
    (0.6, 0.3)
    """
    return -1 * v[1], v[0]


def norm2(v):
    """
    calcul la norme 2 du vecteur en paramètre

    :param v: vecteur
    :return: norme 2
    """
    n=0
    for x in v:
        n+=x**2
    return n


def norm(v):
    """
    calcul la norme d'un vecteur

    :param v: vecteur
    :return: norme
    """
    return np.sqrt(norm2(v))


def unitary_vec(v):
    """
    Calcul le vecteur unitaire/normalisé du vecteur en paramètre

    :param v: vecteur
    :return: vecteur unitaire
    """
    n=norm(v)
    assert(n!=0)
    return tuple([a/n for a in v])

############################################################################

if __name__ == "__main__":
    import doctest
    doctest.testmod()