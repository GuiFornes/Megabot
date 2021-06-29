# L'ensemble des distances sont exprimées en mm : segments de patte et élongations des vérins

import numpy as np
from numpy.linalg import inv
from qpsolvers import solve_qp
from kinetic import get_leg_points_V1_V2


FL=0 # front left leg
FR=1 # front right leg
RL=2 # rear left leg
RR=3 # rear right leg

ALL_LEGS=(FL,FR,RL,RR)

BODY_FRAME=1.0

LEGS={FL:{'origin':(-BODY_FRAME/2.0,BODY_FRAME/2.0,0),
          'lengths':{'ao':135.0,'bo':120.0,'bcx':290.0,'bcy':60.0,
                     'ae':500.0,'de':100.0,'ef':450.0,'fg':300.0,
                     'fh':200.0,'gi':520.0,'bf':600.0,'gj':1055.0,
                     'yaw_a':700.0,'yaw_b':50.0,'yaw_c':280.0},
          'verins':[485, 575, 515],
          'og':1},
      FR:{'origin':(BODY_FRAME/2.0,BODY_FRAME/2.0,0),
          'lengths':{'ao':130.0,'bo':120.0,'bcx':300.0,'bcy':60.0,
                     'ae':500.0,'de':100.0,'ef':445.0,'fg':285.0,
                     'fh':200.0,'gi':500.0,'bf':603.0,'gj':1035.0,
                     'yaw_a':700.0,'yaw_b':55.0,'yaw_c':280.0},
          'verins': [485, 575, 515],
          'og':1},
      RL:{'origin':(-BODY_FRAME/2.0,-BODY_FRAME/2.0,0),
          'lengths':{'ao':130.0,'bo':120.0,'bcx':295.0,'bcy':60.0,
                     'ae':495.0,'de':100.0,'ef':450.0,'fg':300.0,
                     'fh':200.0,'gi':515.0,'bf':600.0,'gj':1055.0,
                     'yaw_a':700.0,'yaw_b':60.0,'yaw_c':280.0},
          'verins': [485, 575, 515],
          'og':1},
      RR:{'origin':(BODY_FRAME/2.0,-BODY_FRAME/2.0,0),
          'lengths':{'ao':130.0,'bo':120.0,'bcx':290.0,'bcy':60.0,
                     'ae':495.0,'de':100.0,'ef':445.0,'fg':300.0,
                     'fh':200.0,'gi':500.0,'bf':600.0,'gj':1045.0,
                     'yaw_a':700.0,'yaw_b':55.0,'yaw_c':280.0},
          'verins': [485, 575, 515],
          'og':1}
      }

######################### Tools for 'Legs' struct ##########################

def set_verins_3(v1, v2, v3, leg_id):
    """
    actualise les valeurs de vérins courantes stockées dans la structure LEGS avec celles entrées en paramètre.
    """
    global LEGS
    LEGS[leg_id]['verins'][0] = v1
    LEGS[leg_id]['verins'][1] = v2
    LEGS[leg_id]['verins'][2] = v3
    return

def set_verins_12(V):
    """actualise les 12 vérins"""
    for i in range(0, len(V), 3):
        set_verins_3(V[i], V[i + 1], V[i + 2], i / 3)

def get_verins_12():
    """retourne les valeurs des 12 vérins"""
    res = []
    for i in range(0, 12, 3):
        res.append(LEGS[i//3]['verins'][0])
        res.append(LEGS[i//3]['verins'][1])
        res.append(LEGS[i//3]['verins'][2])
    return res

def get_verins_3(leg_id):
    """retourne les valeurs des 3 vérins pour une jambe"""
    return LEGS[leg_id]['verins']

def on_ground(leg_id):
    """renseigne le caractère 'au sol' d'une patte"""
    global LEGS
    LEGS[leg_id]['og'] = 1

def stand_up(leg_id):
    """renseigne le caractère levé d'une patte"""
    global LEGS
    LEGS[leg_id]['og'] = 0

def get_og(leg_id):
    return LEGS[leg_id]['og']

x_A = - LEGS[FL]['lengths']['ao']
y_A = 0
x_B = LEGS[FL]['lengths']['bo']
y_B = 0

AB = LEGS[FL]['lengths']['ao'] + LEGS[FL]['lengths']['bo']
AC = np.sqrt((AB + LEGS[FL]['lengths']['bcx']) ** 2 + LEGS[FL]['lengths']['bcy'] ** 2)
BC = np.sqrt(LEGS[FL]['lengths']['bcx'] ** 2 + LEGS[FL]['lengths']['bcy'] ** 2)
AE = LEGS[FL]['lengths']['ae']
AD = AE - LEGS[FL]['lengths']['de']
BF = LEGS[FL]['lengths']['bf']
FH = LEGS[FL]['lengths']['fh']
BH = BF - FH
FG = LEGS[FL]['lengths']['fg']
EF = LEGS[FL]['lengths']['ef']
EG = EF + FG
GI = LEGS[FL]['lengths']['gi']
GJ = LEGS[FL]['lengths']['gj']

KO = LEGS[FL]['lengths']['yaw_c']
LM = LEGS[FL]['lengths']['yaw_b']
MO = LEGS[FL]['lengths']['yaw_a']
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

def robot_ref_to_leg(point, leg_id):
    point = MR[leg_id].T @ point
    point -= L
    return point

def leg_ref_to_robot(point, leg_id):
    point += L
    return MR[leg_id] @ point

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

    >>> cos_angle_to_v3(np.cos(np.pi/4), LEGS[0]['lengths']) - al_kashi_longueur(KO, LO, np.pi/4 - np.arccos(MO/LO))
    0.0
    >>> v3_to_cos_angle(cos_angle_to_v3(np.cos(np.pi/4), LEGS[0]['lengths']), LEGS[0]['lengths']) - np.cos(np.pi/4) < 0.0000001
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

    >>> v3_to_cos_angle(500, LEGS[0]['lengths']) - np.cos(al_kashi_angle(LO, KO, 500) + np.arccos(MO/LO))
    0.0
    >>> cos_angle_to_v3(v3_to_cos_angle(500, LEGS[0]['lengths']), LEGS[0]['lengths']) - 500 < 0.0000001
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