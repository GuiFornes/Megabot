# L'ensemble des distances sont exprimées en mm : segments de patte et élongations des vérins

import numpy as np
from numpy.linalg import inv
from qpsolvers import solve_qp
import kinetic as kin
from upg_tools import *
from upg_jacobian import *
# from com import wait_move, tell_controlers

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

######################### Tools for 'Legs' struct ###########################

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

############################### DIRECT #####################################

def direct_18(V):
    return 0

def direct_12(V):
    """
    Retourne les positions des extrémités des 4 pattes correspondant aux élongations V des vérins
    Les distances sont exprimées en mm et les coordonnées sont exprimées dans le référentiel du robot
    """
    R = []
    for i in range(4):
        R = np.append(R, direct_robot(V[i*3], V[i*3 + 1], V[i*3 + 2], i))
    return R

def direct_robot(v1, v2, v3, leg_id):
    """
    Retourne les positions x, y, z du bout de la patte en fonctions de v1, v2, v3 dans le référentiel du robot
    Se base sur le modèle direct de Julien
    """
    lpl = LEGS[leg_id]['lengths']
    X, Z = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)['J']
    calpha = v3_to_cos_angle(v3, lpl)
    Pos = np.array([X * np.sqrt(1 - calpha ** 2) * 1000, X * calpha * 1000, Z * 1000])
    return MR[leg_id] @ (Pos + L)

def direct_leg(v1, v2, v3):
    """
    Retourne les positions x, y, z du bout de la patte en fonctions de v1, v2, v3 dans le référentiel de la patte FL
    Se base sur le modèle direct de Julien
    """
    lpl = LEGS[kin.FL]['lengths']
    X, Z = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)['J']
    calpha = v3_to_cos_angle(v3, lpl)
    x = X * np.sqrt(1 - calpha ** 2) * 1000
    y = X * calpha * 1000
    z = Z * 1000
    return x, y, z

############################## INDIRECT ####################################

# def move(traj_center, V):
#     V0 = V
#     L = [V0]
#     X0 = direct_12(V0)

def move_12(traj, V, solved=True):
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
        X0 = direct_12(V0)
        dX = traj[i] - X0
        J = gen_jacob_12(V)
        if solved:  # Utilisation du solveur
            P = inv(J).T @ inv(J)
            q = - inv(J).T @ dX
            lb = np.array([450.0, 450.0, 450.0,
                           450.0, 450.0, 450.0,
                           450.0, 450.0, 450.0,
                           450.0, 450.0, 450.0]) - V0
            ub = np.array([650.0, 650.0, 650.0,
                           650.0, 650.0, 650.0,
                           650.0, 650.0, 650.0,
                           650.0, 650.0, 650.0]) - V0
            dV = solve_qp(P, q, lb=lb, ub=ub)
        else:  # Utilisation de la jacobienne sans solveur
            dV = J @ dX
        V0 = V0 + dV
        set_verins_12(V)
        for v in V0: assert 449.9 < v < 650.1, 'Elongation de vérin invalide'
        R.append(V0)
    return R

def move_leg(traj, v1, v2, v3, leg_id, display=False, upgrade=False, solved=True):
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
        lpl = LEGS[leg_id]['lengths']
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
        if solved:  # Utilisation du solveur
            P = inv(J).T @ inv(J)
            q = - inv(J).T @ dX
            lb = np.array([450.0 - v1, 450.0 - v2, 450.0 - v3])
            ub = np.array([650.0 - v1, 650.0 - v2, 650.0 - v3])
            dV = solve_qp(P, q, lb=lb, ub=ub)
        else:  # Utilisation de la jacobienne sans solveur
            dV = J @ dX
        v1, v2, v3 = v1 + dV[0], v2 + dV[1], v3 + dV[2]
        set_verins_3(v1, v2, v3, leg_id)
        R.append((v1, v2, v3))
    return R

################################ MOVE ######################################

def draw_circle_12(n, r, V):
    """
    retourne une trajectoire circulaire pour les 4 pattes
    @param n:
    @param r:
    @param V:
    @return:
    """
    traj_FL = draw_circle(r, n, V[0], V[1], V[2], 0)
    traj_FR = draw_circle(r, n, V[3], V[4], V[5], 1)
    traj_RL = draw_circle(r, n, V[6], V[7], V[8], 2)
    traj_RR = draw_circle(r, n, V[9], V[10], V[11], 3)
    traj = []
    for i in range(n):
        t = traj_FL[i]
        t = np.append(t, MR[1] @ traj_FR[i])
        t = np.append(t, MR[2] @ traj_RL[i])
        t = np.append(t, MR[3] @ traj_RR[i])
        traj.append(t)
    return traj

def draw_circle(r, n, v1, v2, v3, leg_id):
    """
    retourne une trajectoire circulaire pour une patte
    @param r:
    @param n:
    @param v1:
    @param v2:
    @param v3:
    @param leg_id:
    @return:
    """
    lpl = LEGS[leg_id]['lengths']
    pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
    X0, Z0 = pts['J'][0] * 1000, pts['J'][1] * 1000
    x0, y0, z0 = d2_to_d3(X0, Z0, v3_to_cos_angle(v3, lpl))
    R = []
    for k in range(n + 1):
        R.append(np.array([x0+ r * np.cos(2 * k * np.pi / n) - r,
                           y0 + r * np.sin(2 * k * np.pi / n),
                           z0]) + L)
    return R

def draw_line_3(v1, v2, v3, dx, dy, dz, n, leg_id):
    """
    retourne trajectoire rectiligne en liste de point pour u certain décalage et une certaine patte
    @param v1:
    @param v2:
    @param v3:
    @param dx:
    @param dy:
    @param dz:
    @param n:
    @param leg_id:
    @return:
    """
    x0, y0, z0 = direct_leg(v1, v2, v3)
    traj = []
    for i in range(n):
        traj.append((x0 + i * dx / n, y0 + i * dy / n, z0 + i * dz / n))
    return traj

def draw_line_12(V, dx, dy, dz, n):
    """
    retourne la discrétisation d'une trajectoire en ligne droite d'un certain décalage pour les 4 pattes à la fois
    @param V: liste des 12 vérins
    @param dx: décalage en x
    @param dy: décalage en y
    @param dz: décalage en z
    @param n: nombre d'étapes dans la traj
    @return: la trajextoire
    """
    traj_FL = draw_line_3(V[0], V[1], V[2], 15, 20, 10, 10, 0)
    traj_FR = draw_line_3(V[3], V[4], V[5], 15, 20, 10, 10, 1)
    traj_RL = draw_line_3(V[6], V[7], V[8], 15, 20, 10, 10, 2)
    traj_RR = draw_line_3(V[9], V[10], V[11], 15, 20, 10, 10, 3)
    traj = []
    for i in range(n):
        t = traj_FL[i]
        t = np.append(t, MR[1] @ traj_FR[i])
        t = np.append(t, MR[2] @ traj_RL[i])
        t = np.append(t, MR[3] @ traj_RR[i])
        traj.append(t)
    return traj
############################################################################

def upg_inverse_kinetic_robot_ref(legs,leg_id,point):
    """
    fonction de test de notre cinématique inverse de la manière de l'implémentation de julien
    dans le but de comparer sur le simulateur. A noter que ce n'est pas dutout une utilisation optimale car on utilise
    pas ici le potentiel de travailler sur la division d'un mouvements en petits écarts.
    @param legs: structure legs inutile ici car on utilise LEGS (qui est global dans ce fichier)
    @param leg_id: id de la jambe (entre 0 et 3)
    @param point: point cible pour le bout de la jambe en question
    @return: les positions de vérins correspondantes
    """
    lpl = LEGS[leg_id]['lengths']
    V = get_verins_3(leg_id)
    print(point)
    x0, y0, z0 = direct_leg(V[0], V[1], V[2])
    print("x0 =", x0, y0, z0)
    # dX = point['x'] - LEGS[leg_id]['origin']['x']
    # dY = point['y'] - LEGS[leg_id]['origin']['y']
    # calpha = v3_to_cos_angle(V[2], lpl)
    # dx, dy, dz = d2_to_d3(dX, dY, calpha)
    dx = point[0]*1000 - x0
    dy = point[1]*1000 - y0
    dz = point[2]*1000 - z0
    print(dx, dy, dz)
    traj = []
    for i in range(10):
        traj.append(np.array([x0 + i * dx/10,
                              y0 + i * dy/10,
                              z0 + i * dz/10]))
    Verins = move_leg(traj, V[0], V[1], V[2], leg_id, upgrade=True, solved=True)
    print(Verins)
    res = [Verins[9][0]/1000, Verins[9][1]/1000, Verins[9][2]/1000]
    return False, res

def upg_init_legs(controlers):
    """
    met à jour la structure LEGS avec les positions des vérins
    @param controlers: structure controllers du fichier static_walk
    """
    print(controlers[FL].la, LEGS[FL]['verins'])
    for l in ALL_LEGS:
        LEGS[l]['verins'][0] = 450 + controlers[l].la[0]['position']
        LEGS[l]['verins'][1] = 450 + controlers[l].la[1]['position']
        LEGS[l]['verins'][2] = 450 + controlers[l].la[2]['position']
    print(LEGS[FL]['verins'])

def do_the_traj():
    """
    créée une trajectoire, calcul les verins correspondant
    puis envoie l'information sur leurs élongations (qui est interceptée par la simulation
    """
    print("ok1")
    V = get_verins_12()
    # for i in ALL_LEGS:
        # wait_move(i, 2)
    #traj = draw_line_12(V, 15, 20, 10, 10)
    traj = draw_circle_12(20, 100, V)
    V = move_12(traj, V, True)
    print(traj[0] ,V[0])
    # for i in range(len(traj)):
        # tell_controlers(V[i])
        # wait_move(list(range(4)), 0.1)
    V = get_verins_12()
    set_verins_12(V)
    return

############################################################################
if __name__ == "__main__":
    import doctest

    doctest.testmod()