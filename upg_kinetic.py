# L'ensemble des distances sont exprimées en mm : segments de patte et élongations des vérins

import numpy as np
from numpy.linalg import inv
from qpsolvers import solve_qp
import kinetic as kin
from upg_tools import *
from upg_jacobian import *

############################### DIRECT #####################################

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
    lpl = kin.LEGS[leg_id]['lengths']
    X, Z = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)['J']
    calpha = v3_to_cos_angle(v3, lpl)
    Pos = np.array([X * np.sqrt(1 - calpha ** 2) * 1000, X * calpha * 1000, Z * 1000])
    return MR[leg_id] @ (Pos + L)

def direct_leg(v1, v2, v3):
    """
    Retourne les positions x, y, z du bout de la patte en fonctions de v1, v2, v3 dans le référentiel de la patte FL
    Se base sur le modèle direct de Julien
    """
    lpl = kin.LEGS[kin.FL]['lengths']
    X, Z = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)['J']
    calpha = v3_to_cos_angle(v3, lpl)
    x = X * np.sqrt(1 - calpha ** 2) * 1000
    y = X * calpha * 1000
    z = Z * 1000
    return x, y, z

############################## INDIRECT ####################################

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
        lpl = kin.LEGS[leg_id]['lengths']
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
        R.append((v1, v2, v3))
    return R

################################ MOVE ######################################

def draw_circle_12(n, r, V):
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
    lpl = kin.LEGS[leg_id]['lengths']
    pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
    X0, Z0 = pts['J'][0] * 1000, pts['J'][1] * 1000
    x0, y0, z0 = d2_to_d3(X0, Z0, v3_to_cos_angle(v3, lpl))
    R = []
    for k in range(n + 1):
        R.append(np.array([x0+ r * np.cos(2 * k * np.pi / n) - r,
                           y0 + r * np.sin(2 * k * np.pi / n),
                           z0]) + L)
    return R

############################################################################

if __name__ == "__main__":
    import doctest

    doctest.testmod()
