# L'ensemble des distances sont exprimées en mm : segments de patte et élongations des vérins
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import inv
from qpsolvers import solve_qp
from upg_tools import *
from upg_jacobian import *
# from com import wait_move, tell_controlers


############################### DIRECT #####################################

def direct_O(X_abs, V):
    """
    cinématique direct dans le référentiel absolue (base flottante)
    :param X_abs: x absolu
    :param V: liste des 12 élongations
    :return: ????
    """
    if LEGS[0]['og'] and LEGS[2]['og']:
        O0 = X_abs[0 : 3] - direct_leg(V[0], V[1], V[2])
        O2 = MR[0] @ X_abs[6 : 9] - direct_leg(V[6], V[7], V[8])
        return 0.5 * (O0 + O2)
    O1 = MR[1] @ X_abs[3 : 6] - direct_leg(V[3], V[4], V[5])
    O3 = MR[3] @ X_abs[9 : 12] - direct_leg(V[9], V[10], V[11])
    return 0.5 * (O1 + O3)

def direct_12(V):
    """
    Retourne les positions des extrémités des 4 pattes correspondant aux élongations V des vérins
    Les distances sont exprimées en mm et les coordonnées sont exprimées dans le référentiel du robot
    :param V: liste des 12 élongations des vérins
    :return: les positions des 4 extrémités de patte
    """
    R = []
    for i in range(4):
        R = np.append(R, direct_robot(V[i*3], V[i*3 + 1], V[i*3 + 2], i))
    return R

def direct_robot(v1, v2, v3, leg_id):
    """
    Retourne les positions x, y, z du bout de la patte en fonctions de v1, v2, v3 dans le référentiel du robot
    Se base sur le modèle direct de Julien
    :param v1: élongation de v1
    :param v2: élongation de v2
    :param v3: élongation de v3
    :param leg_id: id de l patte
    :return: coordonnées du bout de la patte
    """
    lpl = LEGS[leg_id]['lengths']
    X, Z = get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)['J']
    calpha = v3_to_cos_angle(v3, lpl)
    Pos = np.array([X * np.sqrt(1 - calpha ** 2) * 1000, X * calpha * 1000, Z * 1000])
    return MR[leg_id] @ (Pos + L)

def direct_leg(v1, v2, v3):
    """
    Retourne les positions x, y, z du bout de la patte en fonctions de v1, v2, v3 dans le référentiel de la patte FL
    Se base sur le modèle direct de Julien
    :param v1: élongation de v1
    :param v2: élongation de v2
    :param v3: élongation de v3
    :return: coord x, y et z du bout de la patte
    """
    lpl = LEGS[FL]['lengths']
    X, Z = get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)['J']
    calpha = v3_to_cos_angle(v3, lpl)
    x = X * np.sqrt(1 - calpha ** 2) * 1000
    y = X * calpha * 1000
    z = Z * 1000
    return x, y, z

############################## INDIRECT ####################################

def move(traj_center, X, Omega, V, angle_chassis=10):
    V0 = V
    X0_abs = X
    L = [V0]
    for i in range(1, len(traj_center)):
        O0 = direct_O(X0_abs, V)
        dO = traj_center[i] - O0
        X0_rel = direct_12(V0)
        M = gen_jacob_alt(V0, X0_rel, Omega[0], Omega[1], Omega[2])
        P = M.T @ M
        q = - M.T @ dO

        # Contraintes sur les pattes
        b = np.zeros(30)
        A = np.zeros((30, 30))
        for j in range(4):
            if LEGS[i]['og']:
                for k in range(3):
                    A[j * 3 + k][j * 3 + k] = 1
        for j in range(3):
            A[27 + j][27 + j] = 1
        lb = np.full(30, -np.inf)
        ub = np.full(30, np.inf)
        # Contraintes sur les vérins
        for j in range(12):
            lb[15 + j] = 450.0 - V0[j]
            ub[15 + j] = 650.0 - V0[j]
        # Contraintes sur les angles
        for j in range(3):
            lb[12 + j] = -angle_chassis
            ub[12 + j] = angle_chassis

        sol = solve_qp(P, q, lb=lb, ub=ub, A=A, b=b)
        X0_abs = X0_abs + sol[0:12]
        V0 = V0 + sol[15:27]
        set_verins_12(V)
        for v in V0: assert 449.9 < v < 650.1, 'Elongation de vérin invalide'
        L.append(V0)
    return L

def move_on_floor(traj_center, max_angle=10):
    """ On se considère en (O,O) du référentiel absolu au début du mouvement"""
    V = get_verins_12()
    Omega = get_omega()
    X_abs = direct_12(V)
    for i in range(4):
        X_abs[3*i + 2] = 0
    L = [V]
    for i in range(1, len(traj_center)):
        O0 = direct_O(X_abs, V)
        dO = np.concatenate((np.zeros(12), traj_center[i] - O0))
        M = gen_new_jacob(V, Omega, direct_12(V))
        P = M.T @ M
        print(P)
        print("##############################################")
        print(inv(P))
        q = - M.T @ dO

        # Contraintes sur les vérins
        lb = np.zeros(15)
        ub = np.zeros(15)
        for j in range(12):
            lb[j] = 450.0 - V[j]
            ub[j] = 650.0 - V[j]

        # Contraintes sur les angles
        for j in range(12, 15):
            lb[j] = - max_angle
            ub[j] = max_angle

        sol = solve_qp(P, q, lb=lb, ub=ub)
        V = V + sol[0:12]
        set_verins_12(V)
        for v in V: assert 449.9 < v < 650.1, 'Elongation de vérin invalide'
        L.append(V)
    return L

def move_rel(traj, V, solved=True):
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
        pts = get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
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
    @param n: précision de la discrétisation (nbr de pts)
    @param r: rayon des cercles
    @param V: liste des 12 élongations de vérins
    @return: liste des pts intermédiaires des 4 trajectoires
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
    @param r: rayon
    @param n: précision de la discrétisation (nbr de points)
    @param v1: élongation de v1
    @param v2: élongation de v2
    @param v3: élongation de v3
    @param leg_id: id de la jambe
    @return: trajectoire en liste de cooredonnées intermédiaires
    """
    lpl = LEGS[leg_id]['lengths']
    pts = get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
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
    retourne trajectoire rectiligne en liste de point pour un certain décalage et une certaine patte
    @param v1: elongation de v1
    @param v2: elongation de v2
    @param v3: elongation de v3
    @param dx: décalage en x
    @param dy: décalage en y
    @param dz: décalage en z
    @param n: précision de la discrétisation (nombre de points intermédiaires)
    @param leg_id: id de la jambe
    @return: trajectoire en liste coordonnées intermédiaire
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
    fonction de test de notre cinématique inverse de la manière de l'implémentation de julien,
    dans le but de comparer sur le simulateur. A noter que ce n'est pas dutout une utilisation optimale
    car on utilise pas ici le potentiel de travailler sur la division d'un mouvements en petits écarts,
    étant donné le fait que l'on renvoie que la position finale.

    Ne MArCHE PAS

    @param legs: structure legs inutile ici car on utilise LEGS (qui est global dans ce fichier)
    @param leg_id: id de la jambe (entre 0 et 3)
    @param point: point cible pour le bout de la jambe en question
    @return: les positions de vérins correspondantes
    """

    lpl = LEGS[leg_id]['lengths']
    V = get_verins_3(leg_id)
    print("pt = ", point, "   leg_id = " ,leg_id)
    x0, y0, z0 = direct_leg(V[0], V[1], V[2])
    print("x0 =", x0, y0, z0)
    # dX = point['x'] - LEGS[leg_id]['origin']['x']
    # dY = point['y'] - LEGS[leg_id]['origin']['y']
    # calpha = v3_to_cos_angle(V[2], lpl)
    # dx, dy, dz = d2_to_d3(dX, dY, calpha)
    point *= np.full(3, 1000)
    leg_ref_point = robot_ref_to_leg(point, leg_id)
    print("leg_ref =", leg_ref_point)
    dx = leg_ref_point[0] - x0
    dy = leg_ref_point[1] - y0
    dz = leg_ref_point[2] - z0
    print("dX = ", dx, dy, dz)
    traj = []
    for i in range(10):
        traj.append(np.array([x0 + i * dx/10,
                              y0 + i * dy/10,
                              z0 + i * dz/10]))
    Verins = move_leg(traj, V[0], V[1], V[2], leg_id, upgrade=True, solved=True)
    res = [Verins[9][0]/1000, Verins[9][1]/1000, Verins[9][2]/1000]
    error = False
    for i in range(3):
        if res[i] == 650 or res[i] == 450:
            error = True
    return error, res

def upg_init_legs(controlers):
    """
    met à jour la structure LEGS avec les positions des vérins
    @param controlers: structure controllers du fichier static_walk
    """
    # print(controlers[FL].la, LEGS[FL]['verins'])
    global LEGS
    for l in ALL_LEGS:
        LEGS[l]['verins'][0] = 450 + controlers[l].la[1]['position']
        # print("UPG !", controlers[l].la[0]['position'], controlers[l].la[1]['position'], controlers[l].la[2]['position'])
        LEGS[l]['verins'][1] = 450 + controlers[l].la[2]['position']
        LEGS[l]['verins'][2] = 450 + controlers[l].la[0]['position']
    print(LEGS[FL]['verins'], "ici")

def get_the_traj():
    """
    créée une trajectoire, calcul les verins correspondant
    puis envoie l'information sur leurs élongations
    """
    V = get_verins_12()
    traj = draw_circle_12(20, 100, V)
    return move_rel(traj, V, True)

def shake_dat_ass_abs(n, amp):
    """
    retourne les élongations de vérins nécéssaires à la réalisation d'un mouvement de haut en bas,
    en travaillant dans le référentiel absolu.
    :param n: précision de la discrétisation (nombre de points)
    :param amp: amplitude du mouvement
    :return: liste de tableau de 12 vérins pour chaques positions intermédiaires du mouvement
    """
    z0 = - get_leg_points_V1_V2(get_verins_12()[0] / 1000, get_verins_12()[1] / 1000, LEGS[0]['lengths'])['J'][1]
    L = np.linspace(-1, 1, n)
    traj_center = [[0, 0, z0 - amp * np.sin(l)] for l in L]
    return move_on_floor(traj_center)

def shake_dat_ass_rel(n, amp):
    """
    retourne les élongations de vérins nécéssaires à la réalisation d'un mouvement de haut en bas,
    en travaillant dans le référentiel relatif.
    :param n: précision de la discrétisation (nombre de points)
    :param amp: amplitude du mouvement
    :return: liste de tableau de 12 vérins pour chaques positions intermédiaires du mouvement
    """
    V = get_verins_12()
    X = direct_12(V)
    L = np.linspace(0, np.pi, n)
    traj = [[X[0], X[1], X[2] + amp * np.sin(l),
                    X[3], X[4], X[5] + amp * np.sin(l),
                    X[6], X[7], X[8] + amp * np.sin(l),
                    X[9], X[10], X[11] + amp * np.sin(l)] for l in L]
    return move_rel(traj, V)

############################## ROTATION ####################################

def min_diff(square=True):
    """
    Retourne la patte la plus proche de sa position de repos en faisant la moyenne des écarts d'élongation
    des vérins ou bien la moyenne des écarts d'élongation au carré des vérins
    """
    V = get_verins_12()
    max_diff = np.inf
    best_leg = 0
    for i in range(4):
        current_diff = 0
        for j in range(3):
            if square:
                current_diff += (V[3 * i + j] - IDLE_POSITION['verins'][j]) ** 2
            else:
                current_diff += V[3 * i + j] - IDLE_POSITION['verins'][j]
        if current_diff < max_diff :
            best_leg = i
            max_diff = current_diff
    return best_leg

def traj_rota():
    leg_ig = min_diff()
    x, y, z = direct_robot(LEGS[leg_ig]['verrins'][0], LEGS[leg_ig]['verrins'][1], LEGS[leg_ig]['verrins'][2], leg_ig)




############################################################################
if __name__ == "__main__":
    import doctest

    doctest.testmod()