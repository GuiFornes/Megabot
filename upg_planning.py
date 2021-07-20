from upg_kinetic import *
from upg_tools import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

MIN_RADIUS = 2500
MAX_RADIUS = 15000

INF = 100000

TRAJ_ACCUR = 50


def cmd_joystick(d, r):
    """
    Traite les données reçues des commandes joystick pour paramétrer les trajectoires des 4 jambes.

    :param d: vecteur unitaire de direction
    :param r: facteur de rotation (-1 < r < 1)
    :return: centre de rotation, direction de déplacement.
    >>> cmd_joystick((1, 0), 0)
    ((0, 100000), (1, 0))
    >>> cmd_joystick((0, 1), 0)
    ((-100000, 0), (0, 1))
    >>> cmd_joystick((1, 0), 1)
    ((0, 2500), (1, 0))
    >>> cmd_joystick((0, 1), -1)
    ((2500, 0), (0, 1))
    >>> cmd_joystick((0, 0), 0)
    ((None, None), (None, None))
    >>> cmd_joystick((np.sqrt(3)/2, 1/2), -0.5)
    ((4375.0, -7577.722283113838), (0.8660254037844386, 0.5))
    """
    if d[0] == 0 and d[1] == 0:
        if r == 0:
            # aucun mouvement
            c, d = (None, None), (None, None)
        # rotation sur lui-même
        else:
            c = (0, 0)
    elif r == 0:
        # marche tout droit
        n = normal_vector(d)
        c = n[0]*INF, n[1]*INF
    else:
        # marche courbe standard
        n = normal_vector(d)
        signe_r = (r > 0) - (r < 0)
        c = signe_r * (n[0] * MAX_RADIUS - signe_r * r * n[0] * (MAX_RADIUS - MIN_RADIUS)), \
            signe_r * (n[1] * MAX_RADIUS - signe_r * r * n[1] * (MAX_RADIUS - MIN_RADIUS))
    return c, d


############################## RELATIVE WAY ##################################


def compute_traj_form_joystick_rel(joystick):
    """
    A partir des valeurs retournées par cmd_joystick, calcul les trajectoires circulaires de chacunes des pattes

    :param joystick: position du centre de rotation, direction
    :return: les quatres trajectoires
    """
    centre,  direction = joystick[0], joystick[1]
    print(joystick)
    traj = []
    r = []
    V = get_verins_12()
    pos = direct_rel_12(V)
    # Computing radii
    for leg in range(4):
        radius = np.sqrt((pos[leg * 3 + 0] - centre[0])**2 + (pos[leg * 3 + 1] - centre[1])**2)
        r.append(radius)
    r_max = max(r)
    n = int(2 * np.pi * r_max / TRAJ_ACCUR)
    for i in range(n - 10):
        L=[]
        for leg in range(4):
            alpha = np.arccos(abs((centre[1] - pos[3 * leg + 1])) / r[leg])
            signe_cx = (centre[0] - pos[leg * 3 + 0]) / abs(centre[0] - pos[leg * 3 + 0])
            signe_cy = (centre[1] - pos[leg * 3 + 1]) / abs(centre[1] - pos[leg * 3 + 1])
            if signe_cx < 0 and signe_cy < 0:
                alpha = + np.pi/2 - alpha
            if signe_cx > 0 and signe_cy < 0:
                alpha = + np.pi/2 + alpha
            if signe_cx < 0 and signe_cy > 0:
                alpha = - np.pi/2 + alpha
            if signe_cx > 0 and signe_cy > 0:
                alpha = - np.pi/2 - alpha
            L = np.append(L, (r[leg] * np.cos((2*i*np.pi)/n + alpha) + centre[0],
                              r[leg] * np.sin((2*i*np.pi)/n + alpha) + centre[1],
                              pos[3*leg + 2]))
        traj.append(L)
    return traj


def is_accessible_rel(leg_id, point):
    """
    Calcule l'accessibilité d'un point dans le repère 3D de la jambe.

    :param leg_id: ID de la jambe
    :param point: point cible
    :return: True ou False + valeur des verins pour atteindre ce point
    """
    lpl = ROBOT['legs'][leg_id]['lengths']
    v1, v2, v3 = 535, 615, 520
    x0, y0, z0 = direct_rel_3(v1, v2, v3, leg_id)
    xt, yt, zt = point
    dx, dy, dz = xt-x0, yt-y0, zt-z0
    traj = []
    n = 40
    for i in range(n):
        traj.append(np.array([x0 + i * dx / n,
                              y0 + i * dy / n,
                              z0 + i * dz / n]))
    Verins = move_leg(traj, v1, v2, v3, leg_id, display=False, upgrade=False, solved=True)
    res = [Verins[n-1][0], Verins[n-1][1], Verins[n-1][2]]
    acces = True
    xf, yf, zf = direct_rel_3(res[0], res[1], res[2], leg_id)
    # print("x : ", x0, xt, traj[n-1][0], xf)
    # print("y : ", y0, yt, traj[n-1][1], yf)
    # print("z : ", z0, zt, traj[n-1][2], zf)
    # print("V = ", res)
    if distance([xt, yt, zt], [xf, yf, zf]) > 50:
        acces = False
    return acces, res


def furthest_accessible_rel(traj, leg_id):
    """
    Compute the maximum step size following the trajectory, depending on the accessible zone for the leg.
    traj should be the trajectory of all elgs (basically returned by compute_traj_from_joystick)
    :param traj: trajectory to follow
    :param leg_id: ID of the leg
    :return: the furthest point accessible from the traj
    """
    v1, v2, v3 = get_verins_3(leg_id)
    xt, yt, zt = traj[0][leg_id*3 + 0:leg_id * 3 + 3]
    lpl = ROBOT['legs'][leg_id]['lengths']
    for i in range(1, len(traj)):
        # get data
        x0, y0, z0 = direct_rel_3(v1, v2, v3, leg_id)
        if distance((x0, y0, z0), (xt, yt, zt)) > 20:  # exit condition
            return i
        xt, yt, zt = traj[i][leg_id * 3 + 0], traj[i][leg_id * 3 + 1], traj[i][leg_id * 3 + 2]
        dX = np.array([xt - x0, yt - y0, zt - z0])
        # solve
        J = ROBOT['legs'][leg_id]['matrix'].T @ gen_jacob_3(v1 / 1000, v2 / 1000, v3 / 1000, np.arccos(v3_to_cos_angle(v3, lpl)), lpl)
        P = J.T @ J
        q = - J.T @ dX
        lb = np.array([450.0 - v1, 450.0 - v2, 450.0 - v3])
        ub = np.array([650.0 - v1, 650.0 - v2, 650.0 - v3])
        dV = solve_qp(P, q, lb=lb, ub=ub)
        v1, v2, v3 = v1 + dV[0], v2 + dV[1], v3 + dV[2]
    return len(traj)


def furthest_accessible_step_all_legs_rel(traj, step_height):
    V = get_verins_12()
    init = traj[0].copy()
    for j in range(step_height):
        traj.insert(0, init.copy())
        for leg in range(4):
            traj[0][leg * 3 + 2] += step_height + 1 - j
    for k in range(step_height, len(traj)):
        for leg in range(4):
            traj[k][leg * 3 + 2] += step_height

    max_step = []
    for leg in range(4):
        step = furthest_accessible_rel(traj, leg) - (step_height + 1)
        max_step.append(step)
    return max_step


############################### ABSOLUTE WAY ##################################


def compute_traj_form_joystick_abs(joystick):
    """

    :param joystick:
    :return:
    """
    centre, direction = joystick[0], joystick[1]
    traj = []
    r = []
    pos = direct_abs(get_verins_12(), get_O(), get_omega())
    # Computing radii
    for leg in range(4):
        radius = distance(pos[leg*3+0:leg*3+3], (centre[0], centre[1], 0))
        r.append(radius)
    r_max = max(r)
    n = int(2 * np.pi * r_max / TRAJ_ACCUR)
    for i in range(n-250):
        L = []
        for leg in range(4):
            alpha = np.arccos(abs((centre[1] - pos[3 * leg + 1])) / r[leg])
            signe_cx = (centre[0] - pos[leg * 3 + 0]) / abs(centre[0] - pos[leg * 3 + 0])
            signe_cy = (centre[1] - pos[leg * 3 + 1]) / abs(centre[1] - pos[leg * 3 + 1])
            if signe_cx < 0 and signe_cy < 0:
                alpha = + np.pi / 2 - alpha
            if signe_cx > 0 > signe_cy:
                alpha = + np.pi / 2 + alpha
            if signe_cx < 0 < signe_cy:
                alpha = - np.pi / 2 + alpha
            if signe_cx > 0 and signe_cy > 0:
                alpha = - np.pi / 2 - alpha
            L = np.append(L, (r[leg] * np.cos((2 * i * np.pi) / n + alpha) + centre[0],
                              r[leg] * np.sin((2 * i * np.pi) / n + alpha) + centre[1],
                              pos[3 * leg + 2]))
        traj.append(L)
    return traj


def furthest_accessible_abs(traj, leg_id, max_omega=10, const_omega=True, reg_val=0.01, step_height=0):
    """
    Compute the maximum step size following the trajectory, depending on the accessible zone for the leg.
    traj should be the trajectory of all legs (basically returned by compute_traj_from_joystick)
    :param traj: trajectory to follow
    :param leg_id: ID of the leg
    :param max_omega: maximum angle allowed for the body
    :param const_omega: bool -> is there a constraint on angle
    :param reg_val: regularisation value
    :param step_height: height of one step
    :return: the furthest point accessible from the traj
    """
    temp = np.zeros(12)
    traj_leg = [temp]*len(traj)
    # traj_leg = np.zeros((len(traj), 12))
    for i in range(len(traj)):
        traj_leg[i] = traj[0].copy()
        traj_leg[i][leg_id*3:leg_id*3+3] = traj[i][leg_id*3:leg_id*3+3]
    init = traj[0].copy()
    for j in range(step_height):
        traj_leg.insert(0, init.copy())
        traj_leg[0][leg_id * 3 + 2] += step_height + 1 - j
    for k in range(step_height, len(traj_leg)):
        traj_leg[k][leg_id * 3 + 2] += step_height
    R = reg_val * np.eye(18)
    V = get_verins_12()
    O = get_O()
    omega = get_omega()
    LV = [V]
    LO = [O]
    LOmega = [omega]
    for i in range(1, len(traj_leg)):
        # Calcul de dX
        X0 = direct_abs(V, O, omega)
        if distance(X0[leg_id * 3:leg_id * 3 + 3], traj_leg[i-1][leg_id * 3:leg_id * 3 + 3]) > 20:  # exit condition
            return i-step_height, (LV, LO, LOmega)
        set_X(X0)
        dX = traj_leg[i] - X0
        # Contraintes
        lb = np.full(18, - np.inf)
        ub = np.full(18, np.inf)
        for j in range(12):
            lb[j] = 450.0 - V[j]
            ub[j] = 650.0 - V[j]
        for j in range(3):
            if const_omega:
                lb[15 + j] = - max_omega * np.pi / 180 - omega[j]
                ub[15 + j] = max_omega * np.pi / 180 - omega[j]
            else:
                lb[15 + j] = - np.pi
                ub[15 + j] = np.pi
        # solve
        M = jacob_dX_to_dV_dO_dOmega(V, omega, direct_rel_12(V))
        P = M.T @ M + R
        q = - M.T @ dX
        sol = solve_qp(P, q, lb=lb, ub=ub)
        # Mise à jour des valeurs réelles
        V = V + sol[0:12]
        O = O + sol[12:15]
        omega = omega + sol[15:18]
        for v in V: assert 449.9 < v < 650.1, 'Elongation de vérin invalide'
        LV.append(V)
        LO.append(O)
        LOmega.append(omega)
    return len(traj_leg)-1, (LV, LO, LOmega)


def furthest_accessible_all_legs_abs(traj, step_height=0):
    max_step, data = [], []
    for leg in range(4):
        step, data_leg = furthest_accessible_abs(traj, leg, step_height=step_height, max_omega=45)
        max_step.append(step)
        data.append((data_leg))
    return max_step, data


################################ WALK -> SYL ##################################


def get_joystick():
    return (1, 0), 0


def move_along(traj, step, leg_id):
    v1, v2, v3 = get_verins_3(leg_id)
    lpl = ROBOT['legs'][leg_id]['lengths']
    verrins = []
    for i in range(1, step):
        # get data
        x0, y0, z0 = direct_rel_3(v1, v2, v3, leg_id)
        xt, yt, zt = traj[i][leg_id * 3 + 0], traj[i][leg_id * 3 + 1], traj[i][leg_id * 3 + 2]
        dX = np.array([xt - x0, yt - y0, zt - z0])
        # solve
        J = np.dot(gen_jacob_3(v1 / 1000, v2 / 1000, v3 / 1000, np.arccos(v3_to_cos_angle(v3, lpl)), lpl),
                   ROBOT['legs'][leg_id]['matrix'])
        P = np.dot(inv(J).T, inv(J))
        q = - np.dot(inv(J).T, dX)
        lb = np.array([450.0 - v1, 450.0 - v2, 450.0 - v3])
        ub = np.array([650.0 - v1, 650.0 - v2, 650.0 - v3])
        dV = solve_qp(P, q, lb=lb, ub=ub)
        v1, v2, v3 = v1 + dV[0], v2 + dV[1], v3 + dV[2]
        verrins.append((v1, v2, v3))
    return verrins


def waaalkkk():
    """not for now yet"""
    d, r = get_joystick()
    traj = compute_traj_form_joystick_rel(cmd_joystick(d, r))
    max_step = []
    for leg in range(4):
        max_step.append(furthest_accessible_rel(traj, leg))
    step = max(max_step)

    # Move front left leg
    move_along(traj, step, FL)


############################################################################

if __name__ == "__main__":
    import doctest
    doctest.testmod()

######################## DEPRECATED #############################

# def cmd_joystick(D, R):
#     """
#     Traite les données reçues du controleur joystick pour les convertir en une translation dirigée et une rotation.
#     Pour le moment la direction ne peut être que devant, derrière, avancer, reculer.
#
#     :param D: couple de coord du joystick de déplacement
#     :param R: facteur de rotation du joystick de rotation
#     :return: position du centre de rotation, vecteur de direction de la marche
#     # >>> cmd_joystick((1, 0), 0)
#     # ((0, 100000), (1, 0))
#     # >>> cmd_joystick((0, 1), 0)
#     # ((-100000, 0), (0, 1))
#     # >>> cmd_joystick((0.51, 0.49), 0)
#     # ((0, 100000), (1, 0))
#     # >>> cmd_joystick((0.49, 0.51), 1)
#     # ((-2500, 0.0), (0, 1))
#     # >>> cmd_joystick((0.49, -0.51), 1)
#     # ((2500, 0.0), (0, -1))
#
#     """
#     if D[0] == 0 and D[1] == 0:
#         if abs(R) < 0.1:
#             # l'immobilité
#             return "erreur, pas de mouvement en entrée"
#         # Rotation sur lui-même
#         return (0, 0), (0, 0)
#
#     # Direction NSEW
#     if abs(D[0]) >= abs(D[1]):
#         direction = (D[0] > 0) - (D[0] < 0), 0
#     else:
#         direction = 0, (D[1] > 0) - (D[1] < 0)
#
#     if abs(R) < 0.1:
#         # Marche rectiligne
#         centre = direction[1] * INF, direction[0] * INF
#         return centre, direction
#
#     signeR = (R > 0) - (R < 0)
#     facteur = signeR * (MAX_RADIUS - R * signeR * (MAX_RADIUS - MIN_RADIUS))
#     print(R, signeR, facteur)
#     centre = facteur * direction[1], facteur * direction[0]
#     return centre, direction
#
# def angle_between_linear(coef1, coef2):
#     pointO = 0, 0
#     pointA = 1, coef1
#     pointB = 1, coef2
#     return al_kashi_angle(distance(pointO, pointA), distance(pointO, pointB), distance(pointA, pointB))
#
# def to_ref_traj_next_step(alpha, dx, dy):
#     return np.array([
#         [np.cos(alpha), - np.sin(alpha), dx],
#         [np.sin(alpha), np.cos(alpha), dy],
#         [0, 0, 1]
#     ])
#
# def compute_trajs(traj):
#     all_trajs, trajFL, trajFR, trajRL, trajRR = [], [], [], [], []
#     V = get_verins_12()
#     Pos = direct_rel_12(V)
#     fl = (Pos[0], Pos[1], Pos[2])
#     fr = (Pos[3], Pos[4], Pos[5])
#     last_coef = (fr[1] - fl[1]) / (fr[0] - fl[0])
#     if last_coef > 0.1:
#         print("Legs position aren't initialized")
#     legs_spacing = (abs(fr[1]) + abs(fl[1])) / 2
#
#     for i in range(1, len(traj)):
#         direction = (traj[1][0] - traj[0][0]) / (traj[1][1] - traj[0][1])
#         current_coef = -1 / direction
#         alpha = angle_between_linear(last_coef, current_coef)
#         m_rota = to_ref_traj_next_step(alpha, traj[1][0] - traj[0][0], traj[1][1] - traj[0][1])
#         step_fl = np.dot(m_rota, np.array([-legs_spacing, 0, 1]))
#         step_fr = np.dot(m_rota, np.array([legs_spacing, 0, 1]))
#         trajFL.append(np.array([step_fl[0], step_fl[1], traj[i][i % 4][2]]))
#         trajFR.append(np.array([step_fr[0], step_fr[1], traj[i][i % 4][2]]))
#
#     # trajRL.append(draw_line_3(V[6], V[7], V[8], Pos[0] - Pos[6], Pos[1] - Pos[7], Pos[2] - Pos[8], 20, 2))
#     # trajRR.append(draw_line_3(V[9], V[10], V[11], Pos[3] - Pos[9], Pos[4] - Pos[10], Pos[5] - Pos[11], 20, 2))
#     # for i in range(len(trajFL) - len(trajRL)):
#     #     trajRL.append(trajFL[i])
#     #     trajRR.append(trajFR[i])
#     #
#     # for i in range(len(trajFL)):
#     #     traj.append(trajFL[i])
#     #     traj.append(trajFR[i])
#     #     traj.append(trajRL[i])
#     #     traj.append(trajRR[i])

