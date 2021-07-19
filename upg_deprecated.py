# L'ensemble des distances sont exprimées en mm : segments de patte et élongations des vérins

import kinetic as kin
from upg_kinetic import *

x_A = - kin.LEGS[kin.FL]['lengths']['ao']
y_A = 0
x_B = kin.LEGS[kin.FL]['lengths']['bo']
y_B = 0

AB = kin.LEGS[kin.FL]['lengths']['ao'] + kin.LEGS[kin.FL]['lengths']['bo']
AC = np.sqrt((AB + kin.LEGS[kin.FL]['lengths']['bcx']) ** 2 + kin.LEGS[kin.FL]['lengths']['bcy'] ** 2)
BC = np.sqrt(kin.LEGS[kin.FL]['lengths']['bcx'] ** 2 + kin.LEGS[kin.FL]['lengths']['bcy'] ** 2)
AE = kin.LEGS[kin.FL]['lengths']['ae']
AD = AE - kin.LEGS[kin.FL]['lengths']['de']
BF = kin.LEGS[kin.FL]['lengths']['bf']
FH = kin.LEGS[kin.FL]['lengths']['fh']
BH = BF - FH
FG = kin.LEGS[kin.FL]['lengths']['fg']
EF = kin.LEGS[kin.FL]['lengths']['ef']
EG = EF + FG
GI = kin.LEGS[kin.FL]['lengths']['gi']
GJ = kin.LEGS[kin.FL]['lengths']['gj']

KO = kin.LEGS[kin.FL]['lengths']['yaw_c']
LM = kin.LEGS[kin.FL]['lengths']['yaw_b']
MO = kin.LEGS[kin.FL]['lengths']['yaw_a']
LO = np.sqrt(LM ** 2 + MO ** 2)

ori = np.array([[1, 1, -1, -1], [1, -1, -1, 1]])  # [[oritentation selon x][orientation selon y]]

V = 460, 565, 500
alpha = al_kashi_angle(AB, AC, BC)
beta = np.arccos(MO / LO)


def dist_2_legs(l1, l2, V):
    """
    Calcule la distance entre 2 pattes pour une disposition des vérins donnée

    :param l1: patte 1
    :param l2: patte 2
    :param V: élongation des vérins
    :return: distance entre les 2 pattes
    """
    pos = direct_rel_12(V)
    return distance(pos[l1 * 3:l1 * 3 + 3], pos[l2 * 3:l2 * 3 + 3])


def get_last_leg(leg_id, V):
    """
    Retourne les coordonnées absolues de la 4ème patte lorsque les 3 autres sont au sol

    :param leg_id: ID de la patte levée
    :param V: élongation des vérins
    :return: coordonnées de la patte (absolue)

    >>> init()
    >>> set_og(0, 0)
    >>> V = get_verins_12()
    >>> V[0] += 30
    >>> pos_rel = direct_rel_12(V)
    >>> abs(distance(get_last_leg(0, V), get_leg_pos(1)) - distance(pos_rel[0:3], pos_rel[3:6])) < 0.001
    True
    >>> abs(distance(get_last_leg(0, V), get_leg_pos(2)) - distance(pos_rel[0:3], pos_rel[6:9])) < 0.001
    True
    >>> abs(distance(get_last_leg(0, V), get_leg_pos(3)) - distance(pos_rel[0:3], pos_rel[9:12])) < 0.001
    True
    >>> V[0] = 650
    >>> pos_rel = direct_rel_12(V)
    >>> abs(distance(get_last_leg(0, V), get_leg_pos(1)) - distance(pos_rel[0:3], pos_rel[3:6])) < 0.001
    True
    >>> abs(distance(get_last_leg(0, V), get_leg_pos(2)) - distance(pos_rel[0:3], pos_rel[6:9])) < 0.001
    True
    >>> abs(distance(get_last_leg(0, V), get_leg_pos(3)) - distance(pos_rel[0:3], pos_rel[9:12])) < 0.001
    True
    >>> set_og(1, 0)
    """
    pts = []
    d = []
    for i in range(4):
        if i != leg_id:
            pts.append(get_leg_pos(i))
            d.append(dist_2_legs(i, leg_id, V))

    C1 = pts[0][0] ** 2 + pts[0][1] ** 2 + pts[0][2] ** 2 - d[0] ** 2 - \
         (pts[1][0] ** 2 + pts[1][1] ** 2 + pts[1][2] ** 2 - d[1] ** 2)

    C2 = pts[0][0] ** 2 + pts[0][1] ** 2 + pts[0][2] ** 2 - d[0] ** 2 - \
         (pts[2][0] ** 2 + pts[2][1] ** 2 + pts[2][2] ** 2 - d[2] ** 2)

    x_A = ((pts[2][1] - pts[0][1]) * C1 + (pts[0][1] - pts[1][1]) * C2) / \
          (2 * (pts[1][1] * pts[2][0] + pts[2][1] * pts[0][0] + pts[0][1] * pts[1][0] -
                (pts[0][1] * pts[2][0] + pts[1][1] * pts[0][0] + pts[2][1] * pts[1][0])))

    y_A = (2 * (pts[2][0] - pts[1][0]) * x_A + C2 - C1) / (2 * (pts[1][1] - pts[2][1]))

    z_A2 = d[2] ** 2 - (x_A - pts[2][0]) ** 2 - (y_A - pts[2][1]) ** 2
    if z_A2 < 0:
        z_A = 0.
        print("Patte dans le sol !!! z_A2 = ", z_A2)
        print("d_AD : ", d[2])
        print("d_ID : ", np.sqrt((x_A - pts[2][0]) ** 2 + (y_A - pts[2][1]) ** 2))
        print("x_A : ", x_A)
        print("x_D : ", pts[2][0])
        print("y_A : ", y_A)
        print("y_D : ", pts[2][1])
    else:
        z_A = np.sqrt(z_A2)

    return np.array([x_A, y_A, z_A])


def old_direct_abs(V):
    """
    Calcule la nouvelle position absolue des pattes à partir de la position absolue initiale des pattes pour
    la phase de déplacement en cours (une phase = un déplacement avec le même groupe de pattes au sol)

    :param V: élongation des vérins (vecteur 12)
    :return: coordonnées absolues des pattes (vecteur 12)

    >>> init()
    >>> set_og(0, 0)
    >>> V = get_verins_12()
    >>> X0_rel_ini = direct_rel_12(V)[0:3]
    >>> X0_abs_ini = old_direct_abs(V)[0:3]
    >>> V[0] += 3
    >>> X0_rel = direct_rel_12(V)[0:3]
    >>> X0_abs = old_direct_abs(V)[0:3]
    >>> abs(distance(X0_rel_ini, X0_rel) - distance(X0_abs_ini, X0_abs)) < 0.001
    True
    >>> V[0] += 30
    >>> X0_rel = direct_rel_12(V)[0:3]
    >>> X0_abs = old_direct_abs(V)[0:3]
    >>> abs(distance(X0_rel_ini, X0_rel) - distance(X0_abs_ini, X0_abs)) < 0.001
    True
    >>> V[0] = 650
    >>> X1_rel = direct_rel_12(V)[3:6]
    >>> X1_abs = old_direct_abs(V)[3:6]
    >>> abs(distance(X0_rel, X1_rel) - distance(X0_abs, X1_abs)) < 0.001
    True
    >>> V[1] = 650
    >>> V[2] = 650
    >>> X1_rel = direct_rel_12(V)[3:6]
    >>> X1_abs = old_direct_abs(V)[3:6]
    >>> abs(distance(X0_rel, X1_rel) - distance(X0_abs, X1_abs)) < 0.001
    True
    >>> set_og(1, 0)
    """
    X = get_X()
    for i in range(4):
        if not get_og(i):
            pos_abs = get_last_leg(i, V)
            X[3 * i:3 * i + 3] = pos_abs
    return X


def direct_O(X_abs, V):
    """
    Cinématique directe dans le référentiel absolu

    Ne marche pas

    :param X_abs: x absolu
    :param V: liste des 12 élongations
    :return: position du centre dans le référentiel absolue
    """
    if ROBOT['legs'][0]['og'] and ROBOT['legs'][2]['og']:
        O0 = X_abs[0: 3] - direct_leg(V[0], V[1], V[2], 0)
        O2 = X_abs[6: 9] - ROBOT['legs'][2]['matrix'] @ direct_leg(V[6], V[7], V[8], 2)
        return 0.5 * (O0 + O2)
    O1 = X_abs[3: 6] - ROBOT['legs'][1]['matrix'] @ direct_leg(V[3], V[4], V[5], 1)
    O3 = X_abs[9: 12] - ROBOT['legs'][3]['matrix'] @ direct_leg(V[9], V[10], V[11], 3)
    return 0.5 * (O1 + O3)


def draw_circle_2(v1, v2, r, n, leg_id, solved=False):
    lpl = kin.LEGS[leg_id]['lengths']
    pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
    X0, Z0 = pts['J'][0] * 1000, pts['J'][1] * 1000
    alpha = np.cos(np.pi / 4)
    res = []

    # Calcul des points du cercle
    Lx = []
    Lz = []
    for k in range(n + 1):
        Lx.append(X0 + r * np.cos(2 * k * np.pi / n) - r)
        Lz.append(Z0 + r * np.sin(2 * k * np.pi / n))

    # Parcours du cercle
    for k in range(1, n + 1):
        # print("POSITIONS ______actual :", X0, Z0,"__________cible :", Lx[k], Lz[k])
        # print("VERINS_________actual :", v1, v2)
        dX = np.array([Lx[k] - X0, Lz[k] - Z0])
        J = gen_MJ(pts, v1 / 1000, v2 / 1000)

        if solved:  # Utilisation du solveur
            P = 2 * J.T @ J
            q = J.T @ (np.array([X0 / 1000, Z0 / 1000]) - np.array([Lx[k] / 1000, Lz[k] / 1000]))
            lb = np.array([(450.0 - v1), (450.0 - v2)]) / 1000
            ub = np.array([(650.0 - v1), (650.0 - v2)]) / 1000
            dV = solve_qp(P, q, lb=lb, ub=ub) * 1000
        else:  # Utilisation de la jacobienne sans solveur
            dV = J @ dX

        v1 += dV[0]
        v2 += dV[1]
        res.append((v1, v2))
        pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
        X0, Z0 = pts['J'][0] * 1000, pts['J'][1] * 1000
    return res


def draw_circle_3(v1, v2, v3, r, n, leg_id, solved=False):
    lpl = kin.LEGS[leg_id]['lengths']
    pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
    X0, Z0 = pts['J'][0] * 1000, pts['J'][1] * 1000
    x0, y0, z0 = d2_to_d3(X0, Z0, v3_to_cos_angle(v3))
    res = []

    # Calcul des points du cercle
    Lx = []
    Ly = []
    Lz = []
    for k in range(n + 1):
        Lx.append(x0 + r * np.cos(2 * k * np.pi / n) - r)
        Ly.append(y0 + r * np.sin(2 * k * np.pi / n))
        Lz.append(z0)

    # Parcours du cercle
    for k in range(1, n + 1):
        # print("POSITIONS ______actual :",x0, y0, z0,"__________cible :", Lx[k], Ly[k], Lz[k])
        # print("VERINS_________actual :", v1, v2, v3)
        dX = np.array([Lx[k] - x0, Ly[k] - y0, Lz[k] - z0])
        alpha0 = np.arccos(v3_to_cos_angle(v3))
        J = gen_jacob_3(v1 / 1000, v2 / 1000, v3 / 1000, alpha0, kin.FL)

        if solved:  # Utilisation du solveur
            dV = np.array([0, 0, 0])
            # P = 2 * J.T @ J
            # q = J.T @ (np.array([X0/1000, Z0/1000]) - np.array([X/1000, Z/1000]))
            # lb = np.array([(450.0 - v1), (450.0 - v2)])/1000
            # ub = np.array([(650.0 - v1), (650.0 - v2)])/1000
            # dV = solve_qp(P, q, lb=lb, ub=ub) * 1000
        else:  # Utilisation de la jacobienne sans solveur
            dV = J @ dX

        v1 += dV[0]
        v2 += dV[1]
        v3 += dV[2]
        res.append((v1, v2, v3))
        pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
        X0, Z0 = pts['J'][0] * 1000, pts['J'][1] * 1000
        x0, y0, z0 = d2_to_d3(X0, Z0, v3_to_cos_angle(v3))
    return res


def draw_line_3(v1, v2, v3, dx, dy, dz, n, leg_id):
    """
    Retourne trajectoire rectiligne en liste de point pour un certain décalage et une certaine patte

    :param v1: elongation de v1
    :param v2: elongation de v2
    :param v3: elongation de v3
    :param dx: décalage en x
    :param dy: décalage en y
    :param dz: décalage en z
    :param n: précision de la discrétisation (nombre de points intermédiaires)
    :param leg_id: id de la jambe
    :return: trajectoire en liste coordonnées intermédiaire
    """
    x0, y0, z0 = direct_leg(v1, v2, v3, FL)
    traj = []
    for i in range(n):
        traj.append((x0 + i * dx / n, y0 + i * dy / n, z0 + i * dz / n))
    return traj


def draw_line_12(V, dx, dy, dz, n):
    """
    Retourne la discrétisation d'une trajectoire en ligne droite d'un certain décalage pour les 4 pattes à la fois

    :param V: liste des 12 vérins
    :param dx: décalage en x
    :param dy: décalage en y
    :param dz: décalage en z
    :param n: nombre d'étapes dans la traj
    :return: la trajextoire
    """
    traj_FL = draw_line_3(V[0], V[1], V[2], 15, 20, 10, 10, 0)
    traj_FR = draw_line_3(V[3], V[4], V[5], 15, 20, 10, 10, 1)
    traj_RL = draw_line_3(V[6], V[7], V[8], 15, 20, 10, 10, 2)
    traj_RR = draw_line_3(V[9], V[10], V[11], 15, 20, 10, 10, 3)
    traj = []
    for i in range(n):
        t = traj_FL[i]
        t = np.append(t, ROBOT['legs'][1]['matrix'] @ traj_FR[i])
        t = np.append(t, ROBOT['legs'][2]['matrix'] @ traj_RL[i])
        t = np.append(t, ROBOT['legs'][3]['matrix'] @ traj_RR[i])
        traj.append(t)
    return traj


def move_4_legs(traj, V, upgrade=False):
    """
    Retourne le tableau des élongations successives des 12 vérins (3 par 3) permettant aux
    extrémités des 4 pattes de suivre les trajectoires qui leur ont été attribuées par traj
    traj : [traj_0, traj_1, traj_2, traj_3], liste des 4 trajectoires suivies par les 4 pattes avec traj_i = [[x0, y0, z0], [x1, y1, z1], ...]
    V : [V_0, V_1, V_2, V_3], liste des élongations initiales des 3 vérins dans chacune des 4 pattes avec V_i = [v1, v2, v3]
    Toutes les longueurs sont en mm, les coordonnées des trajectoires sont exprimées dans le repère du robot
    """
    # Calcul des élongations de chacunes des pattes
    Ver = []
    print(traj)
    for i in range(4):
        Ver.append(move_leg(traj[i], V[i][0], V[i][1], V[i][2], i, upgrade=upgrade, solved=True))
    # Mise sous le format attendu
    R = []
    for k in range(len(Ver[0])):
        r = [Ver[0][k],
             Ver[1][k],
             Ver[2][k],
             Ver[3][k]]
        R.append(r)
    return R


def solve_indirect_cyl(x, y, z, x0, y0, z0, v1, v2, v3, lpl, pts):
    """
    Retourne les élongations (v1, v2, v3) permettant de placer le bout de la patte en (x, y, z)
    en minimisant les erreurs d'élongation sur v1 et v2 (utilise Jacob_2)
    Toutes les longueurs sont en mm (entrées comme sorties)
    """
    X, Z, calpha = d3_to_d2(x, y, z)
    Xt = np.array([X, Z]) / 1000
    new_v3 = cos_angle_to_v3(calpha, lpl)

    X, Z, calpha = d3_to_d2(x0, y0, z0)
    X0 = np.array([X, Z]) / 1000

    J = gen_MJ(pts, v1 / 1000, v2 / 1000)

    P = 2 * J.T @ J
    q = J.T @ (X0 - Xt)
    lb = np.array([(450.0 - v1) / 1000, (450.0 - v2) / 1000])
    ub = np.array([(650.0 - v1) / 1000, (650.0 - v2) / 1000])
    dV = solve_qp(P, q, lb=lb, ub=ub)

    return v1 + dV[0] * 1000, v2 + dV[1] * 1000, new_v3


def solve_indirect_cart(x, y, z, x0, y0, z0, v1, v2, v3, lpl, pts):
    """
    Retourne les élongations (v1, v2, v3) permettant de placer le bout de la patte en (x, y, z)
    en minimisant les erreurs d'élongation sur v1, v2 et v3 (utilise Jacob_3)
    Toutes les longueurs sont en mm (entrées comme sorties)
    """
    Xt = np.array([x, y, z]) / 1000
    X0 = np.array([x0, y0, z0]) / 1000

    J = gen_jacob_3(v1 / 1000, v2 / 1000, v3 / 1000, np.arctan(x0 / y0), lpl)

    P = J.T @ J
    q = Xt.T @ J
    lb = np.array([450.0 - v1, 450.0 - v2, 450.0 - v3]) / 1000
    ub = np.array([650.0 - v1, 650.0 - v2, 650.0 - v3]) / 1000
    dV = solve_qp(P, q, lb=lb, ub=ub)

    return v1 + dV[0] * 1000, v2 + dV[1] * 1000, v3 + dV[2] * 1000


def normalized_move_xyz(x, y, z, v1, v2, v3, dstep, p, eps, leg_id):
    """
    Retourne la liste des élongations successives des vérins permettant de placer le bout de la patte en (x,y,z)
    Effectue une normalisation du déplacement
    """
    L = []
    c = 0

    pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, kin.LEGS[kin.FL]['lengths'])
    X, Z = pts['J'][0] * 1000, pts['J'][1] * 1000
    x0, y0, z0 = d2_to_d3(X, Z, v3_to_cos_angle(v3))
    dist = distance(x - x0, y - y0, z - z0)
    while dist > eps and c < 300:
        c += 1
        U = np.array([(x - x0), (y - y0), (z - z0)])
        U = dstep / 100 * U  # / np.linalg.norm(U)**p
        dx, dy, dz = U[0], U[1], U[2]

        v1, v2, v3 = solve_indirect_cart(x0 + dx, y0 + dy, z0 + dz, x0, y0, z0, v1, v2, v3, leg_id, pts)
        L.append((v1, v2, v3))

        pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, kin.LEGS[kin.FL]['lengths'])
        X, Z = pts['J'][0] * 1000, pts['J'][1] * 1000
        x1, y1, z1 = x0 + dx, y0 + dy, z0 + dz
        x0, y0, z0 = d2_to_d3(X, Z, v3_to_cos_angle(v3))

        dist = distance(x - x0, y - y0, z - z0)
    return L


def distance_euclidienne(xi, yi, xj, yj):
    return np.sqrt((xj - xi) ** 2 + (yj - yi) ** 2)


def direct_v1(v1, v2):
    theta1 = alpha + al_kashi_angle(AD, AC, v1)
    x_E = x_A + AE * np.cos(theta1)
    y_E = y_A + AE * np.sin(theta1)
    EB = distance_euclidienne(x_E, y_E, x_B, y_B)

    theta2 = al_kashi_angle(EF, EB, BF)
    theta3 = al_kashi_angle(AE, EB, AB)

    beta = theta2 + theta3 - (np.pi - theta1)
    x_F = x_E + EF * np.cos(beta)
    y_F = y_E + EF * np.sin(beta)
    x_G = x_E + EG * np.cos(beta)
    y_G = y_E + EG * np.sin(beta)

    x_H = (FH * x_B + BH * x_F) / (FH + BH)
    y_H = (FH * y_B + BH * y_F) / (FH + BH)

    GH = distance_euclidienne(x_G, x_H, y_G, y_H)

    theta4 = al_kashi_angle(FG, GH, FH)
    theta5 = al_kashi_angle(GI, GH, v2)
    theta6 = np.pi - (theta4 + theta5 + beta)

    return x_G + GJ * np.cos(theta6), y_G - GJ * np.sin(theta6)


def deltas(theta1, theta2, dstep):
    '''
    Fonction auxiliaire de move_xyz
    Retourne les valeurs de dx, dy et dz en fonction des angles theta1 et theta2 et de la distance d'un pas
    '''
    deltaX = dstep * np.cos(theta2)
    dz = dstep * np.sin(theta2)
    dx = deltaX * np.cos(theta1)
    dy = deltaX * np.sin(theta1)
    return dx, dy, dz


def angle_to_v3(angle):
    '''
    Fonction auxiliaire de move_xyz
    Retourne l'élongation de v3 en fonction de l'angle de la patte au chassis

    >>> angle_to_v3(np.pi/4) - al_kashi_longueur(KO, LO, np.pi/4 - beta)
    0.0
    >>> v3_to_angle(angle_to_v3(np.pi/4)) - np.pi/4 < 0.0000001
    True
    '''
    return np.sqrt(LO ** 2 + KO ** 2 - 2 * LO * KO * np.cos(angle - beta))


def v3_to_angle(v3):
    '''
    Fonction auxiliaire de move_xyz
    Retourne l'angle de la patte au chassis en fonction de l'élongation de v3

    >>> v3_to_angle(500) - (al_kashi_angle(LO, KO, 500) + beta)
    0.0
    >>> angle_to_v3(v3_to_angle(500)) - 500 < 0.0000001
    True
    '''
    return np.arccos((KO ** 2 + LO ** 2 - v3 ** 2) / (2 * KO * LO)) + beta
