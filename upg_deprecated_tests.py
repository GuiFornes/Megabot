from upg_tests import *
import time

def test_jacob_2(v1, v2, dstep_x, dstep_y):
    """
    Retourne l'erreur relative en x et en y de l'application de la Jacobienne du modèle indirect plan
    Prend en argument l'élongation des vérins et la distance de déplacement selon x et y en m
    """
    pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, ROBOT['legs'][FR]['lengths'])
    x, y = pts['J'][0], pts['J'][1]
    Jacob = gen_MJ(pts, v1 / 1000, v2 / 1000)

    Dpos = np.array([dstep_x / 1000, dstep_y / 1000])
    DV = Jacob @ Dpos

    pts = kin.get_leg_points_V1_V2(v1 / 1000 + DV[0], v2 / 1000 + DV[1], ROBOT['legs'][FR]['lengths'])
    # print("position initiale : ", x * 1000, y * 1000)

    new_x, new_y = pts['J']
    # print("nouvelle position : ", new_x * 1000, new_y * 1000)
    print("on s'est déplacé de ", (new_x - x) * 1000, " mm en x et de ", (new_y - y) * 1000, " mm en y")

    err_pos = np.array([new_x - (x + dstep_x / 1000), new_y - (y + dstep_y / 1000)])
    # print("erreur de position en x : ", err_pos[0] * 1000, " mm")
    # print("erreur de position en y : ", err_pos[1] * 1000, " mm")
    err_rel = np.array([None, None])
    if dstep_x != 0: err_rel[0] = err_pos[0] * 1000 / dstep_x
    if dstep_y != 0: err_rel[1] = err_pos[1] * 1000 / dstep_y
    if dstep_x != 0: print("erreur relative en x : ", err_rel[0])
    if dstep_y != 0: print("erreur relative en y : ", err_rel[1])
    print("\n")

    return err_rel


def test_normalized_move_xyz(x0, y0, z0, x, y, z, dstep, p, eps, leg_id):
    """
    Trace la trajectoire de la patte du robot entre 2 points de l'espace en suivant normalized_move_xyz
    """
    # Théorique
    Xt = [x0, x]
    Yt = [y0, y]
    Zt = [z0, z]

    # Positions des vérins
    L = normalized_move_xyz(x, y, z, V[0], V[1], V[2], dstep, p, eps, leg_id)
    V1 = [v[0] for v in L]
    V2 = [v[1] for v in L]
    V3 = [v[2] for v in L]
    T = np.linspace(0, 1, len(L))

    # Positions du bout de la patte
    Pos = list(map(direct_leg, [v[0] for v in L], [v[1] for v in L], [v[2] for v in L]))
    Xp = [p[0] for p in Pos]
    Yp = [p[1] for p in Pos]
    Zp = [p[2] for p in Pos]

    # Tracé de la position du bout de la patte
    fig = plt.figure()
    ax = fig.gca(projection='3d')  # Affichage en 3D
    ax.plot(Xt, Yt, Zt, label='Théorique')  # Tracé de la courbe 3D
    ax.scatter(Xp, Yp, Zp, label='Positions',
               marker='d')  # Tracé des points 3D
    plt.title("Trajectoire patte FL avec (dstep, p, eps) = (" + str(dstep) + ", " + str(p) + ", " + str(eps) + ")")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xbound(400, 800)
    ax.set_ybound(400, 800)
    ax.set_zbound(-300, -700)
    plt.tight_layout()
    plt.show()

    # Tracé des positions des vérins
    plt.plot(T, V1, label='V1')
    plt.plot(T, V2, label='V2')
    plt.plot(T, V3, label='V3')
    plt.title("Elongations des vérins dans le mouvement")
    plt.show()


def test_A(dX, dZ, dalpha, v1, v2, v3):
    lpl = ROBOT['legs'][FL]['lengths']
    # alpha = np.cos(v3_to_cos_angle(v3))
    # print (alpha)
    # pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
    # A = mat_A(pts, v1, v2, v3, alpha)
    #
    # dPos = np.array([dX/1000, dZ/1000, dalpha])
    # dV = A @ dPos
    #
    # new_v1, new_v2 = v1 + dV[0]*1000, v2 + dV[1]*1000
    # new_pts = kin.get_leg_points_V1_V2(new_v1 / 1000, new_v2 / 1000, lpl)
    #
    # new_alpha = np.arccos(v3_to_cos_angle(v3 + dV[2]*1000))
    # print (dV)
    # print(new_alpha)
    # err_rel = (abs(new_pts['J'][0] - pts['J'][0])*1000 - dX)/dX, (abs(new_pts['J'][1] - pts['J'][1])*1000 - dZ)/dZ, \
    #           (abs(new_alpha - alpha) - dalpha)/dalpha

    alpha = np.arccos(v3_to_cos_angle(v3, lpl))
    pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
    dPos = np.array([dX / 1000, dZ / 1000, dalpha])

    # en utilisant mat A :
    A = mat_A(pts, v1, v2, v3, alpha)
    dVa = A @ dPos
    v1a = v1 + dVa[0] * 1000
    v2a = v2 + dVa[1] * 1000
    v3a = v3 + dVa[2] * 1000
    print(direct_leg(v1a, v2a, v3a))

    # en séparant jacob 2D et v3
    J = gen_MJ(pts, v1, v2)
    dPos2d = np.array([dPos[0], dPos[1]])
    dVj = J @ dPos2d
    v1j = v1 + dVj[0] * 1000
    v2j = v2 + dVj[1] * 1000
    v3j = cos_angle_to_v3(np.cos(alpha + dalpha))
    print(direct_leg(v1j, v2j, v3j))


def test_jacob_2_direct():
    """
    Comparaison des resultats en J du modèle direct plan utilisant une jacobienne et de celui utilisant la géométrie
    """
    leg_id = kin.FL
    V = np.array([495, 585, 515])
    x0, y0, z0 = direct_leg(V[0], V[1], V[2])
    print("position de départ :", x0, y0, z0)
    deltaV = np.array([5, -5, 5])
    V += deltaV
    deltaV2 = np.array([deltaV[0], deltaV[1]])
    Jacob = gen_jacob_direct(kin.get_leg_points_V1_V2(V[0] / 1000, V[1] / 1000, ROBOT['legs'][leg_id]['lengths']),
                             V[0] / 1000, V[1] / 1000)
    deltaX = Jacob @ deltaV2
    calpha = v3_to_cos_angle(V[2])
    deltax = np.array([
        deltaX[0] * np.sqrt(1 - calpha ** 2),
        deltaX[0] * calpha,
        deltaX[1]
    ])
    print("position selon jacob_direct :", x0 + deltax[0], y0 + deltax[1], z0 + deltax[2])
    print("position selon direct Julien :", direct_leg(V[0], V[1], V[2]))
    print("")


def test_precision_jacobienne_en_direct(v1, v2, dv1, dv2, leg_id):
    """
    Test des précisions pas à pas pour chaque point de la patte en suivant une méthode géométrique exacte et une methode
    approchée utilisant la jacobienne
    @param v1: position initiale du vérin 1
    @param v2: position initiale du vérin 2
    @param dv1: déplacement du vérin 1
    @param dv2: déplacement du vérin 2
    @param leg_id: id de la patte
    @return:
    """
    pts = kin.get_leg_points_V1_V2(v1, v2, ROBOT['legs'][leg_id]['lengths'])
    x_E, y_E = pts['E']
    x_F, y_F = pts['F']
    x_G, y_G = pts['G']
    x_H, y_H = pts['H']
    x_I, y_I = pts['I']
    print('initial D :', pts['D'])
    pts_end = kin.get_leg_points_V1_V2((v1 + dv1), (v2 + dv2), ROBOT['legs'][leg_id]['lengths'])
    A = distance_3_points(pts['D'], pts['A'], pts['C'])
    B = np.array([0, 2 * v1])
    M_D = inv(A) @ B
    deltaD = M_D * dv1
    D = pts['D'] + deltaD
    print("real D :", pts_end['D'])
    print("jacob D :", D, "\n")

    M_E = (ROBOT['legs'][FL]['lengths']['ae'] / (
            ROBOT['legs'][FL]['lengths']['ae'] - ROBOT['legs'][FL]['lengths']['de'])) * M_D
    deltaE = M_E * dv1
    E = pts['E'] + deltaE
    print('initial E :', pts['E'])
    print("real E :", pts_end['E'])
    print("jacob E :", E, "\n")

    A = distance_3_points(pts['F'], pts['E'], pts['B'])
    B = np.array([
        [2 * (x_F - x_E), 2 * (y_F - y_E)],
        [0, 0]])
    M_F = (inv(A) @ B) @ M_E
    deltaF = M_F * dv1
    F = pts['F'] + deltaF
    print('initial F:', pts['F'])
    print("real F:", pts_end['F'])
    print("jacob F", F, "\n")

    M_G = ((ROBOT['legs'][FL]['lengths']['ef'] + ROBOT['legs'][FL]['lengths']['fg']) / ROBOT['legs'][FL]['lengths'][
        'ef']) * M_F
    deltaG = M_G * dv1
    G = pts['G'] + deltaG
    print('initial G:', pts['G'])
    print("real G:", pts_end['G'])
    print("jacob G", G, "\n")

    M_H = (BH / BF) * M_F
    deltaH = M_H * dv1
    H = pts['H'] + deltaH
    print('initial H:', pts['H'])
    print("real H:", pts_end['H'])
    print("jacob H", H, "\n")

    A = distance_3_points(pts['I'], pts['G'], pts['H'])
    B = np.array([
        [2 * (x_I - x_G), 2 * (y_I - y_G)],
        [0, 0]])
    C = np.array([
        [0, 0],
        [2 * (x_I - x_H), 2 * (y_I - y_H)]])
    D = np.array([0, 2 * v2])
    V1 = inv(A) @ (B @ M_G + C @ M_H)
    V2 = inv(A) @ D
    M_I = np.array([
        [V1[0], V2[0]],
        [V1[1], V2[1]]])
    deltaI = M_I @ np.array([dv1, dv2])
    I = pts['I'] + deltaI
    print('initial I:', pts['I'])
    print("real I:", pts_end['I'])
    print("jacob I", I, "\n")

    M_J = (ROBOT['legs'][FL]['lengths']['gj'] / ROBOT['legs'][FL]['lengths']['gi']) * M_I
    deltaJ = M_J @ np.array([dv1, dv2])
    J = pts['J'] + deltaJ
    print('initial J:', pts['J'])
    print("real J:", pts_end['J'])
    print("jacob J", J, "\n")


def test_comparaison_minimize_vs_jacob_indirect(v1, v2, dx, dy):
    """
    test de comparaison des méthodes de cinématique indirecte
    """
    leg_id = FL
    lpl = ROBOT['legs'][leg_id]['lengths']
    pts = kin.get_leg_points_V1_V2(v1, v2, lpl)
    x0, y0 = pts['J']
    dX = np.array([dx, dy])
    x, y = x0 + dx, y0 + dy
    print("initialement v1, v2 =", v1, v2, "_________position en x, y =", x0, y0, "_______cible x, y =", x, y)

    r = kin.inverse_kinetic(x, y, lpl)
    v1_minimize, v2_minimize = r[0], r[1]
    print("selon minimize, V =", r[0], r[1], "________________________________position atteinte en x, y =",
          kin.get_leg_points_V1_V2(r[0], r[1], lpl)['J'])

    J = gen_MJ(pts, v1, v2)
    dV = J @ dX
    print("selon la methode jacobienne, V =", v1 + dV[0], v2 + dV[1], "___________________position atteinte en x, y =",
          kin.get_leg_points_V1_V2(v1 + dV[0], v2 + dV[1], lpl)['J'])

    P = 2 * J.T @ J
    q = J.T @ (np.array([x0, y0]) - np.array([x, y]))
    lb = np.array([(0.450 - v1), (0.450 - v2)])
    ub = np.array([(0.650 - v1), (0.650 - v2)])
    dV = solve_qp(P, q, lb=lb, ub=ub)
    print("en appliquant un solveur sur la jacobienne, V =", v1 + dV[0], v2 + dV[1], "_____position atteinte en x, y =",
          kin.get_leg_points_V1_V2(v1 + dV[0], v2 + dV[1], lpl)['J'])
    print("\n")


def gen_jacob_direct(pts, v1, v2):
    """
    Retourne la Jacobienne correspondant au modèle cinématique indirect dans le plan de la patte
    Prend en argument la position des points de la patte et l'élongation des verrins en m

    >>> gen_MJ(kin.get_leg_points_V1_V2(0.495, 0.585, ROBOT['legs'][FL]['lengths']), ROBOT['legs'][FL]['lengths'], 0.495, 0.585) @ np.array([0, 0])
    array([0., 0.])
    """
    x_E, y_E = pts['E']
    x_F, y_F = pts['F']
    x_G, y_G = pts['G']
    x_H, y_H = pts['H']
    x_I, y_I = pts['I']

    A = distance_3_points(pts['D'], pts['A'], pts['C'])
    B = np.array([0, 2 * v1])
    M_D = inv(A) @ B

    M_E = (ROBOT['legs'][FL]['lengths']['ae'] / (
            ROBOT['legs'][FL]['lengths']['ae'] - ROBOT['legs'][FL]['lengths']['de'])) * M_D

    A = distance_3_points(pts['F'], pts['E'], pts['B'])
    B = np.array([
        [2 * (x_F - x_E), 2 * (y_F - y_E)],
        [0, 0]])
    M_F = (inv(A) @ B) @ M_E

    M_G = ((ROBOT['legs'][FL]['lengths']['ef'] + ROBOT['legs'][FL]['lengths']['fg']) / ROBOT['legs'][FL]['lengths'][
        'ef']) * M_F

    M_H = (BH / BF) * M_F

    A = distance_3_points(pts['I'], pts['G'], pts['H'])
    B = np.array([
        [2 * (x_I - x_G), 2 * (y_I - y_G)],
        [0, 0]])
    C = np.array([
        [0, 0],
        [2 * (x_I - x_H), 2 * (y_I - y_H)]])
    D = np.array([0, 2 * v2])
    V1 = inv(A) @ (B @ M_G + C @ M_H)
    V2 = inv(A) @ D
    M_I = np.array([
        [V1[0], V2[0]],
        [V1[1], V2[1]]])

    return (ROBOT['legs'][FL]['lengths']['gj'] / ROBOT['legs'][FL]['lengths']['gi']) * M_I


def make_a_penalty(v1, v2, d, eps, leg_id):
    lpl = ROBOT['legs'][leg_id]['lengths']
    pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
    X0, Z0 = pts['J'][0] * 1000, pts['J'][1] * 1000
    alpha = np.cos(np.pi / 4)
    res = []

    # drawing the shoot trajectoire
    Lx = np.zeros(d // eps + 1)
    Lz = np.zeros(d // eps + 1)
    for k in range(d // eps + 1):
        X = X0 + k * eps
        Z = Z0
        Lx[k], Lz[k] = X, Z
    print(Lx, Lz)
    # shooting that penalty
    for k in range(1, d // eps + 1):
        X, Z = Lx[k], Lz[k]
        print("POSITION ______actual :", X0, Z0, "__________cible :", X, Z)
        print("VERINS_________actual :", v1, v2)
        dX = np.array([X - X0, Z - Z0])
        J = gen_MJ(pts, v1 / 1000, v2 / 1000)
        P = 2 * J.T @ J
        q = J.T @ (np.array([X0 / 1000, Z0 / 1000]) - np.array([X / 1000, Z / 1000]))
        lb = np.array([(450.0 - v1), (450.0 - v2)]) / 1000
        ub = np.array([(650.0 - v1), (650.0 - v2)]) / 1000
        dV = J @ dX
        # print(dV)
        # dV = solve_qp(P, q, lb=lb, ub=ub)
        # print(dV)
        v1 += dV[0]
        v2 += dV[1]
        res.append((v1, v2))
        pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
        X0, Z0 = pts['J'][0] * 1000, pts['J'][1] * 1000
    return res


def test_circle_2(v1, v2, r, n, leg_id):
    """
    Trace la trajectoire de la patte du robot réalisant des petits cercles
    """
    # Positions des vérins
    lpl = ROBOT['legs'][FL]['lengths']
    L = draw_circle_2(v1, v2, r, n, leg_id)
    V1 = [v[0] for v in L]
    V2 = [v[1] for v in L]
    T = np.linspace(0, 1, len(L))

    # Positions du bout de la patte
    Pos = list(map(lambda v1, v2, lpl: kin.get_leg_points_V1_V2(v1, v2, lpl)['J'], [v[0] / 1000 for v in L],
                   [v[1] / 1000 for v in L], [lpl for v in L]))
    Xp = [p[0] * 1000 for p in Pos]
    Yp = [p[1] * 1000 for p in Pos]

    # Tracé de la position du bout de la patte
    fig = plt.figure()
    # Cercle théorique
    pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
    X0, Z0 = pts['J'][0] * 1000, pts['J'][1] * 1000
    Lx = np.zeros(n + 1)
    Lz = np.zeros(n + 1)
    for k in range(n + 1):
        X = X0 + r * np.cos(2 * k * np.pi / n) - r
        Z = Z0 + r * np.sin(2 * k * np.pi / n)
        Lx[k], Lz[k] = X, Z
    plt.plot(Lx, Lz)

    plt.scatter(Xp, Yp, label='Positions',
                marker='.', s=3, c='red')  # Tracé des points 3D
    plt.title("Trajectoire patte FL rotation circulaire dand le plan XZ")

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')
    plt.show()

    # Tracé des positions des vérins
    plt.plot(T, V1, label='V1')
    plt.plot(T, V2, label='V2')
    plt.title("Elongations des vérins dans le mouvement")
    # plt.plot.set_xlabel('L (en mm)')
    # plt.plot.set_ylabel('T')
    plt.show()


def draw_move_4_legs(traj, V, upgrade=False):
    """
    Trace la trajectoire des extrémités des pattes du robot suivant traj avec move_4_legs
    """
    # Trajectoires
    Xt0 = [p[0] for p in traj[0]]
    Yt0 = [p[1] for p in traj[0]]
    Zt0 = [p[2] for p in traj[0]]
    Xt1 = [p[0] for p in traj[1]]
    Yt1 = [p[1] for p in traj[1]]
    Zt1 = [p[2] for p in traj[1]]
    Xt2 = [p[0] for p in traj[2]]
    Yt2 = [p[1] for p in traj[2]]
    Zt2 = [p[2] for p in traj[2]]
    Xt3 = [p[0] for p in traj[3]]
    Yt3 = [p[1] for p in traj[3]]
    Zt3 = [p[2] for p in traj[3]]

    # Elongations des vérins
    Ver = move_4_legs(traj, V, upgrade=upgrade)
    V01 = [v[0][0] for v in Ver]
    V02 = [v[0][1] for v in Ver]
    V03 = [v[0][2] for v in Ver]
    V11 = [v[1][0] for v in Ver]
    V12 = [v[1][1] for v in Ver]
    V13 = [v[1][2] for v in Ver]
    V21 = [v[2][0] for v in Ver]
    V22 = [v[2][1] for v in Ver]
    V23 = [v[2][2] for v in Ver]
    V31 = [v[3][0] for v in Ver]
    V32 = [v[3][1] for v in Ver]
    V33 = [v[3][2] for v in Ver]
    T = np.linspace(0, 1, len(Ver))

    # Positions du bout de la patte
    Pos0 = list(map(direct_rel_3, [v for v in V01], [v for v in V02], [v for v in V03], [kin.FL for v in Ver]))
    Pos1 = list(map(direct_rel_3, [v for v in V11], [v for v in V12], [v for v in V13], [kin.FR for v in Ver]))
    Pos2 = list(map(direct_rel_3, [v for v in V21], [v for v in V22], [v for v in V23], [kin.RL for v in Ver]))
    Pos3 = list(map(direct_rel_3, [v for v in V31], [v for v in V32], [v for v in V33], [kin.RR for v in Ver]))
    Xp0 = [p[0] for p in Pos0]
    Yp0 = [p[1] for p in Pos0]
    Zp0 = [p[2] for p in Pos0]
    Xp1 = [p[0] for p in Pos1]
    Yp1 = [p[1] for p in Pos1]
    Zp1 = [p[2] for p in Pos1]
    Xp2 = [p[0] for p in Pos2]
    Yp2 = [p[1] for p in Pos2]
    Zp2 = [p[2] for p in Pos2]
    Xp3 = [p[0] for p in Pos3]
    Yp3 = [p[1] for p in Pos3]
    Zp3 = [p[2] for p in Pos3]

    fig = plt.figure()
    ax = fig.gca(projection='3d')  # Affichage en 3D
    ax.plot(Xt0, Yt0, Zt0, label='Théorique', c='coral')  # Tracé de la courbe théorique
    ax.scatter(Xp0, Yp0, Zp0, label='Positions', marker='.', s=3, c='red')  # Tracé des points réels
    ax.plot(Xt1, Yt1, Zt1, label='Théorique', c='cyan')  # Tracé de la courbe théorique
    ax.scatter(Xp1, Yp1, Zp1, label='Positions', marker='.', s=3, c='blue')  # Tracé des points réels
    ax.plot(Xt2, Yt2, Zt2, label='Théorique', c='deeppink')  # Tracé de la courbe théorique
    ax.scatter(Xp2, Yp2, Zp2, label='Positions', marker='.', s=3, c='purple')  # Tracé des points réels
    ax.plot(Xt3, Yt3, Zt3, label='Théorique', c='chartreuse')  # Tracé de la courbe théorique
    ax.scatter(Xp3, Yp3, Zp3, label='Positions', marker='.', s=3, c='green')  # Tracé des points réels
    plt.title("Trajectoire du bout de la patte dans l'espace")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xbound(-2000, 2000)
    ax.set_ybound(-2000, 2000)
    ax.set_zbound(2000, -2000)
    plt.show()

    plt.plot(T, V01, label='V01')
    plt.plot(T, V02, label='V02')
    plt.plot(T, V03, label='V03')
    plt.plot(T, V11, label='V11')
    plt.plot(T, V12, label='V12')
    plt.plot(T, V13, label='V13')
    plt.plot(T, V21, label='V21')
    plt.plot(T, V22, label='V22')
    plt.plot(T, V23, label='V23')
    plt.plot(T, V31, label='V31')
    plt.plot(T, V32, label='V32')
    plt.plot(T, V33, label='V33')
    plt.title("Elongations des vérins dans le mouvement")
    plt.show()


def draw_move_leg(traj, v1, v2, v3, leg_id, upgrade=False, solved=False, anim=False):
    """
        Trace la trajectoire de la patte du robot suivant traj avec move_leg
    """
    lpl = ROBOT['legs'][leg_id]['lengths']
    # Trajectoire
    Xt = [p[0] - 500 for p in traj]
    Yt = [p[1] - 500 for p in traj]
    Zt = [p[2] for p in traj]

    # Elongations des vérins
    Ver = move_leg(traj, v1, v2, v3, leg_id, upgrade=upgrade, solved=solved)
    V1 = [v[0] for v in Ver]
    V2 = [v[1] for v in Ver]
    V3 = [v[2] for v in Ver]
    T = np.linspace(0, 1, len(Ver))

    # Positions du bout de la patte
    Pos = list(map(direct_leg, [v[0] for v in Ver], [v[1] for v in Ver], [v[2] for v in Ver]))
    Xp = [p[0] for p in Pos]
    Yp = [p[1] for p in Pos]
    Zp = [p[2] for p in Pos]
    print(Xp[0])
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')  # Affichage en 3D
    # ax.plot(Xt, Yt, Zt, label='Théorique')  # Tracé de la courbe théorique
    # ax.scatter(Xp, Yp, Zp, label='Positions', marker='.', s=3, c='red')  # Tracé des points réels
    # plt.title("Trajectoire du bout de la patte dans l'espace")
    # ax.set_xlabel('X')
    # ax.set_ylabel('Y')
    # ax.set_zlabel('Z')
    # ax.set_xbound(0, 1000)
    # ax.set_ybound(0, 1000)
    # ax.set_zbound(0, -1000)
    # plt.show()

    # plt.plot(T, V1, label='V1' )
    # plt.plot(T, V2, label='V2' )
    # plt.plot(T, V3, label='V3' )
    # plt.title("Elongations des vérins dans le mouvement")
    # plt.show()

    # ErrX ,ErrY, ErrZ = [], [], []
    # for k in range(len(Ver)):
    #     ErrX.append(Xp[k]-Xt[k])
    #     ErrY.append(Yp[k]-Yt[k])
    #     ErrZ.append(Zp[k]-Zt[k])
    # plt.plot(ErrX, label="x error")
    # plt.plot(ErrY, label="y error")
    # plt.plot(ErrZ, label="z error")
    # plt.legend()
    # plt.ylabel('error in the movement (mm)')
    # plt.xlabel('steps')
    # plt.show()

    Pos2 = list(map(lambda v1, v2: kin.get_leg_points_V1_V2(v1, v2, lpl), [v[0] / 1000 for v in Ver],
                    [v[1] / 1000 for v in Ver]))
    for i in range(len(traj)):
        Pos2[i] = dict(map(lambda kv: (kv[0], d2_to_d3(kv[1][0] * 1000, kv[1][1] * 1000, v3_to_cos_angle(V3[i], lpl))),
                           Pos2[i].items()))
        # Pos2[i] = dict(map(lambda kv: (kv[0], kv[1]*1000), Pos2[i].items()))

    # print(Pos2[0]['A'][0])
    if anim:
        fig = plt.figure()
        for i in range(len(traj)):
            ax = fig.add_subplot(111, projection='3d')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            ax.set_xbound(-500, 500)
            ax.set_ybound(-500, 500)
            ax.set_zbound(-500, 500)
            ax.scatter(xs=0, ys=0, zs=0, s=40)  # point 0
            ax.scatter(xs=Xp[i], ys=Yp[i], zs=Zp[i], s=40, marker='o', c='red')  # point J
            ax.plot(xs=[0, Pos2[i]['A'][0]], ys=[0, Pos2[i]['A'][1]], zs=[0, Pos2[i]['A'][2]], c='black')  # seg OA
            ax.plot(xs=[Pos2[i]['A'][0], Pos2[i]['E'][0]], ys=[Pos2[0]['A'][1], Pos2[i]['E'][1]],
                    zs=[Pos2[i]['A'][2], Pos2[i]['E'][2]], c='black')  # segment AE
            ax.plot(xs=[Pos2[i]['E'][0], Pos2[i]['G'][0]], ys=[Pos2[0]['E'][1], Pos2[i]['G'][1]],
                    zs=[Pos2[i]['E'][2], Pos2[i]['G'][2]], c='black')  # segment EG
            ax.plot(xs=[Pos2[i]['G'][0], Pos2[i]['J'][0]], ys=[Pos2[0]['G'][1], Pos2[i]['J'][1]],
                    zs=[Pos2[i]['G'][2], Pos2[i]['J'][2]], c='black')  # segment GJ
            ax.plot(xs=[Pos2[i]['A'][0], Pos2[i]['B'][0]], ys=[Pos2[0]['A'][1], Pos2[i]['B'][1]],
                    zs=[Pos2[i]['A'][2], Pos2[i]['B'][2]], c='black')  # segment AB
            ax.plot(xs=[Pos2[i]['B'][0], Pos2[i]['F'][0]], ys=[Pos2[0]['B'][1], Pos2[i]['F'][1]],
                    zs=[Pos2[i]['B'][2], Pos2[i]['F'][2]], c='black')  # segment BF
            ax.plot(xs=[Pos2[i]['H'][0], Pos2[i]['I'][0]], ys=[Pos2[0]['H'][1], Pos2[i]['I'][1]],
                    zs=[Pos2[i]['H'][2], Pos2[i]['I'][2]], c='red')  # vérin v2
            ax.plot(xs=[Pos2[i]['B'][0], Pos2[i]['C'][0]], ys=[Pos2[0]['B'][1], Pos2[i]['C'][1]],
                    zs=[Pos2[i]['B'][2], Pos2[i]['C'][2]], c='black')  # segment BC
            ax.plot(xs=[Pos2[i]['C'][0], Pos2[i]['D'][0]], ys=[Pos2[0]['C'][1], Pos2[i]['D'][1]],
                    zs=[Pos2[i]['C'][2], Pos2[i]['D'][2]], c='red')  # vérin v1

            for j in range(i):
                pointJ = ax.scatter(xs=Xp[j], ys=Yp[j], zs=Zp[j], s=20, marker='o', c='grey')
            plt.pause(0.05)
            plt.gcf().clear()


def test_circle_3(v1, v2, v3, r, n, leg_id):
    """
    Trace la trajectoire de la patte du robot réalisant des petits cercles
    """
    # Positions des vérins
    lpl = ROBOT['legs'][kin.FL]['lengths']
    L = draw_circle_3(v1, v2, v3, r, n, leg_id)
    V1 = [v[0] for v in L]
    V2 = [v[1] for v in L]
    V3 = [v[2] for v in L]
    T = np.linspace(0, 1, len(L))

    # Positions du bout de la patte
    Pos = list(map(direct_leg, [v[0] for v in L], [v[1] for v in L], [v[2] for v in L]))
    Xp = [p[0] for p in Pos]
    Yp = [p[1] for p in Pos]
    Zp = [p[2] for p in Pos]

    pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
    X0, Z0 = pts['J'][0] * 1000, pts['J'][1] * 1000
    x0, y0, z0 = d2_to_d3(X0, Z0, v3_to_cos_angle(v3))
    res = []

    # Calcul des points du cercle théorique
    Lx = []
    Ly = []
    Lz = []
    for k in range(n + 1):
        Lx.append(x0 + r * np.cos(2 * k * np.pi / n) - r)
        Ly.append(y0 + r * np.sin(2 * k * np.pi / n))
        Lz.append(z0)

    # Tracé de la position du bout de la patte
    fig = plt.figure()
    ax = fig.gca(projection='3d')  # Affichage en 3D
    ax.plot(Lx, Ly, Lz, label='Théorique')  # Tracé de la courbe 3D
    ax.scatter(Xp, Yp, Zp, label='Positions',
               marker='.', s=3, c='red')  # Tracé des points 3D
    plt.title("Trajectoire cercle dans l'espace")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xbound(200, 700)
    ax.set_ybound(400, 900)
    ax.set_zbound(-200, -700)
    plt.tight_layout()
    plt.show()

    # Tracé des positions des vérins
    plt.plot(T, V1, label='V1')
    plt.plot(T, V2, label='V2')
    plt.plot(T, V3, label='V3')
    plt.title("Elongations des vérins dans le mouvement")
    plt.show()


def test_penalty_move_XZ(v1, v2, d, eps, leg_id):
    """
    Trace la trajectoire de la patte du robot réalisant des petits cercles
    """
    # Positions des vérins
    lpl = ROBOT['legs'][FL]['lengths']
    L = make_a_penalty(v1, v2, d, eps, leg_id)
    V1 = [v[0] for v in L]
    V2 = [v[1] for v in L]
    T = np.linspace(0, 1, len(L))

    # Positions du bout de la patte
    Pos = list(map(lambda v1, v2, lpl: kin.get_leg_points_V1_V2(v1, v2, lpl)['J'], [v[0] / 1000 for v in L],
                   [v[1] / 1000 for v in L], [lpl for v in L]))
    Xp = [p[0] * 1000 for p in Pos]
    Yp = [p[1] * 1000 for p in Pos]

    # Tracé de la position du bout de la patte
    fig = plt.figure()
    # drawing the shoot trajectoire
    pts = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, lpl)
    X0, Z0 = pts['J'][0] * 1000, pts['J'][1] * 1000
    Lx = np.zeros(d // eps + 1)
    Lz = np.zeros(d // eps + 1)
    for k in range(d // eps + 1):
        X = X0 + k * eps
        Z = Z0
        Lx[k], Lz[k] = X, Z

    plt.plot(Lx, Lz)

    plt.scatter(Xp, Yp, label='Positions',
                marker='.', c='red', s=3)  # Tracé des points 3D
    plt.title("Trajectoire rectiligne de la patte FL dans le plan XZ")

    plt.xlabel('X')

    plt.ylabel('Y')
    plt.axis('equal')
    plt.show()

    # Tracé des positions des vérins
    plt.plot(T, V1, label='V1')
    plt.plot(T, V2, label='V2')
    plt.title("Elongations des vérins dans le mouvement")
    # plt.plot.set_xlabel('L (en mm)')
    # plt.plot.set_ylabel('T')
    plt.show()


def test_compute_traj(R, D, abs=False):
    init()
    if abs:
        traj = compute_traj_from_joystick_abs_equals_nb_points(cmd_joystick(R, D))
    else:
        traj = compute_traj_form_joystick_rel(cmd_joystick(R, D))

    V = get_verins_12()

    # Trajectoires
    Xt0 = [p[0] for p in traj]
    Yt0 = [p[1] for p in traj]
    Zt0 = [p[2] for p in traj]
    Xt1 = [p[3] for p in traj]
    Yt1 = [p[4] for p in traj]
    Zt1 = [p[5] for p in traj]
    Xt2 = [p[6] for p in traj]
    Yt2 = [p[7] for p in traj]
    Zt2 = [p[8] for p in traj]
    Xt3 = [p[9] for p in traj]
    Yt3 = [p[10] for p in traj]
    Zt3 = [p[11] for p in traj]

    # Tracé des trajectoires
    fig = plt.figure()
    ax = fig.gca(projection='3d')  # Affichage en 3D
    ax.scatter(Xt0, Yt0, Zt0, label='Théorique', c='coral')  # Tracé de la courbe théorique
    ax.scatter(Xt1, Yt1, Zt1, label='Théorique', c='cyan')  # Tracé de la courbe théorique
    ax.scatter(Xt2, Yt2, Zt2, label='Théorique', c='deeppink')  # Tracé de la courbe théorique
    ax.scatter(Xt3, Yt3, Zt3, label='Théorique', c='chartreuse')  # Tracé de la courbe théorique
    plt.title("Trajectoire du bout de la patte dans l'espace")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    if abs:
        pos = direct_abs(V, get_O(), get_omega())
        ax.plot((-500, -500, 500, 500, -500), (-500, 500, 500, -500, -500), (580, 580, 580, 580, 580))
        ax.plot((500, pos[0]), (500, pos[1]), (580, pos[2]), c='red')
        ax.plot((500, pos[3]), (-500, pos[4]), (580, pos[5]), c='blue')
        ax.plot((-500, pos[6]), (500, pos[7]), (580, pos[8]), c='purple')
        ax.plot((-500, pos[9]), (-500, pos[10]), (580, pos[11]), c='green')
        ax.set_xbound(-2000, 2000)
        ax.set_ybound(-2000, 2000)
        ax.set_zbound(-500, 1500)
    else:
        pos = direct_rel_12(V)
        ax.plot((-500, -500, 500, 500, -500), (-500, 500, 500, -500, -500))
        ax.plot((500, pos[0]), (500, pos[1]), (0, pos[2]), c='red')
        ax.plot((500, pos[3]), (-500, pos[4]), (0, pos[5]), c='blue')
        ax.plot((-500, pos[6]), (500, pos[7]), (0, pos[8]), c='purple')
        ax.plot((-500, pos[9]), (-500, pos[10]), (0, pos[11]), c='green')
        ax.set_xbound(-2000, 2000)
        ax.set_ybound(-2000, 2000)
        ax.set_zbound(-1000, 1000)
    for i in range(4):
        ax.scatter(traj[0][i * 3 + 0], traj[0][i * 3 + 1], traj[0][i * 3 + 2], c='black', s=60)
    ax.set_xbound(-2000, 2000)
    ax.set_ybound(-2000, 2000)
    ax.set_zbound(2000, -2000)
    plt.show()


def test_furthest_pos(D, R):
    init()
    traj = compute_traj_form_joystick_rel(cmd_joystick(D, R))
    max_step = []
    for leg in range(4):
        step = furthest_accessible_rel(traj, leg)
        max_step.append(step)
    print("maximum point of the trajectory reached for each leg : ", max_step)
    min_value = min(max_step)
    leg = max_step.index(min_value)
    step_length = distance((traj[0][leg * 3 + 0], traj[0][leg * 3 + 1], traj[0][leg * 3 + 1]),
                           (traj[min_value][leg * 3 + 0], traj[min_value][leg * 3 + 1], traj[min_value][leg * 3 + 1]))
    print("resultant step length is ", step_length, " for leg n° ", leg)
    comp_rel(traj)


def test_furthest_all_legs(D, R, step_height=0, abs=False):
    init()
    if abs:
        traj = compute_traj_from_joystick_abs_equals_nb_points(cmd_joystick(D, R))
        # print("traj 0 :", traj[0])
        max_step, data = furthest_accessible_all_legs_abs(traj, step_height)
        min_value = min(max_step)
        leg = max_step.index(min_value)
        # print("traj max_step :", traj[max_step])
        print("maximum point of the trajectory reached for each leg : ", max_step)
        step_length = distance((traj[0][leg*3], traj[0][leg*3+1],
                                traj[0][leg*3+2]),
                               (traj[min_value][leg*3], traj[min_value][leg*3+2],
                                traj[min_value][leg*3+2]))
        print("resultant step length is ", step_length, " for leg n° ", leg)
        for leg in range(4):
            draw_abs(data[leg][0], data[leg][1], data[leg][2], traj=traj)
    else:
        traj = compute_traj_form_joystick_rel(cmd_joystick(D, R))
        max_step = furthest_accessible_step_all_legs_rel(traj, step_height)
        print("maximum point of the trajectory reached for each leg : ", max_step)
        min_value = min(max_step)
        leg = max_step.index(min_value)
        step_length = distance((traj[step_height+1][leg * 3 + 0], traj[step_height+1][leg * 3 + 1],
                                traj[step_height+1][leg * 3 + 1]),
                               (traj[min_value+step_height+1][leg * 3 + 0], traj[min_value+step_height+1][leg * 3 + 1],
                                traj[min_value+step_height+1][leg * 3 + 1]))
        print("resultant step length is ", step_length, " for leg n° ", leg)
        # comp_rel(traj)
    print("\n")


############################################################################

t_move = 0
# test de normalized_move_xyz
if t_move:
    x0, y0, z0 = direct_leg(V[0], V[1], V[2], 0)
    print(x0, y0, z0)
    test_normalized_move_xyz(x0, y0, z0, 750, y0, z0, 1, 1, 5, kin.FL)
    test_normalized_move_xyz(x0, y0, z0, x0, 700, z0, 1, 1, 5, kin.FL)
    test_normalized_move_xyz(x0, y0, z0, x0, y0, -550, 1, 1, 5, kin.FL)
    test_normalized_move_xyz(x0, y0, z0, 750, 600, -400, 1, 1, 5, kin.FL)

t_jacob = 0
# test de l'erreur relative de position en appliquant la Jacobienne
if t_jacob:
    test_jacob_2(495, 585, 0, 0)
    test_jacob_2(505, 585, 0, 0)
    test_jacob_2(515, 585, 0, 0)

    test_jacob_2(495, 585, 10, 0)
    test_jacob_2(505, 585, 10, 0)
    test_jacob_2(515, 585, 10, 0)

    # test_jacob_2(495, 585, -10, 0)
    # test_jacob_2(505, 585, -10, 0)
    # test_jacob_2(515, 585, -10, 0)

    test_jacob_2(495, 585, 0, 10)
    test_jacob_2(505, 585, 0, 10)
    test_jacob_2(515, 585, 0, 10)

    # test_jacob_2(495, 585, 0, -10)
    # test_jacob_2(505, 585, 0, -10)
    # test_jacob_2(515, 585, 0, -10)

    # test_jacob_2(495, 585, 10, 10)
    # test_jacob_2(505, 585, 10, 10)
    # test_jacob_2(515, 585, 10, 10)

    # test_jacob_2(495, 585, 10, -10)
    # test_jacob_2(505, 585, 10, -10)
    # test_jacob_2(515, 585, 10, -10)

    # test_jacob_2(495, 585, -10, 10)
    # test_jacob_2(505, 585, -10, 10)
    # test_jacob_2(515, 585, -10, 10)

    # test_jacob_2(495, 585, -10, -10)
    # test_jacob_2(505, 585, -10, -10)
    # test_jacob_2(515, 585, -10, -10)

t_a = 0
if t_a:
    print(test_A(1, -1, 0.1, V[0], V[1], V[2]))

t_jacob_direct = 0
if t_jacob_direct:
    test_jacob_2_direct()
    test_precision_jacobienne_en_direct(0.495, 0.555, 0.015, -0.015, kin.FL)

test_comp_indirect = 0
if test_comp_indirect:
    test_comparaison_minimize_vs_jacob_indirect(0.535, 0.615, 0.001, -0.0015)
    test_comparaison_minimize_vs_jacob_indirect(0.535, 0.615, 0.01, -0.015)
    test_comparaison_minimize_vs_jacob_indirect(0.535, 0.615, 0.1, -0.15)
    test_comparaison_minimize_vs_jacob_indirect(0.535, 0.615, -0.01, +0.015)

t_artificial_inverse = 0
if t_artificial_inverse:
    traj_FL = draw_circle_rel_3(200, 200, 550, 600, 515, 0)
    traj_FR = draw_circle_rel_3(200, 200, 550, 600, 515, 1)
    for i in range(len(traj_FR)):
        traj_FR[i] = ROBOT['legs'][1]['matrix'] @ traj_FR[i]
    traj_RL = draw_circle_rel_3(200, 200, 550, 600, 515, 2)
    for i in range(len(traj_RL)):
        traj_RL[i] = ROBOT['legs'][2]['matrix'] @ traj_RL[i]
    traj_RR = draw_circle_rel_3(200, 200, 550, 600, 515, 3)
    for i in range(len(traj_RR)):
        traj_RR[i] = ROBOT['legs'][3]['matrix'] @ traj_RR[i]
    traj_4_legs = np.array([traj_FL, traj_FR, traj_RL, traj_RR])
    V = np.array([[550, 600, 515],
                  [550, 600, 515],
                  [550, 600, 515],
                  [550, 600, 515]])
    draw_move_4_legs(traj_4_legs, V)

t_different_moves = 0
if t_different_moves:
    # draw_move_leg(draw_circle_rel_3(150, 100, 550, 600, 515, kin.FL), 550, 600, 515, kin.FL, upgrade=False, solved=True)
    # draw_move_leg(draw_circle_rel_3(150, 200, 550, 600, 515, kin.FL), 550, 600, 515, kin.FL, upgrade=False, solved=False)
    draw_move_leg(draw_circle_rel_3(150, 100, 550, 600, 515, kin.FL), 550, 600, 515, kin.FL, upgrade=True, solved=True,
                  anim=False)
    # draw_move_leg(draw_circle_rel_3(150, 100, 550, 600, 515, kin.FL), 550, 600, 515, kin.FL, upgrade=True, solved=False)
    # test_circle_2(450, 500, 200, 200, kin.FL)
    # test_circle_3(550, 600, 515, 200, 200, kin.FL)
    # test_penalty_move_XZ(450, 500, 500, 10, kin.FL)

t_comp = 0
# comparaison en temps de nos modèles directs (v1 / Julien)
if t_comp:
    t = time.time()
    for i in range(10000):
        direct_v1(495, 585)
    t1 = time.time() - t

    t = time.time()
    for i in range(10000):
        kin.get_leg_points_V1_V2(495 / 1000, 585 / 1000, ROBOT['legs'][FR]['lengths'])['J']
    t2 = time.time() - t

    print("direct_v1 prend ", t1 * 100, " us")
    print("L'algo de Julien prend ", t2 * 100, " us")

    print("direct v1 retourne : ", direct_v1(495, 585))
    print("L'algo de Julien retourne :",
          kin.get_leg_points_V1_V2(495 / 1000, 585 / 1000, ROBOT['legs'][FR]['lengths'])['J'])

t_compute_traj = 1
if t_compute_traj:
    test_compute_traj((0, 0), 1, abs=True)
    test_compute_traj((1, 0), 1, abs=True)
    # test_compute_traj((0, 1), -1, abs=True)
    # test_compute_traj((0, -1), 0, abs=True)
    test_compute_traj((np.sqrt(3)/2, 1/2), 1, abs=True)

t_furthest = 0
if t_furthest:
    # test_furthest_pos((1, 0), 1)
    # test_furthest_pos((np.sqrt(3)/2, 1/2), 1)
    test_furthest_all_legs((np.sqrt(3)/2, 1/2), 1, 200, abs=True)
    # test_furthest_all_legs((1, 0), 1, 0)
    # test_furthest_all_legs((1, 0), 1, 200, abs=True)
    # test_furthest_all_legs((1, 0), 1, 25)
    # test_furthest_all_legs((1, 0), 1, 250, abs=True)
