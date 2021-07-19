import numpy as np
from numpy.linalg import inv
from numpy import dot
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d  # Fonction pour la 3D
import matplotlib.animation as animation
import time

import kinetic
from upg_kinetic import *
from upg_tools import *
from upg_deprecated import *
from upg_jacobian import *
from upg_planning import *
from upg_mass_center import *


############################# FONCTIONS DE TEST ################################


def draw_abs(LV, LO, LOmega):
    """
    Trace la trajectoire des extrémités des pattes du robot dans le repère absolu
    """
    # Elongations des vérins
    V1 = [v[0] for v in LV]
    V2 = [v[1] for v in LV]
    V3 = [v[2] for v in LV]
    V4 = [v[3] for v in LV]
    V5 = [v[4] for v in LV]
    V6 = [v[5] for v in LV]
    V7 = [v[6] for v in LV]
    V8 = [v[7] for v in LV]
    V9 = [v[8] for v in LV]
    V10 = [v[9] for v in LV]
    V11 = [v[10] for v in LV]
    V12 = [v[11] for v in LV]
    T = np.linspace(0, 1, len(LV))

    # Positions des bouts de patte
    Pos = list(map(direct_abs, LV, LO, LOmega))
    Xp0 = [p[0] for p in Pos]
    Yp0 = [p[1] for p in Pos]
    Zp0 = [p[2] for p in Pos]
    Xp1 = [p[3] for p in Pos]
    Yp1 = [p[4] for p in Pos]
    Zp1 = [p[5] for p in Pos]
    Xp2 = [p[6] for p in Pos]
    Yp2 = [p[7] for p in Pos]
    Zp2 = [p[8] for p in Pos]
    Xp3 = [p[9] for p in Pos]
    Yp3 = [p[10] for p in Pos]
    Zp3 = [p[11] for p in Pos]

    # Tracé des trajectoires
    fig = plt.figure()
    ax = fig.gca(projection='3d')  # Affichage en 3D
    ax.scatter(Xp0, Yp0, Zp0, label='Positions', marker='.', s=3, c='red')  # Tracé des points réels
    ax.scatter(Xp1, Yp1, Zp1, label='Positions', marker='.', s=3, c='blue')  # Tracé des points réels
    ax.scatter(Xp2, Yp2, Zp2, label='Positions', marker='.', s=3, c='purple')  # Tracé des points réels
    ax.scatter(Xp3, Yp3, Zp3, label='Positions', marker='.', s=3, c='green')  # Tracé des points réels
    plt.title("Trajectoire du bout de la patte dans l'espace")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xbound(-2000, 2000)
    ax.set_ybound(-2000, 2000)
    ax.set_zbound(2000, -2000)
    plt.show()

    # Tracé des élongations des vérins au cours du temps
    plt.plot(T, V1, label='V1')
    plt.plot(T, V2, label='V2')
    plt.plot(T, V3, label='V3')
    plt.plot(T, V4, label='V4')
    plt.plot(T, V5, label='V5')
    plt.plot(T, V6, label='V6')
    plt.plot(T, V7, label='V7')
    plt.plot(T, V8, label='V8')
    plt.plot(T, V9, label='V9')
    plt.plot(T, V10, label='V10')
    plt.plot(T, V11, label='V11')
    plt.plot(T, V12, label='V12')
    plt.title("Elongations des vérins dans le mouvement")
    plt.show()


def draw_rel(V):
    # Positions des bouts de patte
    Pos = list(map(direct_rel_12, V))
    Xp0 = [p[0] for p in Pos]
    Yp0 = [p[1] for p in Pos]
    Zp0 = [p[2] for p in Pos]
    Xp1 = [p[3] for p in Pos]
    Yp1 = [p[4] for p in Pos]
    Zp1 = [p[5] for p in Pos]
    Xp2 = [p[6] for p in Pos]
    Yp2 = [p[7] for p in Pos]
    Zp2 = [p[8] for p in Pos]
    Xp3 = [p[9] for p in Pos]
    Yp3 = [p[10] for p in Pos]
    Zp3 = [p[11] for p in Pos]

    # Tracé des trajectoires
    fig = plt.figure()
    ax = fig.gca(projection='3d')  # Affichage en 3D
    ax.scatter(Xp0, Yp0, Zp0, label='Positions', marker='.', s=3, c='red')  # Tracé des points réels
    ax.scatter(Xp1, Yp1, Zp1, label='Positions', marker='.', s=3, c='blue')  # Tracé des points réels
    ax.scatter(Xp2, Yp2, Zp2, label='Positions', marker='.', s=3, c='purple')  # Tracé des points réels
    ax.scatter(Xp3, Yp3, Zp3, label='Positions', marker='.', s=3, c='green')  # Tracé des points réels
    plt.title("Trajectoire du bout de la patte dans l'espace")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xbound(-2000, 2000)
    ax.set_ybound(-2000, 2000)
    ax.set_zbound(2000, -2000)
    plt.show()

    # Tracé des élongations des vérins au cours du temps
    T = np.linspace(0, 1, len(V))
    plt.plot(T, [v[0] for v in V], label='V1')
    plt.plot(T, [v[1] for v in V], label='V2')
    plt.plot(T, [v[2] for v in V], label='V3')
    plt.plot(T, [v[3] for v in V], label='V4')
    plt.plot(T, [v[4] for v in V], label='V5')
    plt.plot(T, [v[5] for v in V], label='V6')
    plt.plot(T, [v[6] for v in V], label='V7')
    plt.plot(T, [v[7] for v in V], label='V8')
    plt.plot(T, [v[8] for v in V], label='V9')
    plt.plot(T, [v[9] for v in V], label='V10')
    plt.plot(T, [v[10] for v in V], label='V11')
    plt.plot(T, [v[11] for v in V], label='V12')
    plt.title("Elongations des vérins dans le mouvement")
    plt.show()


def comp_rel(traj):
    """
    Trace la trajectoire des extrémités des pattes du robot suivant traj avec move_rel
    """
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

    # Elongations des vérins
    Ver = move_rel(traj, V)
    V1 = [v[0] for v in Ver]
    V2 = [v[1] for v in Ver]
    V3 = [v[2] for v in Ver]
    V4 = [v[3] for v in Ver]
    V5 = [v[4] for v in Ver]
    V6 = [v[5] for v in Ver]
    V7 = [v[6] for v in Ver]
    V8 = [v[7] for v in Ver]
    V9 = [v[8] for v in Ver]
    V10 = [v[9] for v in Ver]
    V11 = [v[10] for v in Ver]
    V12 = [v[11] for v in Ver]
    T = np.linspace(0, 1, len(Ver))

    # Positions des bouts de patte
    Pos = list(map(direct_rel_12, Ver))
    Xp0 = [p[0] for p in Pos]
    Yp0 = [p[1] for p in Pos]
    Zp0 = [p[2] for p in Pos]
    Xp1 = [p[3] for p in Pos]
    Yp1 = [p[4] for p in Pos]
    Zp1 = [p[5] for p in Pos]
    Xp2 = [p[6] for p in Pos]
    Yp2 = [p[7] for p in Pos]
    Zp2 = [p[8] for p in Pos]
    Xp3 = [p[9] for p in Pos]
    Yp3 = [p[10] for p in Pos]
    Zp3 = [p[11] for p in Pos]

    # Tracé des trajectoires
    fig = plt.figure()
    ax = fig.gca(projection='3d')  # Affichage en 3D
    ax.plot(Xt0, Yt0, Zt0, label='Théorique', c='coral')  # Tracé de la courbe théorique
    ax.scatter(Xt0[0], Yt0[0], Zt0[0], c='black', s=60)
    ax.scatter(Xp0, Yp0, Zp0, label='Positions', marker='.', s=3, c='red')  # Tracé des points réels
    ax.plot(Xt1, Yt1, Zt1, label='Théorique', c='cyan')  # Tracé de la courbe théorique
    ax.scatter(Xt1[0], Yt1[0], Zt1[0], c='black', s=60)
    ax.scatter(Xp1, Yp1, Zp1, label='Positions', marker='.', s=3, c='blue')  # Tracé des points réels
    ax.plot(Xt2, Yt2, Zt2, label='Théorique', c='deeppink')  # Tracé de la courbe théorique
    ax.scatter(Xt2[0], Yt2[0], Zt2[0], c='black', s=60)
    ax.scatter(Xp2, Yp2, Zp2, label='Positions', marker='.', s=3, c='purple')  # Tracé des points réels
    ax.plot(Xt3, Yt3, Zt3, label='Théorique', c='chartreuse')  # Tracé de la courbe théorique

    ax.scatter(Xt3[0], Yt3[0], Zt3[0], c='black', s=60)
    ax.scatter(Xt3[6], Yt3[6], Zt3[6], c='blue', s=60)
    ax.scatter(Xp3[6], Yp3[6], Zp3[6], c='red', s=60)
    ax.scatter(Xt1[0], Yt1[0], Zt1[0], c='black', s=60)
    ax.scatter(Xt1[6], Yt1[6], Zt1[6], c='blue', s=60)
    ax.scatter(Xp1[6], Yp1[6], Zp1[6], c='red', s=60)

    ax.scatter(Xp3, Yp3, Zp3, label='Positions', marker='.', s=3, c='green')  # Tracé des points réels
    plt.title("Trajectoire du bout de la patte dans l'espace")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xbound(-2000, 2000)
    ax.set_ybound(-2000, 2000)
    ax.set_zbound(2000, -2000)
    plt.show()

    # Tracé des élongations des vérins au cours du temps
    plt.plot(T, V1, label='V1')
    plt.plot(T, V2, label='V2')
    plt.plot(T, V3, label='V3')
    plt.plot(T, V4, label='V4')
    plt.plot(T, V5, label='V5')
    plt.plot(T, V6, label='V6')
    plt.plot(T, V7, label='V7')
    plt.plot(T, V8, label='V8')
    plt.plot(T, V9, label='V9')
    plt.plot(T, V10, label='V10')
    plt.plot(T, V11, label='V11')
    plt.plot(T, V12, label='V12')
    plt.title("Elongations des vérins dans le mouvement")
    plt.show()


def test_zone_accessible(leg_id):
    fig = plt.figure()
    ax = fig.gca(projection='3d')  # Affichage en 3D
    for x in range(300, 1701, 100):
        for y in range(300, 1701, 100):
            for z in range(-100, -851, -100):
                if not is_accessible(leg_id, (x, y, z))[0]:
                    x=x #ax.scatter(x, y, z, marker='.', s=30, c='grey')
                else:
                    ax.scatter(x, y, z, marker='.', s=160, c='green')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    # ax.set_xbound(-1000, 2000)
    # ax.set_ybound(-1000, 2000)
    # ax.set_zbound(-1600, -200)
    plt.title("points accessibles pour la patte")
    plt.show()
    middle = (1080, 1065, -530)
    print(is_accessible(leg_id, middle))

def test_zone_accessible_abs(leg_id):
    init()
    V = get_verins_12()
    x, y, z = direct_abs(V)

    fig = plt.figure()
    ax = fig.gca(projection='3d')  # Affichage en 3D


def test_compute_traj(R, D):
    traj = compute_traj_form_joystick(cmd_joystick(R, D))
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
    ax.plot(Xt0, Yt0, Zt0, label='Théorique', c='coral')  # Tracé de la courbe théorique
    ax.plot(Xt1, Yt1, Zt1, label='Théorique', c='cyan')  # Tracé de la courbe théorique
    ax.plot(Xt2, Yt2, Zt2, label='Théorique', c='deeppink')  # Tracé de la courbe théorique
    ax.plot(Xt3, Yt3, Zt3, label='Théorique', c='chartreuse')  # Tracé de la courbe théorique
    plt.title("Trajectoire du bout de la patte dans l'espace")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.plot((-500, -500, 500, 500, -500), (-500, 500, 500, -500, -500))
    pos = direct_rel_12(V)
    ax.plot((500, pos[0]), (500, pos[1]), (0, pos[2]), c='red')
    ax.plot((500, pos[3]), (-500, pos[4]), (0, pos[5]), c='blue')
    ax.plot((-500, pos[6]), (500, pos[7]), (0, pos[8]), c='purple')
    ax.plot((-500, pos[9]), (-500, pos[10]), (0, pos[11]), c='green')
    ax.set_xbound(-2000, 2000)
    ax.set_ybound(-2000, 2000)
    ax.set_zbound(-1000, 1000)
    for i in range(4):
        ax.scatter(traj[0][i*3+0], traj[0][i*3+1], traj[0][i*3+2], c='black', s=60)
    plt.show()


def test_furthest_pos(D, R):
    traj = compute_traj_form_joystick(cmd_joystick(D, R))
    max_step = []
    for leg in range(4):
        step = furthest_accessible(traj, leg)
        max_step.append(step)
    print("maximum point of the trajectory reached for each leg : ", max_step)
    min_value = min(max_step)
    leg = max_step.index(min_value)
    step_length = distance((traj[0][leg * 3 + 0], traj[0][leg * 3 + 1], traj[0][leg * 3 + 1]),
                           (traj[min_value][leg * 3 + 0], traj[min_value][leg * 3 + 1], traj[min_value][leg * 3 + 1]))
    print("resultant step length is ", step_length, " for leg n° ", leg)
    comp_rel(traj)


def test_furthest_all_legs(D, R):
    traj = compute_traj_form_joystick(cmd_joystick(D, R))
    max_step = furthest_accessible_real_step_all_legs(traj, 30)
    print("maximum point of the trajectory reached for each leg : ", max_step)
    min_value = min(max_step)
    leg = max_step.index(min_value)
    step_length = distance((traj[0][leg * 3 + 0], traj[0][leg * 3 + 1], traj[0][leg * 3 + 1]),
                           (traj[min_value][leg * 3 + 0], traj[min_value][leg * 3 + 1], traj[min_value][leg * 3 + 1]))
    print("resultant step length is ", step_length, " for leg n° ", leg)
    comp_rel(traj)


def test_get_last_leg():
    init()
    z_ini = direct_rel_12(get_verins_12())[2]
    set_verins_3(600, 500, 500, 0)
    pos_rel = direct_rel_12(get_verins_12())
    pos_abs = [pos_rel[0], pos_rel[1], - (z_ini - pos_rel[2])]
    set_leg_pos(pos_abs, 0)
    set_og(0, 0)
    for i in range(4):
        print(ROBOT['legs'][i]['pos_abs'])
    print(get_last_leg(0))


################################# TESTS ####################################

t_init = 0
if t_init:
    init()
    print(get_X())
    print(get_og(0))
    print(get_og(1))
    print(get_og(2))
    print(get_og(3))
    print(get_O())
    print(get_omega())

t_accessible = 0
if t_accessible:
    init()
    test_zone_accessible(FL)

t_compute_traj = 1
if t_compute_traj:
    init()
    # test_compute_traj((0, 0), 1)
    # test_compute_traj((1, 0), 1)
    # test_compute_traj((0, 1), -1)
    # test_compute_traj((0, -1), 0)
    # test_compute_traj((np.sqrt(3)/2, 1/2), 1)
    test_furthest_pos((1, 0), 1)
    test_furthest_all_legs((1, 0), 1)

t_abs_1 = 1
if t_abs_1:
    init()
    # set_og(0, 0)
    traj = traj_abs_sin_1(200, 100, 0)
    LV, LO, LOmega = move_abs_one_leg(traj, 0)
    draw_abs(LV, LO, LOmega)

t_abs_4 = 0
if t_abs_4:
    init()
    # set_og(0, 0)
    traj = traj_abs_sin_4(50, 100, 0)
    LV, LO, LOmega = move_abs_all_legs(traj)
    draw_abs(LV, LO, LOmega)

    traj = traj_abs_sin_4(300, 300, 0)
    LV, LO, LOmega = move_abs_all_legs(traj)
    draw_abs(LV, LO, LOmega)

    traj = traj_abs_sin_4(300, 600, 0)
    LV, LO, LOmega = move_abs_all_legs(traj, max_omega=45)
    draw_abs(LV, LO, LOmega)

    traj = traj_abs_sin_4(300, 600, 3)
    LV, LO, LOmega = move_abs_all_legs(traj, max_omega=45)
    draw_abs(LV, LO, LOmega)

    traj = M_letter(100, direct_abs(get_verins_12(), get_O(), get_omega()))
    LV, LO, LOmega = move_abs_all_legs(traj, max_omega=45)
    draw_abs(LV, LO, LOmega)

t_rel = 0
if t_rel:
    init()
    traj = draw_circle_rel_12(30, 200, get_verins_12())
    comp_rel(traj)
    draw_rel(shake_dat_ass_rel(30, 200))

t_com = 0
if t_com == 1:
    init()
    print(proj_center_of_mass_robot())