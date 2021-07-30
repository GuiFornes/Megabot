import numpy as np

from upg_tools import *
from upg_kinetic import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

v_init = [650, 650, 650]

def experiment_init():
    """ initie les valeurs des vérins dans la structure robot """
    V = [600] * 12
    # print(V[0:3])
    set_verins_12(V)


def elongation_to_delta_verins(V):
    """
    à partir d'une liste de longueur de vérins, retourne le delta de l'élongation par rapport à la position précédente.

    :param V: liste des élongations
    :return: liste des delta de déplacement des vérins
    """
    dV = []
    dV.append([0, 0, 0])
    for i in range(1 ,len(V)):
        dV.append([V[i][0]-V[i-1][0], V[i][1]-V[i-1][1], V[i][2]-V[i-1][2]])
    return dV


def delta_verins_to_victor_array(dV):
    """
    fonction de mise en forme pour envoyer les informations pour les expérimentations sur le MegaBot,
    i.e. en séparant la valeur de son signe.

    :param dV: delta d'élongation des vérins
    :return: delta absolu, liste des signes des delta
    """
    deltaV = []
    sens = []
    for i in range(len(dV)):
        sens_temp = []
        deltaV_temp = []
        for j in range(3):
            deltaV_temp.append(round(abs(dV[i][j]), 1))
            if dV[i][j] != 0 and dV[i][j]/abs(dV[i][j]) < 0:
                sens_temp.append(1)
            else:
                sens_temp.append(0)
        sens.append(sens_temp)
        deltaV.append(deltaV_temp)
    return deltaV, sens


def make_a_circle_leg(leg_id, r, n):
    """
    Créé une trajectoire circulaire pour une jambe et calcul les valeurs de vérins correspondant

    :param leg_id: ID de la patte
    :param r: rayon (en mm)
    :param n: nombre de points dans la trajectoire
    :return: liste des verins, delta absolu des vérin et signe des vérins
    """
    experiment_init()
    v1, v2, v3 = get_verins_3(leg_id)
    circle = draw_circle_rel_3(r, n, v1, v2, v3, leg_id)
    V = move_leg(circle, v1, v2, v3, leg_id)
    dV = elongation_to_delta_verins(V)
    # dV[0] = [v1 - v_init[0], v2 - v_init[1], v3 - v_init[2]]
    deltaV, sens = delta_verins_to_victor_array(dV)
    print("\n")
    for i in range(len(dV)):
        print(dV[i], deltaV[i], sens[i])
    print("\n")
    print(deltaV, sens)
    return V, deltaV, sens


def make_a_line(src, dst, dstep):
    """
    Créé une trajectoire linéaire pour une jambe

    :param src: position initiale
    :param dst: position cible
    :param dstep: distance entre 2 points de la trajectoire
    :return: trajectoire discrétisée
    """
    dist = distance(src, dst)
    traj = []
    for i in range(int(dist/dstep) + 1):
        traj.append(src + (dst - src) * i / int(dist/dstep))
    #print(traj)
    return traj


def make_a_M(leg_id, dstep):
    """
    Créé une trajectoire en forme de M pour une jambe et calcul les valeurs de vérins correspondant

    :param leg_id: ID de la patte
    :param dstep: distance entre 2 points de la trajectoire
    :return: liste des verins, delta absolu des vérin et signe des vérins
    """
    init()
    v1, v2, v3 = 475, 550, 550
    set_verins_3(v1, v2, v3, FL)
    v1, v2, v3 = get_verins_3(leg_id)
    x, y, z = direct_rel_3(v1, v2, v3, leg_id)
    traj = make_a_line(np.array([x, y, z]), np.array([x, y, z+200]), dstep)
    traj += make_a_line(np.array([x, y, z+200]), np.array([x-50, y+50, z+100]), dstep)
    traj += make_a_line(np.array([x-50, y+50, z+100]), np.array([x-100, y+100, z+200]), dstep)
    traj += make_a_line(np.array([x-100, y+100, z+200]), np.array([x-100, y+100, z]), dstep)
    traj += make_a_line(np.array([x-100, y+100, z]), np.array([x, y, z]), dstep)
    # for i in range(len(traj)):
    #     print(traj[i])
    V = move_leg(traj, v1, v2, v3, leg_id)
    dV = elongation_to_delta_verins(V)
    deltaV, sens = delta_verins_to_victor_array(dV)
    for i in range(len(deltaV)):
        print(deltaV[i])
    for i in range(len(deltaV)):
        print(sens[i])
    return V, deltaV, sens


def test_move_along_traj(circle=1):
    """
    trace l'évolution de la position de la patte en fonction de la liste des vérins entrée en paramètre

    :param circle: booléen, trace un cercle ou un 'M'
    :return:
    """
    if circle:
        V, deltaV, sens = make_a_circle_leg(FL, 200, 20)
    else:
        V, deltaV, sens = make_a_M(FL, 20)
    for i in range(len(V)):
        V[i] = V[i][0], V[i][1], V[i][2], 550, 550, 550, 550, 550, 550, 550, 550, 550
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


if __name__ == "__main__":
    # make_a_line(np.array([1100, 1100, 1100]), np.array([1100, 1100, 1200]), 10)
    # make_a_M(FL, 20)
    test_move_along_traj(circle=1)
