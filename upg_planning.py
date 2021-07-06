from upg_kinetic import *
from upg_tools import *
import numpy as np

MIN_RADIUS = 2500
MAX_RADIUS = 10000

INF = 100000


def cmd_joystick(D, R):
    """
    Traite les données reçues du controleur joystick pour les convertir en une translation dirigée et une rotation

    :param D: couple de coord du joystick de déplacement
    :param R: facteur de rotation du joystick de rotation
    :return: position du centre, direction de marche
    >>> cmd_joystick((1, 0), 0)
    ((0, 100000), (1, 0))
    >>> cmd_joystick((0, 1), 0)
    ((100000, 0), (0, 1))
    >>> cmd_joystick((0.51, 0.49), 0)
    ((0, 100000), (1, 0))
    >>> cmd_joystick((0.49, 0.51), 1)
    ((2499.999999999999, 0.0), (0, 1))

    """
    if D[0] == 0 and D[1] == 0:
        if abs(R) < 0.1:
            # l'immobilité
            return "erreur, pas de mouvement en entrée"
        # Rotation sur lui-même
        return (0, 0), (0, 0)

    # Direction NSEW
    if abs(D[0]) >= abs(D[1]):
        direction = (D[0] > 0) - (D[0] < 0), 0
    else:
        direction = 0, (D[1] > 0) - (D[1] < 0)

    if abs(R) < 0.1:
        # Marche rectiligne
        centre = direction[1] * INF, direction[0] * INF
        return centre, direction

    signeD = (R > 0) - (R < 0)
    facteur = signeD * -7500 / 0.9 * abs(R) + 10000 + 750 / 0.9
    centre = facteur * direction[1], facteur * direction[0]
    return centre, direction


############################################################################
if __name__ == "__main__":
    import doctest
    doctest.testmod()

######################## DEPRECATED #############################
def length(a, b):
    return distance(b[0], b[1], x2=a[0], y2=a[1])

def angle_between_linear(coef1, coef2):
    pointO = 0, 0
    pointA = 1, coef1
    pointB = 1, coef2
    return al_kashi_angle(length(pointO, pointA), length(pointO, pointB), length(pointA, pointB))

def to_ref_traj_next_step(alpha, dx, dy):
    return np.array([
        [np.cos(alpha), - np.sin(alpha), dx],
        [np.sin(alpha), np.cos(alpha), dy],
        [0, 0, 1]
    ])

def compute_trajs(traj):
    all_trajs, trajFL, trajFR, trajRL, trajRR = [], [], [], [], []
    V = get_verins_12()
    Pos = direct_12(V)
    fl = (Pos[0], Pos[1], Pos[2])
    fr = (Pos[3], Pos[4], Pos[5])
    last_coef = (fr[1] - fl[1]) / (fr[0] - fl[0])
    if last_coef > 0.1:
        print("Legs position aren't initialized")
    legs_spacing = (abs(fr[1]) + abs(fl[1])) / 2

    for i in range(1, len(traj)):
        direction = (traj[1][0] - traj[0][0]) / (traj[1][1] - traj[0][1])
        current_coef = -1 / direction
        alpha = angle_between_linear(last_coef, current_coef)
        m_rota = to_ref_traj_next_step(alpha, traj[1][0] - traj[0][0], traj[1][1] - traj[0][1])
        step_fl = m_rota @ np.array([-legs_spacing, 0, 1])
        step_fr = m_rota @ np.array([legs_spacing, 0, 1])
        trajFL.append(np.array([step_fl[0], step_fl[1], traj[i][i % 4][2]]))
        trajFR.append(np.array([step_fr[0], step_fr[1], traj[i][i % 4][2]]))

    trajRL.append(draw_line_3(V[6], V[7], V[8], Pos[0] - Pos[6], Pos[1] - Pos[7], Pos[2] - Pos[8], 20, 2))
    trajRR.append(draw_line_3(V[9], V[10], V[11], Pos[3] - Pos[9], Pos[4] - Pos[10], Pos[5] - Pos[11], 20, 2))
    for i in range(len(trajFL) - len(trajRL)):
        trajRL.append(trajFL[i])
        trajRR.append(trajFR[i])

    for i in range(len(trajFL)):
        traj.append(trajFL[i])
        traj.append(trajFR[i])
        traj.append(trajRL[i])
        traj.append(trajRR[i])

