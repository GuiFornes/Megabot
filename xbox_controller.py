import pygame.joystick
import time
from upg_planning import *

def get_joystick():
    d = (- pygame.joystick.Joystick(0).get_axis(1), - pygame.joystick.Joystick(0).get_axis(0))
    r = - pygame.joystick.Joystick(0).get_axis(3)
    return d, r

if __name__ == "__main__":
    print(pygame.joystick.get_init())
    pygame.joystick.init()
    print(pygame.joystick.get_init())
    print(pygame.joystick.get_count())
    pygame.joystick.Joystick(0).init()
    while True:
        d, r = get_joystick()
        print(d, r)
        # traj = compute_traj_form_joystick_abs(d, r)
        if  0 == 1:# Tracé
            fig = plt.figure()
            ax = fig.gca(projection='3d')  # Affichage en 3D
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
            ax.plot(Xt0, Yt0, Zt0, label='Théorique', c='coral', ms=0.1)  # Tracé de la courbe théorique
            ax.plot(Xt1, Yt1, Zt1, label='Théorique', c='cyan', ms=0.1)  # Tracé de la courbe théorique
            ax.plot(Xt2, Yt2, Zt2, label='Théorique', c='deeppink', ms=0.1)  # Tracé de la courbe théorique
            ax.plot(Xt3, Yt3, Zt3, label='Théorique', c='chartreuse', ms=0.1)  # Tracé de la courbe théorique
            plt.title("Trajectoire du bout de la patte dans l'espace")
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            ax.set_xbound(-2000, 2000)
            ax.set_ybound(-2000, 2000)
            ax.set_zbound(2000, -2000)
            plt.show()
        time.sleep(1)
    # pygame.joystick.Joystick(0).quit()
