import kinetic as kin
import numpy as np

weigths = {}
for p in 'ABCDEFGHIJO':
    weigths[p]={}
weigths['Q']={}
weigths['R']={}
weigths['A']['B']=2.0
weigths['A']['C']=1.6
weigths['A']['E']=3.6
weigths['E']['G']=4.0
weigths['B']['F']=3.45
weigths['G']['J']=3.65
weigths['D']['C']=5.3
weigths['H']['I']=5.3
weigths['R']['Q']=32
weigths['O']['Q']=80
tweigths = {}
for a in weigths:
    for b in weigths[a]:
        weigths[b][a]=weigths[a][b]

def leg_center_of_mass(v1, v2, weigths):
    points = kin.get_leg_points_V1_V2(v1, v2)
    lweigths = []
    total_weigth = 0
    application_points = []
    G=0
    for i in range(len(points)):
        pt1 = points.keys()[i]
        for j in range(i+1, len(points)):
            pt2 = points.keys()[j]
            if pt1 in weigths and pt2 in weigths['pt']:
                lweigths.append((weigths[pt1]['pt2'], pt1, pt2))
    for w, pt1, pt2 in lweigths:
        m = np.mean(pt1, pt2)
        application_points.append((m, w))
    for m, w in application_points:
        G += m
        total_weigth += w
    G[0] /= total_weigth
    return G, total_weigth