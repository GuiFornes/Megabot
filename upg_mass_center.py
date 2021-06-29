import kinetic as kin
from upg_tools import *
import numpy as np

# body_weights = {}
# for p in 'OP':
#     body_weights[p]={}
# body_weights['O']['O']=32
# body_weights['O']['P']=80
# for a in body_weights:
#     for b in body_weights[a]:
#         body_weights[b][a]=body_weights[a][b]

leg_weights = {}
for p in 'ABCDEFGHIJO':
    leg_weights[p]={}
leg_weights['A']['B']=2.0
leg_weights['A']['C']=1.6
leg_weights['A']['E']=3.6
leg_weights['E']['G']=4.0
leg_weights['B']['F']=3.45
leg_weights['G']['J']=3.65
leg_weights['D']['C']=5.3
leg_weights['H']['I']=5.3
tweights = {}
for a in leg_weights:
    for b in leg_weights[a]:
        leg_weights[b][a]=leg_weights[a][b]

def leg_center_of_mass(v1, v2, leg_id):
    """
    Return the center of mass of the leg leg_id in its referencial and the equivalent mass
    """
    points = kin.get_leg_points_V1_V2(v1 / 1000, v2 / 1000, kin.LEGS[leg_id]['lengths'])
    lweights = []
    total_weight = 0
    application_points = []
    G=0
    for i in range(len(points)):
        pt1 = list(points)[i]
        for j in range(i+1, len(points)):
            pt2 = list(points)[j]
            if pt1 in leg_weights and pt2 in leg_weights[pt1]:
                lweights.append((leg_weights[pt1][pt2], pt1, pt2))
    for w, pt1, pt2 in lweights:
        m = np.mean([points[pt1], points[pt2]], axis=0)
        application_points.append((m, w))
    for m, w in application_points:
        G += m * w
        total_weight += w
    G /= total_weight
    return G * 1000, total_weight

def body_center_of_mass(passenger=True, passenger_weight=80.0):
    """
    Return the center of mass of the body of the MegaBot in its referencial and the equivalent mass
    Can take a passenger in argument
    """
    if passenger:
        total_weight = 32.0 + passenger_weight
        return [0.0, 0.0, 750.0 * passenger_weight] / total_weight, total_weight
    return [0.0, 0.0], 32.0

V = get_verins_12()
print(leg_center_of_mass(V[0], V[1], 0))