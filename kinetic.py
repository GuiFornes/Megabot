#!/usr/bin/python
from tkinter import *

import math
import sys
import numpy

from geo import *

FL=0 # front left leg
FR=1 # front right leg
RL=2 # rear left leg
RR=3 # rear right leg

ALL_LEGS=(FL,FR,RL,RR)

BODY_FRAME=1.0

# leg angle is relative to x axis

LEGS={FL:{'origin':(-BODY_FRAME/2.0,BODY_FRAME/2.0,0),                          \
          'angle':3.0*math.pi/4.0,                                              \
          'angle_orig':3.0*math.pi/4.0,                                         \
          'matrix': [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],                  \
          'imatrix':[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],                  \
          'lengths':{'ao':135.0,'bo':120.0,'bcx':290.0,'bcy':60.0,              \
                     'ae':500.0,'de':100.0,'ef':450.0,'fg':300.0,               \
                     'fh':200.0,'gi':520.0,'bf':600.0,'gj':1055.0,              \
                     'yaw_a':700.0,'yaw_b':50.0,'yaw_c':280.0}               }, \
      FR:{'origin':(BODY_FRAME/2.0,BODY_FRAME/2.0,0),          \
          'angle':math.pi/4.0,           \
          'angle_orig':math.pi/4.0,      \
          'matrix': [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],                  \
          'imatrix':[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],                  \
          'lengths':{'ao':130.0,'bo':120.0,'bcx':300.0,'bcy':60.0,              \
                     'ae':500.0,'de':100.0,'ef':445.0,'fg':285.0,               \
                     'fh':200.0,'gi':500.0,'bf':603.0,'gj':1035.0,              \
                     'yaw_a':700.0,'yaw_b':55.0,'yaw_c':280.0}               }, \
      RL:{'origin':(-BODY_FRAME/2.0,-BODY_FRAME/2.0,0),        \
          'angle':5.0*math.pi/4.0,       \
          'angle_orig': 5.0*math.pi/4.0, \
          'matrix': [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],                  \
          'imatrix':[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],                  \
          'lengths':{'ao':130.0,'bo':120.0,'bcx':295.0,'bcy':60.0,              \
                     'ae':495.0,'de':100.0,'ef':450.0,'fg':300.0,               \
                     'fh':200.0,'gi':515.0,'bf':600.0,'gj':1055.0,              \
                     'yaw_a':700.0,'yaw_b':60.0,'yaw_c':280.0}               }, \
      RR:{'origin':(BODY_FRAME/2.0,-BODY_FRAME/2.0,0),         \
          'angle':7.0*math.pi/4.0,       \
          'angle_orig':7.0*math.pi/4.0,  \
          'matrix': [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],                  \
          'imatrix':[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],                  \
          'lengths':{'ao':130.0,'bo':120.0,'bcx':290.0,'bcy':60.0,              \
                     'ae':495.0,'de':100.0,'ef':445.0,'fg':300.0,               \
                     'fh':200.0,'gi':500.0,'bf':600.0,'gj':1045.0,              \
                     'yaw_a':700.0,'yaw_b':55.0,'yaw_c':280.0}               }  \
}

ROBOT={'LEGS':LEGS}




LEG_PARTS_LENGTHS={'ao':130.0,'bo':120.0,'bcx':290.0,'bcy':60.0,'ae':500.0,'de':100.0,'ef':450.0,'fg':300.0,'fh':200.0,'gi':520.0,'bf':600.0,'gj':1030.0}


def update_transformation_matrices_for(matrix,imatrix,angle,origin):
    """
    Update transformation matrices (rotations) of the struct LEGS

    :param matrix: rotation matrix that will be modified
    :param imatrix: inverse matrix which will be too
    :param angle: angle between leg and body-frame
    :param origin: position of origin of the leg
    :return:
    """
    m = [ [math.cos(angle) , 0 , -math.sin(angle) , origin[X] ],
          [math.sin(angle) , 0 , -math.cos(angle) , origin[Y] ],
          [              0 , 1 ,                0 , origin[Z] ],
          [              0 , 0 ,                0 , 1 ] ]
    
    im = [ [math.cos(-angle) , -math.sin(-angle) , 0 , -origin[X] * math.cos(-angle) + origin[Y] * math.sin(-angle)],
           [               0 ,                0 , 1 , -origin[Z] ],
           [math.sin(angle)  , -math.cos(angle) , 0 , origin[X] * math.sin(-angle) + origin[Y] * math.cos(-angle) ],
           [               0 ,                0 , 0 , 1 ] ]
    for l in range(len(matrix)):
        for c in range(len(matrix[0])):
            matrix[l][c]=m[l][c]
            imatrix[l][c]=im[l][c]
                   



    
# compute transformation matrices
def update_transformation_matrices(legs,leg_list=ALL_LEGS):
    """update all matrices of LEGS calling its sub-function above"""
    for l in leg_list:
        update_transformation_matrices_for(legs[l]['matrix'],legs[l]['imatrix'],legs[l]['angle'],legs[l]['origin'])



update_transformation_matrices(ROBOT['LEGS'])


def get_leg_points(V1,V2,AO,BO,BCx,BCy,AE,DE,EF,FG,FH,GI,BF,GJ):
    """Annexe function of get_leg_points_V1_V2 that should not be call alone"""
    #print 'get leg width V1: ',V1,' and V2: ',V2
    points={'O':(0.0,0.0),'A':(-AO,0.0),'B':(BO,0),'C':(BCx+BO,BCy)} #,'Q':(-0.8,0),'R':(-0.8-0.79,-0.49)}

    points['D'],p2= compute(points['A'],points['C'], AE-DE , V1) 
    points['E'] = one_line(points['A'], points['D'] , AE)
    points['F'],p2= compute(points['E'],points['B'], EF , BF)
    points['H'] = one_line(points['B'], points['F'] , BF-FH)
    points['G'] = one_line(points['E'], points['F'] , EF+FG)
    points['I'],p2= compute(points['G'],points['H'], GI , V2)
    points['J'] = one_line(points['G'], points['I'] , GJ)
    return points

def get_leg_points_V1_V2(v1,v2,lpl):
    """
    Return the list of coords of each points of the leg in the leg referentiel (2 dimensions)

    :param v1: elongation of v1
    :param v2: elongation of v2
    :param lpl: list o all constant lengths
    :return: list of all coords
    """
    return get_leg_points(v1,
                          v2,
                          lpl['ao']/1000.0,
                          lpl['bo']/1000.0,
                          lpl['bcx']/1000.0,
                          lpl['bcy']/1000.0,
                          lpl['ae']/1000.0,
                          lpl['de']/1000.0,
                          lpl['ef']/1000.0,
                          lpl['fg']/1000.0,
                          lpl['fh']/1000.0,
                          lpl['gi']/1000.0,
                          lpl['bf']/1000.0,
                          lpl['gj']/1000.0)

# return v3 in mm
def delta_angle_to_linear_actuator_v3(delta,lpl):
    """
    Compute the elongation needed to have delta angle between the leg and the body-frame

    :param delta: angle between leg and frame
    :param lpl: list of all constant lengths
    :return: elongation of v3 (in mm)
    """
    na=math.pi/4.0+delta - math.asin(lpl['yaw_b']/lpl['yaw_a'])
    v3=math.sqrt(lpl['yaw_c']**2+lpl['yaw_a']**2-2*lpl['yaw_c']*lpl['yaw_a']*math.cos(na))
    return v3

# v3 suppose to be in mm
def v3_to_delta_angle(v3,lpl):
    """
    Compute the angle between the leg and the body-frame

    :param v3: elongation of v3 (in mm)
    :param lpl: list of all constant lengths
    :return: angle in degrees
    """
    na=math.acos((lpl['yaw_c']**2+lpl['yaw_a']**2-v3**2)/(2*lpl['yaw_c']*lpl['yaw_a']))
    return -(na+math.asin(lpl['yaw_b']/lpl['yaw_a'])-math.pi/4.0)

#v3 suppose to be in m
def update_angle(legs,leg,v3):
    """
    Update the angle of 'leg' in the struct LEGS

    :param legs: struct LEGS
    :param leg: leg id
    :param v3: elongation of v3 in m
    :return:
    """
    legs[leg]['angle']=legs[leg]['angle_orig'] + v3_to_delta_angle(v3*1000.0,legs[leg]['lengths'])
    update_transformation_matrices(legs,(leg,))



    
from scipy.optimize import minimize


def inverse_kinetic_fun(v,dx,dy,lpl):
    """ complementary function to give to 'inverse_kinetic' """
    for k in v:
        if k<0.449 or k>0.651:
            print("!!!!! out of bounds",k)
            return None
    p=get_leg_points_V1_V2(v[0],v[1],lpl)
    return norm(p2v(p['J'],(dx,dy)))


def inverse_kinetic(x,y,lpl):
    """
    compute the 2 dimensions approximation of inverse kinematic

    :param x: asked x position
    :param y: asked y position
    :param lpl: all the constant lengths of leg
    :return: v1 and v2 elongations, all final points, deviation from position
    """
    res=minimize(lambda v,dx,dy: inverse_kinetic_fun(v,dx,dy,lpl),(0.55,0.55),args=(x,y),\
                 #method='Nelder-Mead',\
                 bounds=((0.45,0.65),(0.45,0.65)),
                 tol=0.0001)
    p=get_leg_points_V1_V2(res.x[0],res.x[1],lpl)
    return (res.x[0],res.x[1],p,norm(p2v(p['J'],(x,y))))



def inverse_kinetic_robot_ref(legs,leg,point):
    """
    approximate inverse kinematic using minimize module to compute the necessary elongation of V
    to place the end of the 'leg' at 'point'.

    :param legs: struct LEGS
    :param leg: leg id
    :param point: destination point
    :return: err:boolean is this movement possible? (tuple):elongation of v1, v2, v3
    """
    # compute angle for v3
    error=False
    dx=point[X]-legs[leg]['origin'][X]
    dy=point[Y]-legs[leg]['origin'][Y]
    angle=math.atan2(dy,dx) # angle between x axis and (leg,point) 
    if angle<0:
        angle+=math.pi*2.0
    v3=delta_angle_to_linear_actuator_v3(legs[leg]['angle_orig'] - angle,legs[leg]['lengths'])
    v3=v3/1000.0
    if v3<0.45 or v3>0.65:
        error=True
        print("impossible angle for v3",leg," delta:",math.degrees(legs[leg]['angle_orig'] - angle),' => v3:',v3)
    legs[leg]['angle']=angle
    update_transformation_matrices(legs,(leg,))    
    point_leg_ref=numpy.matmul(legs[leg]['imatrix'],point+[1]).tolist()
    r=inverse_kinetic(point_leg_ref[X],point_leg_ref[Y],legs[leg]['lengths'])
    if r[3]>0.1:
        error=True
        print("too far!")
    print(r[0], r[1], v3)
    return error,(r[0],r[1],v3) # v3 -> A

def robot_ref_leg_point(legs,leg,v1,v2,v3):
    """
    Compute the position of the end of the 'leg' in the ref of the robot using direct kinematic

    :param legs: struct legs
    :param leg: leg id
    :param v1: v1 elongation
    :param v2: v2 elongation
    :param v3: v3 elongation
    :return: (tuple) position of the leg
    """
    p=get_leg_points_V1_V2(v1,v2,legs[leg]['lengths'])
    # p['J'] is the x,y coords in leg referencial
    a=(p['J'][X],p['J'][Y],0,1)
    update_angle(legs,leg,v3)
    b=numpy.matmul(legs[leg]['matrix'],a).tolist()
    return (b[X],b[Y],b[Z])

def robot_ref_leg_points(legs,leg,v1,v2,v3):
    """
    Compute the position of all points of the 'leg' in the ref of the robot using direct cinematic

    :param legs: struct legs
    :param leg: leg id
    :param v1: v1 elongation
    :param v2: v2 elongation
    :param v3: v3 elongation
    :return: (list of tuple) position of all points of the leg
    """
    p=get_leg_points_V1_V2(v1,v2,legs[leg]['lengths'])
    # p['J'] is the x,y coords in leg referencial
    update_angle(legs,leg,v3)
    q={}
    for i in list(p.keys()):
        a=(p[i][X],p[i][Y],0,1)
        b=numpy.matmul(legs[leg]['matrix'],a).tolist()
        q[i]=(b[X],b[Y],b[Z])
    return q
