from math import sqrt

X=0
Y=1
Z=2


def norm2(v):
    n=0
    for x in v:
        n+=x**2
    return n
def norm(v):
    return sqrt(norm2(v))

def normed(v):
    n=norm(v)
    assert(n!=0)
    return tuple([a/n for a in v])

def dist(a,b):
    return norm(((b[X]-a[X]),(b[Y]-a[Y])))

def p2v(a,b):
    return tuple ([b[i] - a[i] for i in range(len(a))])

def plus(v1,v2):
    return tuple ([ v1[i] + v2[i] for i in range(len(v1)) ])

def prod(v1,a):
    return tuple([i*a for i in v1])
    

def one_line(centre, point , distance):
    cp=p2v(centre,point)
    d=normed(cp)
    return plus(centre,prod(d,distance))

def scalaire(v1,v2):
    n=0
    for i in range(len(v1)):
        n+=v1[i]*v2[i]
    return n

def vectoriel(v1,v2):
    return v1[X]*v2[Y] - v1[Y]*v2[X]

def vectoriel3D(v1,v2):
    return (
        v1[Y]*v2[Z] - v1[Z]*v2[Y],
        v1[Z]*v2[X] - v1[X]*v2[Z],
        v1[X]*v2[Y] - v1[Y]*v2[X]
    )
def compute(xxx_todo_changeme, xxx_todo_changeme1,da,db):
    (xa,ya) = xxx_todo_changeme
    (xb,yb) = xxx_todo_changeme1
    xa=float(xa)
    ya=float(ya)
    xb=float(xb)
    yb=float(yb)
    da=float(da)
    db=float(db)
    Dab = sqrt((xb-xa)**2 + (yb-ya)**2)
    xp=(da**2-db**2 + Dab**2) / (2.0*Dab)
    t=da**2-xp**2
    if t<0:
        return ()
    if t==0:
        yp1=0
        yp2=0
    else:
        yp1=sqrt(t)
        yp2=-yp1
    i = normed( ((xb-xa), (yb-ya)) )
    x1=xa + xp * i[0]  - yp1*i[1]
    y1=ya + xp * i[1]  + yp1*i[0]
    x2=xa + xp * i[0]  - yp2*i[1]
    y2=ya + xp * i[1]  + yp2*i[0]
    return ( (x1,y1) , (x2,y2) )

def middle(a,b):
    return (a[X]+(b[X]-a[X])/2,a[Y]+(b[Y]-a[Y])/2)


def distance_point_to_line(p1,p2,p):
    return math.fabs((p2[Y] - p1[Y]) * p[X] - (p2[X] - p1[X])* p[Y] + p2[X]*p1[Y] - p2[Y]*p1[X]) / math.sqrt((p2[Y]-p1[Y])**2 + (p2[X] - p1[X])**2)
    

def get_point_projection_on_plane(p1,p2,p3,m):
    #print p1,p2,p3,m
    normale=normed(
        vectoriel3D(
            p2v(p1,p2),
            p2v(p1,p3)
        )
    )
    p=plus(
        m,
        prod(
            normale,
            -scalaire(p2v(p1,m),normale)
        )
    )
    # check if it is in triangle
    w1 = ( p1[X]*(p3[Y]-p1[Y]) + (p[Y] - p1[Y])*(p3[X]-p1[X])-p[X]*(p3[Y]-p1[Y])) \
         / \
         ( (p2[Y] - p1[Y])*(p3[X]-p1[X]) - (p2[X]-p1[X])*(p3[Y]-p1[Y]) )
    w2 = ( p[Y] - p1[Y] - w1*(p2[Y]-p1[Y]) ) / (p3[Y]-p1[Y])

    d1 = distance_point_to_line(p1,p2,p)
    d2 = distance_point_to_line(p1,p3,p)
    d3 = distance_point_to_line(p2,p3,p)
    
    return p , w1 >=0 and w2 >=0 and w1+w2<=1.0 , min(d1,d2,d3)           
    
