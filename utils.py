from numpy import pi,cos,sin
import numpy as np
import matplotlib.pyplot as plt

def IntersectionLineSegments(p1,p2,p3,p4):
    # https://math.stackexchange.com/questions/3176543/intersection-point-of-2-lines-defined-by-2-points-each
    n = p2-p1
    m = p3-p4
    p = p3-p1

    D = n[0]*m[1]-n[1]*m[0]
    if np.abs(D) < 1e-6:
        return [np.nan, np.nan]
    Qx = m[1]*p[0]-m[0]*p[1]
    t = Qx/D
    
    intPt = p1+t*(p2-p1)

    return intPt

def IntersectionLineSegments2(line1, line2):
    # https://math.stackexchange.com/questions/3176543/intersection-point-of-2-lines-defined-by-2-points-each
    p1 = np.array(line1.point1)
    p2 = np.array(line1.point2)
    p3 = np.array(line2.point1)
    p4 = np.array(line2.point2)
    n = p2-p1
    m = p3-p4
    p = p3-p1

    D = n[0]*m[1]-n[1]*m[0]
    if np.abs(D) < 1e-6:
        return (None, None)
    Qx = m[1]*p[0]-m[0]*p[1]
    t = Qx/D
    
    intPt = p1+t*(p2-p1)
    if CheckPtLiesInsideLineSeg(p1, p2, intPt) and CheckPtLiesInsideLineSeg(p3, p4, intPt):
        return intPt
    else:
        return (None, None)

    
def CheckPtLiesInsideLineSeg(a,b,c):
# check if point c lies inside line segment ab, assume abc are collinear

    if np.dot(c-a,b-c)>0:
        return True
    else:
        return False
    
    
def RotateCoords(coords, theta):
    R = np.matrix([[cos(theta), -sin(theta)],[sin(theta), cos(theta)]])  
    coords = np.matrix([coords[0], coords[1]])
    corods = coords.T
    coordsRot = R*corods
    
    return [coordsRot[0,0], coordsRot[1,0]]

# def PlotCircle(C, r,col):
#     alVec = np.linspace(0,2*pi,1000)

#     tc_x = C[0]+r*cos(alVec)
#     tc_y = C[1]+r*sin(alVec)
#     plt.plot(tc_x, tc_y, col,zorder=0) 
    
def PlotArc(C, r, phis,col):

    alVec = np.linspace(phis[0],phis[1],100)

    tc_x = C[0]+r*cos(alVec)
    tc_y = C[1]+r*sin(alVec)
    plt.scatter(C[0], C[1], marker='x', color=col)
    plt.plot(tc_x, tc_y, col) 
    

def Angdiff(ti, tf ):
   
    ti = np.mod(ti, 2*pi)
    tf = np.mod(tf, 2*pi)
    
    if (InInt(ti, tf, 0)):
        diff = tf+(2*pi-ti)
    else:
        diff = tf-ti
        
    diff = np.abs(diff)        
    
    return diff
    
def InInt(lb, ub, t ):   
    # Checks if t is in the interval (lb, ub), interval goes ccw from lb to ub
    if t is None or np.isnan(t):
        # print('invalid value for t: '+str(t))
        return False
    lb = np.mod(lb, 2*pi)
    ub = np.mod(ub, 2*pi)
    t = np.mod(t, 2*pi) 
    if lb==ub:
        print("improper interval, lb and ub cannot be the same")
        return False
    elif (lb > ub):
        if (t >= lb or t <= ub): 
            return True
        elif np.abs(t-lb)<1e-6 or np.abs(t-ub)<1e-6:
            return True
        else:
            return False
    else:
        if (t >= lb and t <= ub):
            return True
        elif np.abs(t-lb)<1e-6 or np.abs(t-ub)<1e-6:
            return True
        else:
            return False
            
def MidAng(lb, ub ):  
    # finds the middle angle between lb and up, going ccw from lb to ub
    lb = np.mod(lb, 2*pi)
    ub = np.mod(ub, 2*pi)
    
    if lb == 0:
        case =1
    elif ub ==0:
        case =2
    elif InInt(lb, ub, 0 ):
        case =2
    else:
        case =1
    
    if case ==1:
        midang = (lb+ub)/2
        midang = np.mod(midang, 2*pi)
        return midang
    elif case==2:
        midang = (lb-2*pi+ub)/2
        midang = np.mod(midang, 2*pi)
        return midang  
def CheckIntrvsIntersect(intr1, intr2):
    
    if InInt(intr1[0], intr1[1], intr2[0] ) or InInt(intr1[0], intr1[1], intr2[1] ):
        return True
    elif InInt(intr2[0], intr2[1], intr1[0] ) or InInt(intr2[0], intr2[1], intr1[1] ):
        return True        
    else:
        return False
          
def CommonAngle(intr1, intr2):
    # Assumes the interval intersects, and returns one common angle
    if InInt(intr1[0], intr1[1], intr2[0] ):
        return intr2[0]
    elif InInt(intr1[0], intr1[1], intr2[1] ):
        return intr2[1]
    elif InInt(intr2[0], intr2[1], intr1[0] ):
        return intr1[0]     
    elif  InInt(intr2[0], intr2[1], intr1[1] ):
        return intr1[1]
    else:
        return None
    
def RotateVec(vec, theta ):
    # rotates the vec in ccw direction for an angle of theta
    # returns tuple
    return (np.cos(theta)*vec[0]-np.sin(theta)*vec[1], np.sin(theta)*vec[0]+np.cos(theta)*vec[1])
            	
def PlotPolygon(vertices, segments, col):

    vertices = np.array(vertices)
    segments = np.array(segments)
    numSegs = np.size(segments, 0)

    for i in range(numSegs):
        seg = segments[i,:]
        vertsSeg = vertices[seg,:]
        plt.plot(vertsSeg[:,0], vertsSeg[:,1], color=col, linewidth=2)

    return

def PlotCircle(C, r,fmt):
    alVec = np.linspace(0,2*pi,1000)

    tc_x = C[0]+r*cos(alVec)
    tc_y = C[1]+r*sin(alVec)
    plt.plot(tc_x, tc_y, color=fmt.color, linewidth=fmt.linewidth, linestyle=fmt.linestyle,zorder=0) 
    plt.scatter(C[0], C[1],marker='x')

    return

def PlotArc(C, r, phis,fmt):

    alVec = np.linspace(phis[0],phis[1],100)

    tc_x = C[0]+r*cos(alVec)
    tc_y = C[1]+r*sin(alVec)
    plt.scatter(C[0], C[1], color=fmt.color, marker=fmt.marker)
    plt.plot(tc_x, tc_y, color=fmt.color, linewidth=fmt.linewidth, linestyle=fmt.linestyle) 

    return

def PlotArrow(p1, hdng, arrLen, fmt):

    # p2 = p1 + arrLen*np.array([cos(hdng), sin(hdng)])
    dx = arrLen*np.cos(hdng)
    dy = arrLen*np.sin(hdng)

    # plt.plot([p1[0],p2[0]], [p1[1],p2[1]], color=fmt.color, linewidth=fmt.linewidth, linestyle=fmt.linestyle)     
    plt.arrow(p1[0],p1[1],dx,dy,head_width=.1*np.sqrt(dx**2+dy**2),color=fmt.color,linewidth=fmt.linewidth, linestyle=fmt.linestyle)
    return

def PlotInterval(pos, intr, fmt):
    PlotArrow(pos, intr[0], fmt.arrLen, fmt)
    PlotArrow(pos, intr[1], fmt.arrLen, fmt)
    
    return
    
def PlotLineSeg(p1, p2, fmt):

    plt.plot([p1[0],p2[0]], [p1[1],p2[1]], color=fmt.color, linewidth=fmt.linewidth, linestyle=fmt.linestyle) 

    return
