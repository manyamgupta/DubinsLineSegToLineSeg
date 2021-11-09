# Author: Satyanarayana Gupta Manyam

from logging import setLoggerClass
import matplotlib.pyplot as plt
from numpy import pi,cos,sin
import numpy as np
from shapely.ctypes_declarations import EXCEPTION_HANDLER_FUNCTYPE
from shapely.geometry import LineString
from shapely.geometry import Point
from shapely.geometry import MultiPoint
from shapely.geometry.polygon import Polygon
from collections import namedtuple
import utils
import time
import dubins
from types import SimpleNamespace

deffmt = SimpleNamespace(color='blue', linewidth=2, linestyle='-', marker='x')
def PlotDubinsPath(dubPath,fmt=deffmt):
    qs, _ = dubPath.sample_many(.02)    
    qs = np.array(qs)
    xs = qs[:, 0]
    ys = qs[:, 1]
    plt.plot(xs, ys, color=fmt.color, linewidth=fmt.linewidth, linestyle=fmt.linestyle) 
    return qs

def PlotDubPathSegments(iniConf, pathMode, segLengths, rho, fmt=deffmt):

    if pathMode != 'None':
        for k in range(len(pathMode)):
            iniConf = PlotSegment(iniConf, segLengths[k], pathMode[k], rho, fmt)

    return

def PlotArc(cntr, rho, angPosEnds, fmt=deffmt):
    # Arcs are defined by the center, radius and the angular positions of the ends
    # The input angPosEnds1/2 shoudl contain two values of the ends, and the arc always goes in ccw from first to second

    startPos = cntr+rho*np.array([np.cos(angPosEnds[0]), np.sin(angPosEnds[0])])
    startConf = [startPos[0], startPos[1], angPosEnds[0]+np.pi/2]
    segLength = np.mod(angPosEnds[1]-angPosEnds[0], 2*np.pi)*rho
    if segLength <= 1e-9: segLength = 2*np.pi*rho
    PlotSegment(startConf, segLength, 'L', rho, fmt)

    return
def PlotSegment(startConf, segLength, segType, rho, fmt=deffmt):

    pt1 = startConf[0:2]
    t1 = startConf[2]
           
    if segType == 'S':
        pt2 = startConf[0:2] + segLength*np.array([cos(t1), sin(t1)])
        plt.plot([pt1[0],pt2[0]], [pt1[1],pt2[1]], color=fmt.color, linewidth=fmt.linewidth, linestyle=fmt.linestyle)
        finalConf = np.array([pt2[0], pt2[1], t1])
    elif segType == 'L' or segType == 'R':

        rotSense = 1 if segType =='L' else -1
        center = pt1 + rho*np.array([cos(t1+rotSense*pi/2), sin(t1+rotSense*pi/2)])

        alVec = np.linspace(t1-rotSense*pi/2 , t1 -rotSense*pi/2 +rotSense*segLength/rho,100)

        tc_x = center[0]+rho*cos(alVec)
        tc_y = center[1]+rho*sin(alVec)
        plt.plot(tc_x, tc_y, color=fmt.color, linewidth=fmt.linewidth, linestyle=fmt.linestyle) 
        t2 = np.mod(t1  + rotSense*segLength/rho, 2*pi)
        finalConf = [tc_x[-1], tc_y[-1], t2]
    else:
        raise Exception('Error: Ineligible turn, each letter can be L or R or S. Let me let you go')


    return finalConf

def PathTypeNum(ptype):

    typesList = ['L', 'R', 'S', 'LR', 'RL', 'LS', 'RS', 'SL', 'SR', 'LSL', 'RSR', 'LSR', 'RSL', 'LRL', 'RLR']
    return typesList.index(ptype)+1

def PathTypeToNum(ptype):

    typesList = ['LSL', 'LSR', 'RSL', 'RSR', 'RLR', 'LRL']
    indx = typesList.index(ptype)
    return indx

def PathNumtoType(pnum):
    # LSL =0; LSR = 1; RSL = 2; RSR = 3; RLR = 4; LRL = 5; 

    typesList = ['LSL', 'LSR', 'RSL', 'RSR', 'RLR', 'LRL']

    return typesList[pnum]

def RotSenseSeq(pathMode):

    rotSeq = []
    for k in range(len(pathMode)):
        w = pathMode[k]
        if w =='L':
            rotSeq.append(1)
        elif w =='R':
            rotSeq.append(-1)
        elif w=='S':
            rotSeq.append(0)
        else:
            raise Exception('Error: Incorrect turn, each letter can be L or R or S')

    return np.array(rotSeq)

def DistPtToLineSeg(pt, lineSeg):
    # perpendicular distance from point to linesegment
    lineSeg  =np.array(lineSeg)    
    A = lineSeg[0]; B = lineSeg[1]
    lenAB = np.linalg.norm(A-B)
    triangleArea = abs((B[0]-A[0])*(A[1]-pt[1]) - (B[1]-A[1])*(A[0]-pt[0]))
    
    return triangleArea/lenAB

def DistPtHdngToLineSeg(pt, hdng, lineSeg):
    # distance from point to linesegment along a given heading
    lineSeg  =np.array(lineSeg)    
    A = lineSeg[0]; B = lineSeg[1]
    PB = B-pt; PA = A-pt
    if np.cross(PB,PA)<0:
        B = lineSeg[0]; A = lineSeg[1]

    h = DistPtToLineSeg(pt, lineSeg) #perpendicual distance
    hdng_perp = np.arctan2(B[1]-A[1], B[0]-A[0]) + np.pi/2
    theta = hdng-hdng_perp

    return h/np.cos(theta)

def LengthOfLR(iniConf, finConf, rho):
        
    C2 = finConf[0:2] + rho*np.array([sin(finConf[2]), -cos(finConf[2])])
    C1 = iniConf[0:2] + rho*np.array([-sin(iniConf[2]), cos(iniConf[2])])

    gamma = np.arctan2(C2[1]-C1[1], C2[0]-C1[0])
    phi1 = np.mod(pi/2 + gamma, 2*pi)
    phi2 = np.mod(phi1-finConf[2], 2*pi)
    lenLR = rho*(phi1+phi2)
    lengths = np.array([lenLR, rho*phi1, rho*phi2])

    return lengths

def IntersectionLineCircle(lineSeg, C, r):

    A = lineSeg[0]; B = lineSeg[1]
    A = np.array(A); B = np.array(B)
    delX = B[0]-A[0]
    delY = B[1]-A[1]
    quadEqnCoefs = np.array([delX**2+delY**2, 2*delX*(A[0]-C[0])+2*delY*(A[1]-C[1]), (A[0]-C[0])**2+(A[1]-C[1])**2-r**2 ])
    lams = np.roots(quadEqnCoefs)

    intPts = []
    for k in range(0,2):
        if np.imag(lams[k]) == 0 and lams[k]>=0 and lams[k]<=1:
            intPt = A + lams[k]*(B-A)
            intPts.append(intPt)

    return intPts


def IntersectionCircleCircle(cntr1, rho1, cntr2, rho2):
    # using the answer from here
    # https://math.stackexchange.com/questions/256100/how-can-i-find-the-points-at-which-two-circles-intersect


    cntr1 = np.array(cntr1)
    cntr2 = np.array(cntr2)
    
    intPts = []
    R = np.linalg.norm(cntr1-cntr2)
    if R>rho1+rho2:
        return intPts
    exp1 = 0.5*(cntr1+cntr2)
    exp2 = ((rho1**2-rho2**2)/R**2/2)*(cntr2-cntr1)
    exp3 = 0.5*( np.sqrt(2*((rho1**2+rho2**2)/R**2) -((rho1**2-rho2**2)**2)/R**4 -1  ) )*np.array([cntr2[1]-cntr1[1], cntr1[0]-cntr2[0]])
    intPts.append(exp1+exp2+exp3)
    intPts.append(exp1+exp2-exp3)

    return intPts

def IntersectionArcARC(cntr1, rho1, angPosEnds1, cntr2, rho2, angPosEnds2):
    # Arcs are defined by the center, radius and the angular positions of the ends
    # The input angPosEnds1/2 shoudl contain two values of the ends, and the arc always goes in ccw from first to second

    intPtsCircles = IntersectionCircleCircle(cntr1, rho1, cntr2, rho2)
    intPts_Arcs = []
    for intPtC in intPtsCircles:
        angPos_ints_arc1 = np.arctan2(intPtC[1]-cntr1[1], intPtC[0]-cntr1[0])
        angPos_ints_arc2 = np.arctan2(intPtC[1]-cntr2[1], intPtC[0]-cntr2[0])

        if utils.InInt(angPosEnds1[0], angPosEnds1[1], angPos_ints_arc1) and utils.InInt(angPosEnds2[0], angPosEnds2[1], angPos_ints_arc2):
            intPts_Arcs.append(intPtC)

    return intPts_Arcs

def LineReflectionXaxis(lineSeg):

    lineSegRefl  = np.array(lineSeg)
    lineSegRefl[0,1] = -lineSeg[0][1]
    lineSegRefl[1,1] = -lineSeg[1][1]
    return lineSegRefl

def PtReflectionXaxis(pt):

    return [pt[0],-pt[1]]

def CheckPtLiesOnLineSeg(pt, lineSeg):

    lineSeg = np.array(lineSeg)
    if abs(np.linalg.norm(lineSeg[0]-pt)+np.linalg.norm(lineSeg[1]-pt) - np.linalg.norm(lineSeg[1]-lineSeg[0])) < 1e-6:
        return True
    else:
        return False

def DubinsInflexionPoints(dubPath):
    
    phi1 = dubPath.segment_length(0)
    phi3 = dubPath.segment_length(2)
    pathLength = dubPath.path_length()
    
    p1 = dubPath.sample(phi1+.0001)
    p2 = dubPath.sample(pathLength-phi3-0.0001)
    
    p1 = np.array([p1[0], p1[1]])
    p2 = np.array([p2[0], p2[1]])
    return p1,p2

def PathDefiningPoints(startConf, pathType, segLengths, rho):

    defPoints = []   
    centerAngPos = [] # Each entry in this list contains the centre of arc and angular positin of the ends, nan for straight lines 
    startPt = startConf[0:2]
    t1 = startConf[2]
    defPoints.append(startPt)
    for k in range(len(pathType)):     
        segType = pathType[k]
        segLength = segLengths[k]
        
        if segType == 'S':
            startPt = startPt[0:2] + segLength*np.array([cos(t1), sin(t1)])
            centerAngPos.append([np.nan, np.nan, np.nan, np.nan])
            
        elif segType == 'L' or segType == 'R':

            rotSense = 1 if segType =='L' else -1
            center = startPt + rho*np.array([cos(t1+rotSense*pi/2), sin(t1+rotSense*pi/2)])
            al1 = t1 -rotSense*pi/2
            al2 = t1 -rotSense*pi/2 +rotSense*segLength/rho
            startPt = center + rho*np.array([np.cos(al2), np.sin(al2)])
            t1 = np.mod(t1  + rotSense*segLength/rho, 2*pi)
            centerAngPos.append([center[0], center[1], al1, al2])

        defPoints.append(startPt)

    return np.array(defPoints), centerAngPos

def CheckArcIntersectsLine(cntr, rho, angPosEnds, lineEnds):
#   Assumes the arc always turns left

    intsFlag = False
    intPts = IntersectionLineCircle(lineEnds, cntr, rho)
    if len(intPts)>0:
        for intPt in intPts:
            al = np.arctan2(intPt[1]-cntr[1], intPt[0]-cntr[0])
            if utils.InInt(angPosEnds[0], angPosEnds[1], al):
                return True
    
    return intsFlag



def DubPathTypeString(pathType):

    if pathType == 0:
        pathString = ['L','S','L']
    elif pathType == 1:
        pathString = ['L','S','R']
    elif pathType == 2:
        pathString = ['R','S','L']
    elif pathType == 3:
        pathString = ['R','S','R']
    elif pathType == 4:
        pathString = ['R','L','R']
    elif pathType == 5:
        pathString = ['L','R','L']
    else:
        raise Exception('Incorrect pathType')


    return pathString

def FindIntPts(startConf, pathType, segLengths, rho, listObs):
    # Assumes the obstacles are circles
    # each entry in the list is tuple or array, first two entries of array are centre of circle, and third is radius

    intPtsAll = []
    defPoints, centerAngPosList = PathDefiningPoints(startConf, pathType, segLengths, rho)
    for k in range(len(pathType)):
        segType = pathType[k]
        
        centerAngPos = centerAngPosList[k]
        if segType == 'S':
            lineEnds = [defPoints[k], defPoints[k+1]]
            for obs in listObs:
                intPts = IntersectionLineCircle(lineEnds, obs[0:2], obs[2])
                if len(intPts) > 0:
                    intPtsAll.extend(intPts)
        else:
            if segType == 'L':
                cntr = centerAngPos[0:2]
                al1 = centerAngPos[2]
                al2 = centerAngPos[3]
            elif segType == 'R':
                cntr = centerAngPos[0:2]
                al1 = centerAngPos[3]
                al2 = centerAngPos[2]
            for obs in listObs:
                ints_arcs = IntersectionArcARC(cntr, rho, [al1, al2], obs[0:2], obs[2], [0, 2*np.pi])
                if len(ints_arcs) > 0:
                    intPtsAll.extend(ints_arcs)          

    return intPtsAll

def CheckDubPathFeas(startConf, pathType, segLengths, rho, listObs):
    # Assumes the obstacles are circles
    # each entry in the list is tuple or array, first two entries of array are centre of circle, and third is radius

    pathFeas = True
    defPoints, centerAngPosList = PathDefiningPoints(startConf, pathType, segLengths, rho)
    for k in range(len(pathType)):
        segType = pathType[k]
        
        centerAngPos = centerAngPosList[k]
        if segType == 'S':
            lineEnds = [defPoints[k], defPoints[k+1]]
            for obs in listObs:
                intPts = IntersectionLineCircle(lineEnds, obs[0:2], obs[2])
                if len(intPts) > 0:
                    return False
        else:
            if segType == 'L':
                cntr = centerAngPos[0:2]
                al1 = centerAngPos[2]
                al2 = centerAngPos[3]
            elif segType == 'R':
                cntr = centerAngPos[0:2]
                al1 = centerAngPos[3]
                al2 = centerAngPos[2]
            for obs in listObs:
                ints_arcs = IntersectionArcARC(cntr, rho, [al1, al2], obs[0:2], obs[2], [0, 2*np.pi])
                if len(ints_arcs) > 0:
                    return False                

    return pathFeas

if __name__ == "__main__":

    plotformat = namedtuple("plotformat","color linewidth linestyle marker")
    rho = 1
    pathStr = 'LSR'
    pathSegLengths = [2, 3, 4]
    iniConfig = [0,0,0]

    finConfig = [4,4,.5]
    obsList = [ [(1,1),(5,2),(6,5),(2,5)] ]
    # dubPath = dubins.shortest_path(iniConfig, finConfig, rho)

    # pathFeas = IsDubPathFeas(dubPath, rho, obsList)
    # print(f"{pathFeas=}")

    # utils.PlotPolygon(obsList[0], np.array([[0, 1], [1, 2], [2, 3], [3, 0]]), plotformat('r',2,'-',''))
    # PlotDubinsPath(dubPath, plotformat('b',2,'-',''))
    # plt.axis('equal')
    # plt.show()

    # cntr1 = Point([2,3])
    # circleBoundary1 = cntr1.buffer(rho).boundary

    # cntr2 = Point([3,4])
    # circleBoundary2 = cntr2.buffer(rho).boundary
    # p1 = Polygon(obsList[0])
    # intersection = circleBoundary1.intersection(circleBoundary2)
    # print(intersection)

    # utils.PlotCircle(cntr1.coords[0], rho, plotformat('b',2,'-',''))
    # utils.PlotCircle(cntr2.coords[0], rho, plotformat('b',2,'-',''))
    # utils.PlotPolygon(obsList[0], np.array([[0, 1], [1, 2], [2, 3], [3, 0]]), plotformat('r',2,'-',''))
    # plt.axis('equal')
    # plt.show()

    # c1 = [2,1]; r1 = 2.5
    # c2 = [-2,-2]; r2 = 3.5
    # # cntr1 = Point(c1)
    # # cntr2 = Point(c2)
    # # circleBoundary1 = cntr1.buffer(r1).boundary
    # # circleBoundary2 = cntr2.buffer(r2).boundary
    # # tic = time.time()
    # # ints_shapely = circleBoundary1.intersection(circleBoundary2)
    # # compTimeShapely = time.time()-tic

    # tic = time.time()
    # angPosEnds1 = [0,2*np.pi]
    # angPosEnds2 = [0,2*np.pi]
    # ints_arcs = IntersectionArcARC(c1, r1, angPosEnds1, c2, r2, angPosEnds2)
    # comTimeMyFunc = time.time()-tic

    # # print(f"{comTimeMyFunc=}")

    # PlotArc(c1, r1, angPosEnds1, plotformat('b',2,'-',''))
    # PlotArc(c2, r2, angPosEnds2, plotformat('g',2,'-',''))

    # if ints_arcs:
    #     intPts = np.array(ints_arcs)
    #     plt.scatter(intPts[:,0],intPts[:,1], marker='x', color = 'r')
    # plt.axis('equal')
    # plt.show()

    # startConf = [2,3,.5]
    # pathType = 'LSR'
    # segLengths = [1, np.pi, 1]
    # rho=1.5
    # defPoints, centerAngPosList = PathDefiningPoints(startConf, pathType, segLengths, rho)
    # PlotDubPathSegments(startConf, pathType, segLengths, rho, plotformat('b',2,'-',''))
    # # plt.scatter(defPoints[:,0], defPoints[:,1], marker= 'x', color='g')
    # for cntrAngPos in centerAngPosList:
    #     pt = cntrAngPos[0:2]+rho*np.array([np.cos(cntrAngPos[2]), np.sin(cntrAngPos[2])])
    #     plt.scatter(pt[0], pt[1], marker= 'x', color='r')
    #     pt = cntrAngPos[0:2]+rho*np.array([np.cos(cntrAngPos[3]), np.sin(cntrAngPos[3])])
    #     plt.scatter(pt[0], pt[1], marker= 'x', color='r')

    # plt.axis('equal')
    # plt.show()

    startConf = [1,-2,.5]
    pathType = 'LSR'
    segLengths = [.4, 5.5, 9]
    rho=2
    listObs = []
    listObs.append([1,2,2])
    listObs.append([8,8,3])
    listObs.append([10,0,4])
    listObs.append([2,9,2])

    # pathFeas = CheckDubPathFeas(startConf, pathType, segLengths, rho, listObs)
    # print(f"{pathFeas=}")
    intPts = FindIntPts(startConf, pathType, segLengths, rho, listObs)
    
    PlotDubPathSegments(startConf, pathType, segLengths, rho, plotformat('b',2,'-',''))
    
    for obs in listObs:
       utils.PlotCircle(obs[0:2], obs[2], plotformat('r',2,'-',''))
    plt.axis('equal')
    plt.show()


