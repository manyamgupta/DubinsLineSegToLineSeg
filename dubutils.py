# Author: Satyanarayana Gupta Manyam

import matplotlib.pyplot as plt
from numpy import pi,cos,sin
import numpy as np
import utils

def PlotDubinsPath(dubPath,fmt):
    qs, _ = dubPath.sample_many(.02)    
    qs = np.array(qs)
    xs = qs[:, 0]
    ys = qs[:, 1]
    plt.plot(xs, ys, color=fmt.color, linewidth=fmt.linewidth, linestyle=fmt.linestyle) 
    return qs

def PlotDubPathSegments(iniConf, pathMode, segLengths, rho, fmt):

    startConf = iniConf
    for k in range(len(pathMode)):
        startConf = PlotSegment(startConf, segLengths[k], pathMode[k], rho, fmt)

    return

def PlotSegment(startConf, segLength, segType, rho, fmt):

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
    lengths = np.array([lenLR, phi1, phi2])

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