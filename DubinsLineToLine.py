# Author: Satyanarayana Gupta Manyam
# Shortest Dubins path from a line and sector to a line and sector

import numpy as np
from numpy.lib.shape_base import expand_dims
from shapely.geometry import Polygon
from shapely.geometry import Point
from collections import namedtuple
import matplotlib.pyplot as plt
import time
import dubins
import utils
import dubutils as du
import DubToLineSeg as dls
import DubinsToPrlgrm as dpgm

def Translate(refPt, pt2):
    # this tranlates teh refPt to origin, and the pt2 to the new position
    return np.array(pt2) - np.array(refPt)

def RotateVec(vec, theta ):
    # rotates the vec in ccw direction for an angle of theta
    rotMat = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    return np.matmul(rotMat, np.array(vec))

def TransformLinesToPrlGrm(line1, line2):
    # This function translates the line1 and line 2, 
    # such that the starting position is (0,0) 
    # and the final feasible set of positions is a paralellogram    
    # output: paralellogram
    # For any path from [0,0] to paralellogram, there exists a path from line1 to line2
    refPt = np.array(line1[0])

    line1_pt2_trans = np.array(line1[1])-refPt
    line2_pt1_trans = np.array(line2[0])-refPt
    line2_pt2_trans = np.array(line2[1])-refPt
    prlGrm = [line2_pt1_trans, line2_pt1_trans-line1_pt2_trans, line2_pt2_trans-line1_pt2_trans, line2_pt2_trans]

    return prlGrm

def LinesToPrlGrm(line1, line2):
    # This does not move the first point of line1 to origin
    # Output is a paralellogram such that any path from line1[0] to paralellogram is same as line1 to line2
    line1 = np.array(line1)
    line2 = np.array(line2)
    delXY_line1 = line1[1]-line1[0]
    line1_pt2 = np.array(line1[1])

    prlGrm = [line2[0], line2[0]-delXY_line1, line2[1]-delXY_line1, line2[1]]

    return prlGrm

def TransformRotate(angle1, line1, line2, sector2):
    # This function translates and rotate the line1 and line 2, 
    # such that the starting position is (0,0) and start heading is 0
    # and the final feasible set of positions is a paralellogram
    prlGrm_rot = []
    prlGrm = TransformLinesToPrlGrm(line1, line2)
    sector2 = [sector2[0]-angle1, sector2[1]-angle1]
    for k in range(4):
        vertex_rot= RotateVec(prlGrm[k], -angle1 )
        prlGrm_rot.append(vertex_rot)
    return prlGrm_rot, sector2

def FindMinPositions(finalPos, line1, line2):

    # inputs are original line1 and line2.
    # finalPos is the final position of the dubins path starting from line1[0]
    # The finalPos may not lie on line2, because that is the output of the point to pralellogram problem
    # This function moves the final position and initial position such that they lie on the lines

    tol = 1e-6
    line1 = np.array(line1)
    line2 = np.array(line2)
    delXY_line1 = line1[1]-line1[0]

    line2_proj = [line2[0]-delXY_line1, line2[1]-delXY_line1]
    posLine1 = [np.nan, np.nan]
    posLine2 = [np.nan, np.nan]
    # utils.PlotLineSeg(line2_proj[0], line2_proj[1], plotformat('g',2,'-',''))
    if du.CheckPtLiesOnLineSeg(finalPos, line2):
        posLine1 = line1[0]
        posLine2 = finalPos
        return posLine1, posLine2
    elif du.CheckPtLiesOnLineSeg(finalPos, line2_proj):
        posLine1 = line1[1]
        posLine2 = finalPos+delXY_line1
        return posLine1, posLine2
    else:
        for i in range(2):
            delXY = line2[i]-finalPos
            if (np.dot(delXY, delXY_line1)/np.linalg.norm(delXY)/np.linalg.norm(delXY_line1)>=1-tol and
            np.dot(delXY, delXY_line1)/np.linalg.norm(delXY)/np.linalg.norm(delXY_line1)<=1+tol):
                posLine2 = line2[i]
                posLine1 = line1[0]+delXY
                return posLine1, posLine2

    prlgrm = LinesToPrlGrm(line1, line2)
    prlgrm_poly = Polygon(prlgrm)

    # print(f"{finalPos=}")
    # utils.PlotLineSeg(line1[0], line1[1], plotformat('g',2,'-',''))
    # utils.PlotLineSeg(line2[0], line2[1], plotformat('g',2,'-',''))
    # utils.PlotParalellogram(prlgrm, plotformat('c',2,'--',''))
    # plt.axis('equal')
    # plt.show()

    if prlgrm_poly.contains(Point(finalPos)):
        line1_ornt = np.arctan2(line1[1][1]-line1[0][1], line1[1][0]-line1[0][0])
        dist_prime = du.DistPtHdngToLineSeg(finalPos, line1_ornt, line2)
        posLine1 = line1[0] + dist_prime*np.array([np.cos(line1_ornt), np.sin(line1_ornt)])
        posLine2 = finalPos + dist_prime*np.array([np.cos(line1_ornt), np.sin(line1_ornt)])
        return posLine1, posLine2

    return posLine1, posLine2

def DubinsLineToLine(line1, sector1, line2, sector2, rho):
    # inputs:
    # line1 and line2: list of two points, each point is a tuple, also can be np.array
    # sector 1 and sector 2: tuples/arrays, each contains sector lower bound and upper bound (in that order)
    # rho: minimum turn radius
    # outputs;
    # minLength: shortest path length
    # minConfStart, minConfGoal: configurations corresponding to shortest dubins path on line 1 and line 2
    # minPathType: type of the shortest path, ex 'LSL', 'SR' etc.
    # minPathSegLengths: length of each of the segments in the shortest path

    minLengthsList = []
    minConfigsList = [] # List of configs: tuple of 4 entries (sector1 heading, final configuration in the rotated system, type of Dubins path, list of lengths of each segment )
    
    minLength = 0
    minConfig = []

    for i in range(2):
        # in the two iterations, the paralellogam is rotated such that the lower boundary and upper boundary
        # of the sector 1 aligns with the x-axis (0 degree heading).

        prlGrm_rot, sector2_rot = TransformRotate(sector1[i], line1, line2, sector2)
        prlGrm_rot_poly = Polygon(prlGrm_rot)
        if prlGrm_rot_poly.contains(Point(0,0)) and utils.InInt(sector2_rot[0], sector2_rot[1], 0):
            minLengthsList.append(0)
            minConfigsList.append((sector1[i], (0,0,0), 'None', [0] ) )
        #find the minimum path from [0,0,0] to the paralellogram with heading as lower boundary of sector2
        minLength1, minConfig1 = dpgm.DubinsToPrlgrm(prlGrm_rot, sector2_rot[0], rho)
        minLength2, minConfig2 = dpgm.DubinsToPrlgrm(prlGrm_rot, sector2_rot[1], rho)
        minLengthsList.append(minLength1)
        minLengthsList.append(minLength2)
        minConfigsList.append( (sector1[i], minConfig1[0], minConfig1[1], minConfig1[2]) )
        minConfigsList.append( (sector1[i], minConfig2[0], minConfig2[1], minConfig2[2])  )

        for k in range(4):
            lineSeg = np.array([prlGrm_rot[k], prlGrm_rot[k+1]])
            #path L from [0,0,0] to the line segment, final heading not specified
            finPosList, finHdngList, lengthsList = dls.LtoLine(lineSeg, rho)
            for j in range(len(lengthsList)):
                if utils.InInt(sector2_rot[0], sector2_rot[1], finHdngList[j]):
                    minLengthsList.append(lengthsList[j])
                    minConfigsList.append( ( sector1[i], (finPosList[j][0], finPosList[j][1], finHdngList[j]), 'L', [lengthsList[j]] ) ) 

            #path R from [0,0,0] to the line segment, final heading not specified
            finPosList, finHdngList, lengthsList = dls.RtoLine(lineSeg, rho)
            for j in range(len(lengthsList)):
                if utils.InInt(sector2_rot[0], sector2_rot[1], finHdngList[j]):
                    minLengthsList.append(lengthsList[j])
                    minConfigsList.append( (sector1[i], (finPosList[j][0], finPosList[j][1], finHdngList[j]), 'R', [lengthsList[j]] ) )

            #LR local minimum
            finPt, finHdng, lengthsMin = dls.LocalMinLR(lineSeg, rho)
            if np.isfinite(finPt[0]):
                if utils.InInt(sector2_rot[0], sector2_rot[1], finHdng ):
                    minLengthsList.append(lengthsMin[0])
                    minConfigsList.append( (sector1[i], (finPt[0], finPt[1], finHdng), 'LR', lengthsMin[1:3] ) )

            #RL local minimum
            finPt, finHdng, lengthsMin = dls.LocalMinRL(lineSeg, rho)
            if np.isfinite(finPt[0]):
                if utils.InInt(sector2_rot[0], sector2_rot[1], finHdng ):
                    minLengthsList.append(lengthsMin[0])
                    minConfigsList.append( (sector1[i], (finPt[0], finPt[1], finHdng), 'RL', lengthsMin[1:3] ) )

            #LS local minimum
            finPt, finHdng, lengthsMin = dls.LocalMinLS(lineSeg, rho)
            if np.isfinite(finPt[0]):
                if utils.InInt(sector2_rot[0], sector2_rot[1], finHdng ):
                    minLengthsList.append(lengthsMin[0])
                    minConfigsList.append( (sector1[i], (finPt[0], finPt[1], finHdng), 'LS', lengthsMin[1:3] ) )

            #RS local minimum
            finPt, finHdng, lengthsMin = dls.LocalMinRS(lineSeg, rho)
            if np.isfinite(finPt[0]):
                if utils.InInt(sector2_rot[0], sector2_rot[1], finHdng ):
                    minLengthsList.append(lengthsMin[0])
                    minConfigsList.append( (sector1[i], (finPt[0], finPt[1], finHdng), 'RS', lengthsMin[1:3] ) )

            #SL local minimum
            finPt, finHdng, lengthsMin = dls.LocalMinSL(lineSeg, rho)
            if np.isfinite(finPt[0]):
                if utils.InInt(sector2_rot[0], sector2_rot[1], finHdng ):
                    minLengthsList.append(lengthsMin[0])
                    minConfigsList.append( (sector1[i], (finPt[0], finPt[1], finHdng), 'SL', lengthsMin[1:3] ) )

            #SR local minimum
            finPt, finHdng, lengthsMin = dls.LocalMinSR(lineSeg, rho)
            if np.isfinite(finPt[0]):
                if utils.InInt(sector2_rot[0], sector2_rot[1], finHdng ):
                    minLengthsList.append(lengthsMin[0])
                    minConfigsList.append( (sector1[i], (finPt[0], finPt[1], finHdng), 'SR', lengthsMin[1:3] ) )
            if i==0:
                # This does not depend on the rotation, and needs to be computed only once
                #S local minimum
                finPt, finHdng, lengthS = dls.LocalMinS(lineSeg, rho)
                if np.isfinite(lengthS):
                    if utils.InInt(sector2_rot[0], sector2_rot[1], finHdng) and utils.InInt(0, sector1[1]-sector1[0], finHdng):
                        minLengthsList.append(lengthS)
                        minConfigsList.append( (sector1[0], (finPt[0], finPt[1], finHdng), 'S', [lengthS] ) )

    minLengthsVec = np.array(minLengthsList)
    minLength = np.min(minLengthsVec)
    minInd = np.argmin(minLengthsVec)
    minConfig = minConfigsList[minInd]

    minPathType = minConfig[2]
    minPathSegLengths = minConfig[3]
    # Rotating the final configuration and final heading back to the original orientation
    finalPos_revRotd = RotateVec([minConfig[1][0], minConfig[1][1] ], minConfig[0] )
    finalHdng_revRotd = minConfig[1][2] + minConfig[0]
    finalPos_moved = finalPos_revRotd + line1[0]
        
    posLine1, posLine2 = FindMinPositions(finalPos_moved, line1, line2)
    if minConfig[2] == 'S':
        minConfStart = [posLine1[0], posLine1[1], finalHdng_revRotd]
        minConfGoal = [posLine2[0], posLine2[1], finalHdng_revRotd]
    else:
        minConfStart = [posLine1[0], posLine1[1], minConfig[0]]
        minConfGoal = [posLine2[0], posLine2[1], finalHdng_revRotd]

    
    # print(f"{finalPos_revRotd=}")
    # print(f"{finalHdng_revRotd=}")
    # print(f"{finalPos_moved=}")

    # prlGrm = TransformLinesToPrlGrm(line1, line2)
    # prlGrm_moved = [prlGrm[0]+line1[0], prlGrm[1]+line1[0], prlGrm[2]+line1[0], prlGrm[3]+line1[0]]
    # plt.figure()
    # utils.PlotLineSeg(line1[0], line1[1], plotformat('b',2,'--',''))
    # utils.PlotLineSeg(line2[0], line2[1], plotformat('b',2,'--',''))
    # utils.PlotParalellogram(prlGrm_moved, plotformat('g',2,'--','x'))
    # # utils.PlotParalellogram(prlGrm, plotformat('g',2,'--','x'))
    # utils.PlotArrow([0,0], 0, 1, plotformat('c',2,'--','x'))

    # utils.PlotParalellogram(prlGrm, plotformat('m',2,'--','x'))
    # # du.PlotDubPathSegments([line1[0][0], line1[0][1], minConfig[0]], minConfig[2], minConfig[3],rho, plotformat('b',2,'-',''))
    # # du.PlotDubPathSegments(minConfStart, minConfig[2], minConfig[3],rho, plotformat('b',2,'-',''))
    # # utils.PlotArrow(minConfStart[0:2], minConfStart[2], 1, plotformat('c',2,'--','x'))
    # # utils.PlotArrow(minConfGoal[0:2], minConfGoal[2], 1, plotformat('c',2,'--','x'))
    # plt.axis('equal')
    # plt.show()

    return minLength, minConfStart, minConfGoal, minPathType, minPathSegLengths
    
def DubLineToLineNum(line1, sector1, line2, sector2, rho):
    # Numerically discretize and compute line-sector to line-sector shortest Dubins path
    
    discs = 31
    line1 = np.array(line1)
    line2 = np.array(line2)
    lengthVec = np.zeros(discs**4)
    configsVec = np.zeros([discs**4,6])
    indx = 0
    for i in range(discs):
        for j in range(discs):
            for k in range(discs):
                for l in range(discs):

                    lam1 = i/(discs-1)
                    lam2 = j/(discs-1)
                    lam3 = k/(discs-1)
                    lam4 = l/(discs-1)

                    p1 = (1-lam1)*line1[0]+lam1*line1[1]
                    t1 = (1-lam2)*sector1[0]+lam2*sector1[1]
                    p2 = (1-lam3)*line2[0]+lam3*line2[1]
                    t2 = (1-lam4)*sector2[0]+lam4*sector2[1]
                    dubPath = dubins.shortest_path([p1[0],p1[1],t1], [p2[0],p2[1],t2], rho)
                    length = dubPath.path_length()
                    lengthVec[indx] = length
                    configsVec[indx, :] = [p1[0],p1[1],t1, p2[0],p2[1],t2]
                    indx += 1
    minLength = np.min(lengthVec)
    minIndx = np.argmin(lengthVec)
    minConfigs = configsVec[minIndx]
    return minLength, minConfigs

    
if __name__ == "__main__":
    plotformat = namedtuple("plotformat","color linewidth linestyle marker")

    line1 = [(0,0), (5,1)]
    line2 = [(2,2), (3,5)]
    sector1 = [.25, 1]
    sector2 = [1, 3]
    rho = 1

    start = time.time()
    minLength, minConfStart, minConfGoal, minPathType, minPathSegLengths = DubinsLineToLine(line1, sector1, line2, sector2, rho)
    computation_time = time.time()-start
    print(f"{minLength=}")
    print(f"{minConfStart=}")
    print(f"{minConfGoal=}")
    print(f"{minPathType=}")
    print(f"{computation_time=}")

    # start = time.time()
    # minLengthNum, minConfigsNum = DubLineToLineNum(line1, sector1, line2, sector2, rho)
    # computation_time_num = time.time()-start
    # print(f"{minLengthNum=}")
    # print(f"{minConfigsNum=}")
    # print(f"{computation_time_num=}")

    

    plt.figure()
    utils.PlotLineSeg(line1[0], line1[1], plotformat('g',2,'-',''))
    utils.PlotLineSeg(line2[0], line2[1], plotformat('g',2,'-',''))
    utils.PlotArrow(line1[0], sector1[0], 1, plotformat('c',2, '-','x'))
    utils.PlotArrow(line1[0], sector1[1], 1, plotformat('c',2,'-','x'))
    utils.PlotArrow(line2[0], sector2[0], 1, plotformat('c',2,'-','x'))
    utils.PlotArrow(line2[0], sector2[1], 1, plotformat('c',2,'-','x'))
    
    if minPathType != 'None':
        du.PlotDubPathSegments(minConfStart, minPathType, minPathSegLengths,rho, plotformat('b',2,'-',''))
        utils.PlotArrow(minConfStart[0:2], minConfStart[2], 1, plotformat('m',2,'dotted','x'))
        utils.PlotArrow(minConfGoal[0:2], minConfGoal[2], 1, plotformat('m',2,'dotted','x'))

    plt.axis('equal')

    plt.show()

    ############### testing the translation and rotation functions ###############
    # pt1 = np.array([-6,-3])
    # pt2 = np.array([-5,8])

    # pt2_trans = Translate(pt1, pt2)

    # pt2_rot = RotateVec(pt2, np.pi/2 )
    # plt.figure()
    # utils.PlotLineSeg([0,0], pt2, plotformat('b',2,'-',''))
    # utils.PlotLineSeg([0,0], pt2_rot, plotformat('g',2,'-',''))

    ############### testing the tfunction TransformLinesToPrlGrm and TransformRotate  ###############

    # prlGrm = TransformLinesToPrlGrm(line1, line2)
    # prlGrm_rot, sector2_rot = TransformRotate(sector1[0], line1, line2, sector2)

    # utils.PlotLineSeg(line1[0], line1[1], plotformat('b',2,'-',''))
    # utils.PlotLineSeg(line2[0], line2[1], plotformat('b',2,'-',''))
    # utils.PlotArrow(line1[0], sector1[0],1,plotformat('c',2,'--',''))
    # utils.PlotArrow(line2[0], sector2[0],1,plotformat('c',2,'--',''))
    # utils.PlotArrow(line2[0], sector2[1],1,plotformat('c',2,'--',''))

    # plt.axis('equal')
    # plt.figure()

    # utils.PlotParalellogram(prlGrm, plotformat('g',2,'--','x'))
    # utils.PlotParalellogram(prlGrm_rot, plotformat('g',2,'--','x'))
    # utils.PlotArrow(prlGrm_rot[0], sector2_rot[0],1,plotformat('c',2,'--',''))
    # utils.PlotArrow(prlGrm_rot[0], sector2_rot[1],1,plotformat('c',2,'--',''))
    # plt.scatter([0],[0],marker='x')
    # utils.PlotArrow([0,0], 0,1,plotformat('c',2,'--',''))
    # plt.axis('equal')
    # plt.show()
