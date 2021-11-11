# Author: Satyanarayana Gupta Manyam
# Shortest Dubins path from a line and sector to a line and sector

from types import SimpleNamespace
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
        vertex_rot= utils.RotateVec(prlGrm[k], -angle1 )
        prlGrm_rot.append(vertex_rot)
    return prlGrm_rot, sector2

def FindMinPositions(finalPos, line1, line2):

    # inputs are original line1 and line2.
    # finalPos is the final position of the dubins path starting from line1[0]
    # The finalPos may not lie on line2, because that is the output of the point to pralellogram problem
    # This function moves the final position and initial position such that they lie on the original lines

    tol = 1e-6
    line1 = np.array(line1)
    line2 = np.array(line2)
    delXY_line1 = line1[1]-line1[0]

    line2_proj = [line2[0]-delXY_line1, line2[1]-delXY_line1]
    posLine1 = [np.nan, np.nan]
    posLine2 = [np.nan, np.nan]
    # utils.PlotLineSeg(line2_proj[0], line2_proj[1], plotformat('g',2,'-',''))
    if du.CheckPtLiesOnLineSeg(finalPos, line2): 
        #if the final pos is on original line2, then line1[0] is the start position
        posLine1 = line1[0]
        posLine2 = finalPos
        return posLine1, posLine2
    elif du.CheckPtLiesOnLineSeg(finalPos, line2_proj):
        #if the final pos is on projected line2, then line1[1] is the start position
        posLine1 = line1[1]
        posLine2 = finalPos+delXY_line1
        return posLine1, posLine2
    else:
        # this is the case when the final pos lies on a line parallel to line 1 and passing through an end point of line2
        for i in range(2):
            delXY = line2[i]-finalPos
            # if (np.dot(delXY, delXY_line1)/np.linalg.norm(delXY)/np.linalg.norm(delXY_line1)>=1-tol and
            # np.dot(delXY, delXY_line1)/np.linalg.norm(delXY)/np.linalg.norm(delXY_line1)<=1+tol):
            if np.abs( np.dot(delXY, delXY_line1)/np.linalg.norm(delXY)/np.linalg.norm(delXY_line1)-1)<=tol:        
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

def PathOneSeg(prlgrm, sector2, rho, pathType):

    lengthsList = []
    configsList = []
    minLength = np.nan
    minConfig = [np.nan, np.nan, np.nan]
    dtpLength, dtpConfig = dpgm.DubinsToPrlgrm(prlgrm, sector2[0], rho, pathType)
    if np.isfinite(dtpLength): 
        lengthsList.append(dtpLength)
        configsList.append( (dtpConfig[0], dtpConfig[1], dtpConfig[2]) )

    dtpLength, dtpConfig = dpgm.DubinsToPrlgrm(prlgrm, sector2[1], rho, pathType)
    if np.isfinite(dtpLength): 
        lengthsList.append(dtpLength)
        configsList.append( (dtpConfig[0], dtpConfig[1], dtpConfig[2]) )

    if pathType != 'S':
        for k in range(4):
            lineSeg = np.array([prlgrm[k], prlgrm[np.mod(k+1,4)]])
            #path L from [0,0,0] to the line segment, final heading not specified
            if pathType == 'L':
                finPosList, finHdngList, seglengths = dls.LtoLine(lineSeg, rho)
            elif pathType == 'R':
                finPosList, finHdngList, seglengths = dls.RtoLine(lineSeg, rho)

            for j in range(len(seglengths)):
                if utils.InInt(sector2[0], sector2[1], finHdngList[j]):
                    lengthsList.append(seglengths[j])
                    configsList.append( ((finPosList[j][0], finPosList[j][1], finHdngList[j]), pathType, [seglengths[j]] ) ) 

    if lengthsList:
        lengthsVec = np.array(lengthsList)
        minLength = np.min(lengthsVec)
        minInd = np.argmin(lengthsVec)
        minConfig = configsList[minInd]
    
    return minLength, minConfig



def PathTwoSegs(prlgrm, sector2, rho, pathType):

    lengthsList = []
    configsList = []
    minLength = np.nan
    minConfig = [np.nan, np.nan, np.nan]
    dtpLength, dtpConfig = dpgm.DubinsToPrlgrm(prlgrm, sector2[0], rho, pathType)
    if np.isfinite(dtpLength): 
        lengthsList.append(dtpLength)
        configsList.append( (dtpConfig[0], dtpConfig[1], dtpConfig[2]) )

    dtpLength, dtpConfig = dpgm.DubinsToPrlgrm(prlgrm, sector2[1], rho, pathType)
    if np.isfinite(dtpLength): 
        lengthsList.append(dtpLength)
        configsList.append( (dtpConfig[0], dtpConfig[1], dtpConfig[2]) )

    #local minimum
    for k in range(4):
        lineSeg = np.array([prlgrm[k], prlgrm[np.mod(k+1,4)]])
        if pathType == 'LS':
            finPt, finHdng, lengthsMin = dls.LocalMinLS(lineSeg, rho)
        elif pathType == 'LR':
            finPt, finHdng, lengthsMin = dls.LocalMinLR(lineSeg, rho)
        elif pathType == 'RS':
            finPt, finHdng, lengthsMin = dls.LocalMinRS(lineSeg, rho)
        elif pathType == 'RL':
            finPt, finHdng, lengthsMin = dls.LocalMinRL(lineSeg, rho)
        if np.isfinite(finPt[0]):
            if utils.InInt(sector2[0], sector2[1], finHdng ):
                lengthsList.append(lengthsMin[0])
                configsList.append( ((finPt[0], finPt[1], finHdng), pathType, lengthsMin[1:3] ) )
    if lengthsList:
        lengthsVec = np.array(lengthsList)
        minLength = np.min(lengthsVec)
        minInd = np.argmin(lengthsVec)
        minConfig = configsList[minInd]

    return minLength, minConfig

def DubinsLineToLineV2(line1, sector1, line2, sector2, rho):
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

    ############################################################################################################
    ######### The paralellogam is rotated such that the 0-degrees is aligned with sector 1 upper_limit #########
    
    prlGrm_rot, sector2_rot = TransformRotate(sector1[1], line1, line2, sector2) 

    prlGrm_rot_poly = Polygon(prlGrm_rot)
    if prlGrm_rot_poly.contains(Point(0,0)) and utils.InInt(sector2_rot[0], sector2_rot[1], 0):
        minLengthsList.append(0)
        minConfigsList.append((sector1[1], (0,0,0), 'None', [0] ) )

    #minimum LSR path from [0,0,0] to the paralellogram with final heading as upper boundary of sector2
    minLength, minConfig = dpgm.DubinsToPrlgrm(prlGrm_rot, sector2_rot[1], rho, 'LSR')
    minLengthsList.append(minLength)
    minConfigsList.append( (sector1[1], minConfig[0], minConfig[1], minConfig[2]) )

    for pType in ['LSL', 'LRL']:
        #minimum LSL or LRL path from [0,0,0] to the paralellogram with final heading as lower boundary of sector2
        minLength, minConfig = dpgm.DubinsToPrlgrm(prlGrm_rot, sector2_rot[0], rho, pType)
        if np.isfinite(minLength):
            minLengthsList.append(minLength)    
            minConfigsList.append( (sector1[1], minConfig[0], minConfig[1], minConfig[2])  )

    for pType in ['LS', 'LR', 'RL']:
        #minimum two segment path from [0,0,0] to the paralellogram with final heading in sector2
        minLength, minConfig = PathTwoSegs(prlGrm_rot, sector2_rot, rho, pType)
        if np.isfinite(minLength):
            minLengthsList.append(minLength)    
            minConfigsList.append( (sector1[1], minConfig[0], minConfig[1], minConfig[2])  )

    for pType in ['L', 'R', 'S']:
    #minimum one segment path from [0,0,0] to the paralellogram with final heading in sector2
        minLength, minConfig = PathOneSeg(prlGrm_rot, sector2_rot, rho, pType)
        if np.isfinite(minLength):
            minLengthsList.append(minLength)    
            minConfigsList.append( (sector1[1], minConfig[0], minConfig[1], minConfig[2])  )


    ############################################################################################################
    ######### The paralellogam is rotated such that the 0-degrees is aligned with sector 1 lower_limit #########
    prlGrm_rot, sector2_rot = TransformRotate(sector1[0], line1, line2, sector2) 
    
    prlGrm_rot_poly = Polygon(prlGrm_rot)
    if prlGrm_rot_poly.contains(Point(0,0)) and utils.InInt(sector2_rot[0], sector2_rot[1], 0):
        minLengthsList.append(0)
        minConfigsList.append((sector1[1], (0,0,0), 'None', [0] ) )

    #minimum RSL path from [0,0,0] to the paralellogram with final heading as lower boundary of sector2
    minLength, minConfig = dpgm.DubinsToPrlgrm(prlGrm_rot, sector2_rot[0], rho, 'RSL')
    minLengthsList.append(minLength)
    minConfigsList.append( (sector1[0], minConfig[0], minConfig[1], minConfig[2]) )

    for pType in ['RSR', 'RLR']:
        # #minimum RSR or RLR path from [0,0,0] to the paralellogram with final heading as upper boundary of sector2
        minLength, minConfig = dpgm.DubinsToPrlgrm(prlGrm_rot, sector2_rot[1], rho, pType)
        if np.isfinite(minLength):
            minLengthsList.append(minLength)
            minConfigsList.append( (sector1[0], minConfig[0], minConfig[1], minConfig[2]) )

    for pType in ['RS', 'RL', 'LR']:
    #minimum two segment path from [0,0,0] to the paralellogram with final heading in sector2
        minLength, minConfig = PathTwoSegs(prlGrm_rot, sector2_rot, rho,pType)
        if np.isfinite(minLength):
            minLengthsList.append(minLength)    
            minConfigsList.append( (sector1[0], minConfig[0], minConfig[1], minConfig[2])  )

    for pType in ['L', 'R', 'S']:
    #minimum one seg path from [0,0,0] to the paralellogram with final heading in sector2
        minLength, minConfig = PathOneSeg(prlGrm_rot, sector2_rot, rho, pType)
        if np.isfinite(minLength):
            minLengthsList.append(minLength)    
            minConfigsList.append( (sector1[0], minConfig[0], minConfig[1], minConfig[2])  )

    #S local minimum
    for k in range(4):
        lineSeg = np.array([prlGrm_rot[k], prlGrm_rot[np.mod(k+1,4)]])    
        finPt, finHdng, lengthS = dls.LocalMinS(lineSeg, rho)
        if np.isfinite(lengthS):
            if utils.InInt(sector2_rot[0], sector2_rot[1], finHdng) and utils.InInt(0, sector1[1]-sector1[0], finHdng):
                minLengthsList.append(lengthS)
                minConfigsList.append( (sector1[0], (finPt[0], finPt[1], finHdng), 'S', [lengthS] ) )

    # plt.figure()
    # prlgrm = LinesToPrlGrm(line1, line2)
    # fmt = SimpleNamespace(color='m',linewidth=2, linestyle='--')
    # utils.PlotParalellogram(prlgrm,fmt)
    # plt.axis('equal')
    # for k in range(len(minLengthsList)):
    #     du.PlotDubPathSegments([line1[0][0], line1[0][1], minConfigsList[k][0]], minConfigsList[k][2], minConfigsList[k][3], rho)
        # print(str(minConfigsList[k][2])+": "+ str(minLengthsList[k])+ ": " + str(minConfigsList[k][3]))

    ############################################################################################################
    ######### The paths with free initial heading and fixed final angle
    # ####### We find that by reversing line 1 and line2 #########

    line1_bkwd = line2
    line2_bkwd = line1
    sector1_bkwd = [sector2[0]+np.pi, sector2[1]+np.pi]
    sector2_bkwd = [sector1[0]+np.pi, sector1[1]+np.pi]

    minLengthsList_bkwd = []
    minConfigsList_bkwd = []

    ######### The paralellogam is rotated such that the 0-degrees is aligned with sector 1 upper_limit #########
    prlGrm_bkwd_rot, sector2_bkwd_rot = TransformRotate(sector1_bkwd[1], line1_bkwd, line2_bkwd, sector2_bkwd) 
    
    for pType in ['LS', 'LR']:
        #minimum backward LS or LR path from [0,0,0] to the paralellogram with final heading in sector2
        minLength_bkwd, minConfig_bkwd = PathTwoSegs(prlGrm_bkwd_rot, sector2_bkwd_rot, rho, pType)
        if np.isfinite(minLength_bkwd):
            minLengthsList_bkwd.append(minLength_bkwd)    
            minConfigsList_bkwd.append( (sector1_bkwd[1], minConfig_bkwd[0], minConfig_bkwd[1], minConfig_bkwd[2])  )

    ######### The paralellogam is rotated such that the 0-degrees is aligned with sector 1 lower_limit #########
    prlGrm_bkwd_rot, sector2_bkwd_rot = TransformRotate(sector1_bkwd[0], line1_bkwd, line2_bkwd, sector2_bkwd) 

    for pType in ['RS', 'RL']:
        #minimum backward Rx path from [0,0,0] to the paralellogram with final heading in sector2
        minLength_bkwd, minConfig_bkwd = PathTwoSegs(prlGrm_bkwd_rot, sector2_bkwd_rot, rho, pType)
        if np.isfinite(minLength_bkwd):
            minLengthsList_bkwd.append(minLength_bkwd)    
            minConfigsList_bkwd.append( (sector1_bkwd[0], minConfig_bkwd[0], minConfig_bkwd[1], minConfig_bkwd[2])  )

    # print("Backward paths:")
    # plt.figure()
    # prlgrm_bkwd = LinesToPrlGrm(line1_bkwd, line2_bkwd)
    # fmt = SimpleNamespace(color='m',linewidth=2, linestyle='--')
    # utils.PlotParalellogram(prlgrm_bkwd,fmt)
    # for k in range(len(minLengthsList_bkwd)):
    #     du.PlotDubPathSegments([line1_bkwd[0][0], line1_bkwd[0][1], minConfigsList_bkwd[k][0]], minConfigsList_bkwd[k][2], minConfigsList_bkwd[k][3], rho)
        # print(str(minConfigsList_bkwd[k][2])+": "+ str(minLengthsList_bkwd[k]) + ": " + str(minConfigsList_bkwd[k][3]) )

    
    minLengthsVec = np.array(minLengthsList)
    minLength = np.min(minLengthsVec)
    minInd = np.argmin(minLengthsVec)
    minConfig = minConfigsList[minInd]
    minPathType = minConfig[2]
    minPathSegLengths = minConfig[3]
    minConfStart, minConfGoal = FindMinConfsSG(minConfig, line1, line2)
    if len(minLengthsList_bkwd)>0:
        if min(minLengthsList_bkwd)< minLength:
            minLengthsVec_bkwd = np.array(minLengthsList_bkwd)
            minLength = np.min(minLengthsVec_bkwd)
            minInd = np.argmin(minLengthsVec_bkwd)
            minConfig_bkwd = minConfigsList_bkwd[minInd]
            minPathType_bkwd = minConfig_bkwd[2]
            minPathSegLengths_bkwd = minConfig_bkwd[3]     
            minConfStart_bkwd, minConfGoal_bkwd = FindMinConfsSG(minConfig_bkwd, line1_bkwd, line2_bkwd)     

            minConfStart = [minConfGoal_bkwd[0], minConfGoal_bkwd[1], np.mod(minConfGoal_bkwd[2]+np.pi, 2*np.pi)]
            minConfGoal = [minConfStart_bkwd[0], minConfStart_bkwd[1], np.mod(minConfStart_bkwd[2]+np.pi, 2*np.pi)]
            minPathType = BackwardPathType(minPathType_bkwd)
            minPathSegLengths = np.flip(minPathSegLengths_bkwd)
   
    return minLength, minConfStart, minConfGoal, minPathType, minPathSegLengths

def BackwardPathType(pType):

    pType_bkwd = ""
    for seg in pType:
        if seg == 'L':
            pType_bkwd += 'R'
        elif seg == 'R':
            pType_bkwd += 'L'
        else:
            pType_bkwd += seg

    return pType_bkwd[::-1]

def FindMinConfsSG(minConfig, line1, line2):
    # Rotating the final configuration and final heading back to the original orientation
    finalPos_revRotd = utils.RotateVec([minConfig[1][0], minConfig[1][1] ], minConfig[0] )
    finalHdng_revRotd = minConfig[1][2] + minConfig[0]
    finalPos_moved = finalPos_revRotd + line1[0]
        
    posLine1, posLine2 = FindMinPositions(finalPos_moved, line1, line2)
    if minConfig[2] == 'S':
        minConfStart = [posLine1[0], posLine1[1], finalHdng_revRotd]
        minConfGoal = [posLine2[0], posLine2[1], finalHdng_revRotd]
    else:
        minConfStart = [posLine1[0], posLine1[1], minConfig[0]]
        minConfGoal = [posLine2[0], posLine2[1], finalHdng_revRotd]

    return minConfStart, minConfGoal
    
if __name__ == "__main__":
    plotformat = namedtuple("plotformat","color linewidth linestyle marker")

    line1 = [(-1,-1), (6,4)]
    line2 = [(2,2), (4, 6)]
    sector1 = [1, 2.5]
    sector2 = [5, 6]
    rho = 2

    start = time.time()
    minLength, minConfStart, minConfGoal, minPathType, minPathSegLengths = DubinsLineToLineV2(line1, sector1, line2, sector2, rho)
    computation_time = time.time()-start
    print(f"{minLength=}")
    print(f"{minConfStart=}")
    print(f"{minConfGoal=}")
    print(f"{minPathType=}")
    print(f"{computation_time=}")


    plt.figure()
    utils.PlotLineSeg(line1[0], line1[1], plotformat('g',2,'-',''))
    utils.PlotLineSeg(line2[0], line2[1], plotformat('g',2,'-',''))
    
    if minPathType != 'None':
        du.PlotDubPathSegments(minConfStart, minPathType, minPathSegLengths,rho, plotformat('b',2,'-',''))
        utils.PlotArrow(minConfStart[0:2], sector1[0], 1, plotformat('c',2, '-','x'))
        utils.PlotArrow(minConfStart[0:2], sector1[1], 1, plotformat('c',2,'-','x'))
        utils.PlotArrow(minConfGoal[0:2], sector2[0], 1, plotformat('c',2,'-','x'))
        utils.PlotArrow(minConfGoal[0:2], sector2[1], 1, plotformat('c',2,'-','x'))
        utils.PlotArrow(minConfStart[0:2], minConfStart[2], 1, plotformat('m',2,'dotted','x'))
        utils.PlotArrow(minConfGoal[0:2], minConfGoal[2], 1, plotformat('m',2,'dotted','x'))

    plt.axis('equal')

    plt.show()

    ############### testing the translation and rotation functions ###############
    # pt1 = np.array([-6,-3])
    # pt2 = np.array([-5,8])

    # pt2_trans = Translate(pt1, pt2)

    # pt2_rot = utils.RotateVec(pt2, np.pi/2 )
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
