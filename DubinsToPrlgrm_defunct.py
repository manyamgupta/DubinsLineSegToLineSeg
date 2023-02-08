# Author: Satyanarayana Gupta Manyam
# Shortest Dubins path to a paralellogram and given final heading
# All functions here assumes start config is [0,0,0]

import numpy as np
from shapely.geometry import Polygon
from shapely.geometry import Point
from collections import namedtuple
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
import utils
import dubins
import dubutils as du
import DubToLineSeg as dls

def DubinsToPrlgrm(prlGrm, finHdng, rho, pathType='Dubins'):
    # inputs:
    # prlgrm: list of vertices of the parallelogram in counterclock-wise order
    # finHdng: final heading at any position inside paralellogram
    # rho: minimum turn radius
    # Assumes start config is [0,0,0]

    minLength = np.nan
    minConfig = [np.nan, np.nan, np.nan]
    prlgrmPoly = Polygon(prlGrm)
    # prlGrm.append(prlGrm[0])
    lengthsVec = []
    configsList = []

    if pathType == 'S' or pathType =='Dubins':
        if prlgrmPoly.contains(Point(0,0)) and np.abs(finHdng)<=1e-8:
            lengthsVec.append(0)
            configsList.append(( (0,0,0), 'None', [0] ) )

    # Path Type L to any point with given final heading
    if pathType == 'L' or pathType =='Dubins':
        finPos, lengthL = dls.PathLtoFinHdng(finHdng, rho)
        if prlgrmPoly.contains(Point(finPos)):
            lengthsVec.append(lengthL)
            configsList.append(( (finPos[0], finPos[1], finHdng), 'L', [lengthL] ) )
    
    # Path Type R to any point with given final heading
    if pathType == 'R' or pathType =='Dubins':
        finPos, lengthR = dls.PathRtoFinHdng(finHdng, rho)    
        if prlgrmPoly.contains(Point(finPos)):
            lengthsVec.append(lengthR)
            configsList.append(( (finPos[0], finPos[1], finHdng), 'R', [lengthR] ) )

    for k in range(4):
                
        lineSeg = np.array([prlGrm[k], prlGrm[np.mod(k+1,4)]])

        # Path Type S to any point on paralellogram edges with given final heading
        if pathType == 'S' or pathType =='Dubins':
            finPos, lengthS = dls.PathStoLine(lineSeg, finHdng )
            if np.isfinite(finPos[0]):
                lengthsVec.append(lengthS)
                configsList.append(( (finPos[0], finPos[1], finHdng), 'S', [lengthS] ) )

        # Path Type LR to any point on paralellogram edges with given final heading
        if pathType == 'LR' or pathType =='Dubins':
            finPosList, lengthsList = dls.LRtoLine(lineSeg, finHdng, rho)   
            for j in range(len(finPosList)):
                finPos = finPosList[j]
                lengthsLR = lengthsList[j]
                lengthsVec.append(lengthsLR[0])
                configsList.append(( (finPos[0], finPos[1], finHdng), 'LR', [lengthsLR[1], lengthsLR[2]] ) )

        # Path Type RL to any point on paralellogram edges with given final heading
        if pathType == 'RL' or pathType =='Dubins':
            finPosList, lengthsList = dls.RLtoLine(lineSeg, finHdng, rho)   
            for j in range(len(finPosList)):
                finPos = finPosList[j]
                lengthsRL = lengthsList[j]
                lengthsVec.append(lengthsRL[0])
                configsList.append(( (finPos[0], finPos[1], finHdng), 'RL', [lengthsRL[1], lengthsRL[2]] ) )
    
        # Path Type LS to any point on paralellogram edges with given final heading
        if pathType == 'LS' or pathType =='Dubins':
            finPos, lengthsLS = dls.LStoLine(lineSeg, finHdng, rho) 
            if np.isfinite(finPos[0]):   
                lengthsVec.append(lengthsLS[0])
                configsList.append(( (finPos[0], finPos[1], finHdng), 'LS', [lengthsLS[1], lengthsLS[2]] ) )

        # Path Type RS to any point on paralellogram edges with given final heading
        if pathType == 'RS' or pathType =='Dubins':
            finPos, lengthsRS = dls.RStoLine(lineSeg, finHdng, rho) 
            if np.isfinite(finPos[0]):   
                lengthsVec.append(lengthsRS[0])
                configsList.append(( (finPos[0], finPos[1], finHdng), 'RS', [lengthsRS[1], lengthsRS[2]] ) )

        # Path Type SL to any point on paralellogram edges with given final heading
        if pathType == 'SL' or pathType =='Dubins':
            finPos, lengthsSL = dls.SLtoLine(lineSeg, finHdng, rho) 
            if np.isfinite(finPos[0]):   
                lengthsVec.append(lengthsSL[0])
                configsList.append(( (finPos[0], finPos[1], finHdng), 'SL', [lengthsSL[1], lengthsSL[2]] ) )

        # Path Type SR to any point on paralellogram edges with given final heading
        if pathType == 'SR' or pathType =='Dubins':
            finPos, lengthsSR = dls.SRtoLine(lineSeg, finHdng, rho) 
            if np.isfinite(finPos[0]):   
                lengthsVec.append(lengthsSR[0])
                configsList.append(( (finPos[0], finPos[1], finHdng), 'SR', [lengthsSR[1], lengthsSR[2]] ) )

        # Path Type LSL_min to any point on paralellogram edges with given final heading
        if pathType == 'LSL' or pathType =='Dubins':
            finPos, lengthLSL_min, pathLSL = dls.LocalMinLSL(lineSeg, finHdng, rho)
            if np.isfinite(finPos[0]):   
                lengthsVec.append(lengthLSL_min)
                configsList.append(( (finPos[0], finPos[1], finHdng), 'LSL', [pathLSL.segment_length(0), pathLSL.segment_length(1), pathLSL.segment_length(2)] ) )
                # du.PlotDubinsPath(pathLSL)

        # Path Type RSR_min to any point on paralellogram edges with given final heading
        if pathType == 'RSR' or pathType =='Dubins':
            finPos, lengthRSR_min, pathRSR = dls.LocalMinRSR(lineSeg, finHdng, rho)
            if np.isfinite(finPos[0]):   
                lengthsVec.append(lengthRSR_min)
                configsList.append(( (finPos[0], finPos[1], finHdng), 'RSR', [pathRSR.segment_length(0), pathRSR.segment_length(1), pathRSR.segment_length(2)] ) )

        # Path Type LSR_min to any point on paralellogram edges with given final heading
        if pathType == 'LSR' or pathType =='Dubins':
            finPos, lengthLSR_min, pathLSR = dls.LocalMinLSR(lineSeg, finHdng, rho)
            if np.isfinite(finPos[0]):   
                lengthsVec.append(lengthLSR_min)
                configsList.append(( (finPos[0], finPos[1], finHdng), 'LSR', [pathLSR.segment_length(0), pathLSR.segment_length(1), pathLSR.segment_length(2)] ) )

        # Path Type RSL_min to any point on paralellogram edges with given final heading
        if pathType == 'RSL' or pathType =='Dubins':
            finPos, lengthRSL_min, pathRSL = dls.LocalMinRSL(lineSeg, finHdng, rho)
            if np.isfinite(finPos[0]):   
                lengthsVec.append(lengthRSL_min)
                configsList.append(( (finPos[0], finPos[1], finHdng), 'RSL', [pathRSL.segment_length(0), pathRSL.segment_length(1), pathRSL.segment_length(2)] ) )

        # Path Type LRL_feas to any point on paralellogram edges with given final heading
        if pathType == 'LRL' or pathType =='Dubins':
            finPosList, lengthsList = dls.LRLFeasLimits(lineSeg, finHdng, rho)
            for j in range(len(finPosList)):
                finPos = finPosList[j]
                lengthsLRL = lengthsList[j]        
                lengthsVec.append(lengthsLRL[0])
                configsList.append(( (finPos[0], finPos[1], finHdng), 'LRL', lengthsLRL[1:4] ) )

        # Path Type RLR_feas to any point on paralellogram edges with given final heading
        if pathType == 'RLR' or pathType =='Dubins':
            finPosList, lengthsList = dls.RLRFeasLimits(lineSeg, finHdng, rho)
            for j in range(len(finPosList)):
                finPos = finPosList[j]
                lengthsRLR = lengthsList[j]        
                lengthsVec.append(lengthsRLR[0])
                configsList.append(( (finPos[0], finPos[1], finHdng), 'RLR', lengthsRLR[1:4] ) )

        # Shortest Dubins to vertices of the paralellogram
        if pathType =='Dubins':
            pathDub = dubins.shortest_path([0,0,0], [prlGrm[k][0], prlGrm[k][1], finHdng ], rho)
            lengthsVec.append(pathDub.path_length())
            configsList.append(( (prlGrm[k][0], prlGrm[k][1], finHdng ), du.PathNumtoType(pathDub.path_type()), [pathDub.segment_length(0), pathDub.segment_length(1), pathDub.segment_length(2)] ) )


        # Path of given type to vertices of the paralellogram
        if pathType in ['LSL','LSR','RSL','RSR','LRL','RLR']:
            pathDub = dubins.path([0,0,0], [prlGrm[k][0], prlGrm[k][1], finHdng ], rho, du.PathTypeToNum(pathType))
            if pathDub != None:
                lengthsVec.append(pathDub.path_length())
                configsList.append(( (prlGrm[k][0], prlGrm[k][1], finHdng ), du.PathNumtoType(pathDub.path_type()), [pathDub.segment_length(0), pathDub.segment_length(1), pathDub.segment_length(2)] ) )

    lengthsVec = np.array(lengthsVec)
    if lengthsVec.size >0:
        minLength = np.min(lengthsVec)
        minInd = np.argmin(lengthsVec)
        minConfig = configsList[minInd]         
    # utils.PlotParalellogram(prlGrm)
    # for i in range(len(lengthsVec)):
    #     du.PlotDubPathSegments([0,0,0], configsList[i][1], configsList[i][2], rho)
    # plt.axis('equal')


    return minLength, minConfig
    

if __name__ == "__main__":

    # prlGrm = [(-2,-1), (3,-1),(4,5),(-1,5)]
    # prlGrm = [(2,-1), (7,-1),(8,5),(3,5)]
    prlGrm =  [(1,-1), (4,-1), (2,3), (-1,3) ]

    finHdng =0
    rho=1
    plotformat = namedtuple("plotformat","color linewidth linestyle marker")

    start = time.time()
    minLength, minConfig = DubinsToPrlgrm(prlGrm, finHdng, rho)
    computation_time = time.time()-start

    print(f"{minLength=}")
    print(f"{minConfig=}")
    print(f"{computation_time=}")

    utils.PlotParalellogram(prlGrm, plotformat('g',2,'-','x'))
    utils.PlotArrow([0,0],0,1, plotformat('c',2,'--','x'))
    du.PlotDubPathSegments([0,0,0],minConfig[1],minConfig[2],rho, plotformat('b',2,'-',''))
    plt.axis('equal')
