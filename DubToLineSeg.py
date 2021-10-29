# Author: Satyanarayana Gupta Manyam
# Shortest Dubins path to a line segment
# All functions here assumes start config is [0,0,0]
# Local minimums are with respect to the lamda, parameter defining position on the line-segment

import numpy as np
from numpy import isfinite, pi,cos,sin
import matplotlib.pyplot as plt
import dubins 
import dubutils as du
import utils
from collections import namedtuple 

def PathLtoFinHdng(finHdng, rho):

    lengthL = np.mod(finHdng, 2*np.pi)*rho
    finPos = np.array([rho*sin(finHdng),rho-rho*cos(finHdng)])
    return finPos, lengthL 

def PathRtoFinHdng(finHdng, rho):

    lengthR = np.mod(-finHdng, 2*np.pi)*rho
    finPos = np.array([-rho*sin(finHdng),-rho+rho*cos(finHdng)])
    return finPos, lengthR

def LocalMinS(lineSeg, rho):
    A = lineSeg[0]; B = lineSeg[1]
    finPt = np.array([np.nan, np.nan])
    finHdng = np.nan
    lengthS = np.nan
    
    vx = B[0]-A[0]
    vy = B[1]-A[1]
    lam = -(A[0]*vx+A[1]*vy)/(vx*vx+vy*vy)
    if lam>=0 and lam <=1:
        finPt = (1-lam)*lineSeg[0]+ lam*lineSeg[1]
        finHdng = np.mod(np.arctan2(finPt[1], finPt[0]), 2*np.pi)
        lengthS = np.linalg.norm(finPt)

    return finPt, finHdng, lengthS

def PathStoLine(lineSeg, finHdng ):
    finPos = [np.nan, np.nan]
    lengthS = np.nan
    if np.abs(finHdng) < 1e-6:
        A = lineSeg[0]; B = lineSeg[1]
        if A[1]*B[1]<=0:
            lam = A[1]/(A[1]-B[1])
            if lam>=0 and lam<=1:
                lengthS = (1-lam)*A[0]+lam*B[0]
                if lengthS >= 0:
                    finPos = [lengthS,0]
                else: 
                    lengthS = np.nan
           
    return finPos, lengthS

def LtoLine(lineSeg, rho):    
    
    lineSeg = np.array(lineSeg)
    center = np.array([0, rho])
    h = du.DistPtToLineSeg(center, lineSeg)
    finPosList = []
    finHdngList = []
    finLengthsList = []
    if h<=rho:
        finPosList = du.IntersectionLineCircle(lineSeg, [0,rho], rho)
        
        for finPos in finPosList:            
            # if du.CheckPtLiesOnLineSeg(finPos, lineSeg):
            finHdng = np.arctan2(finPos[1]-rho, finPos[0])+pi/2
            finHdngList.append(finHdng)
            finLengthsList.append(rho*np.mod(finHdng, 2*np.pi))
                
    return finPosList, finHdngList, finLengthsList

def RtoLine(lineSeg, rho):

    # assumes initial configuration is [0,0,0]
    lineSeg_rfl = du.LineReflectionXaxis(lineSeg)
    finPosList_rfl, finHdngList_rfl, finLengthsList = LtoLine(lineSeg_rfl, rho)
    finPosList = []
    finHdngList = []
    for k in range(len(finPosList_rfl)):
        finPosList.append(du.PtReflectionXaxis(finPosList_rfl[k]))
        finHdngList.append( np.mod(-finHdngList_rfl[k], 2*np.pi ) )
        
    return finPosList, finHdngList, finLengthsList

def LRtoLine(lineSeg, finHdng, rho):
    # Computes path from (0,0,0) to any point on the lineSef (AB), with a given final heading (finHdng)
    
    A = lineSeg[0]; B = lineSeg[1]
    A = np.array(A); B = np.array(B)
    h = du.DistPtToLineSeg(np.array([0,rho]), lineSeg)
    if h>3*rho:
        finPtsList = []
        lengthsList = []
    else:
        delX = B[0]-A[0]
        delY = B[1]-A[1]
        T1 = A[0]+rho*sin(finHdng)
        T2 = A[1]-rho*cos(finHdng)-rho
        quadEqnCoeffs = np.array([delX**2 + delY**2, 2*T1*delX+2*T2*delY, T1**2+T2**2-4*rho*rho])
        lams = np.roots(quadEqnCoeffs)

        finPtsList = []
        lengthsList = []
        for k in range(0,2):
            if np.imag(lams[k]) == 0 and lams[k]>=0 and lams[k]<=1:
                finPt = A + lams[k]*(B-A)
                # if du.CheckPtLiesOnLineSeg(finPt, lineSeg):
                finPtsList.append(finPt)
                lengths = du.LengthOfLR(np.array([0,0,0]), np.array([finPt[0], finPt[1], finHdng]), rho)
                lengthsList.append(lengths)
        
    return finPtsList, lengthsList

def RLtoLine(lineSeg, finHdng, rho):

    lineSeg_rfl = du.LineReflectionXaxis(lineSeg)
    finHdng_rfl = np.mod(-finHdng, 2*np.pi)
    finPosList_rfl, lengthsList = LRtoLine(lineSeg_rfl, finHdng_rfl, rho)
    
    finPosList=[]
    for finPos in finPosList_rfl:
        finPosList.append(du.PtReflectionXaxis(finPos))
    
    return finPosList, lengthsList 

def LocalMinLR(lineSeg, rho):
    # finds the local minimum of the LR pat
    # return value lengthsMinLR contains LR path length, segment 1 and segment 2 lenght (in that order)
    A = lineSeg[0]; B = lineSeg[1]
    finPt = np.array([np.nan, np.nan])
    finHdng = np.nan
    lengthsMinLR = np.array([np.nan, np.nan, np.nan])
    h = du.DistPtToLineSeg(np.array([0,rho]), lineSeg)

    if h <= 3*rho:
        alpha = np.arccos(h/(3*rho))
        crossprodCACB = np.cross(A-np.array([0,rho]), B-np.array([0,rho]))
        if crossprodCACB < 0:
            finHdng = np.arctan2((B[1]-A[1]),(B[0]-A[0])) + alpha
        else:
            finHdng = np.arctan2((A[1]-B[1]),(A[0]-B[0])) + alpha
        finPtsList, lengthsList = LRtoLine(lineSeg, finHdng, rho)
        for k in range(len(lengthsList)):
            if lengthsList[k][2]/rho < pi-.0000000001:
                lengthsMinLR = lengthsList[k]
                finPt = finPtsList[k]
                break

    return finPt, finHdng, lengthsMinLR

def LocalMinRL(lineSeg, rho):

    lineSeg_rfl = du.LineReflectionXaxis(lineSeg)
    finPt, finHdng, lengthLocalMinRL = LocalMinLR(lineSeg_rfl, rho)

    if np.isfinite(finPt[0]):
        finPt = du.PtReflectionXaxis(finPt)
        finHdng = np.mod(-finHdng,np.pi)

    return finPt, finHdng, lengthLocalMinRL 

def LStoLine(lineSeg, finHdng, rho):
    A = lineSeg[0]; B = lineSeg[1]
    finHdng = np.mod(finHdng, 2*np.pi)
    C = np.array([0,rho])
    crossprodCACB = np.cross(A-C, B-C)
    if crossprodCACB > 0:
        A = lineSeg[1]; B = lineSeg[0]
    lengthsLS = [np.nan, np.nan, np.nan]
    finPos = [np.nan, np.nan]
    # angleCB = np.arctan2(B[1]-C[1], B[0]-C[0])
    # angleCA = np.arctan2(A[1]-C[1], A[0]-C[0])
    # alpha = np.arcsin(rho/np.linalg.norm(C-A))
    # beta = np.arcsin(rho/np.linalg.norm(C-B))
    # finHdngRange = (angleCB+beta, angleCA+alpha)
    lam = (A[1]*cos(finHdng)-A[0]*sin(finHdng)+rho-rho*cos(finHdng))/( (B[0]-A[0])*sin(finHdng)-(B[1]-A[1])*cos(finHdng) )
    # if utils.InInt(finHdngRange[0], finHdngRange[1], finHdng):
    if lam>=0 and lam<=1:
        finPos = (1-lam)*A + lam*B
        C2 = finPos + rho*np.array([-sin(finHdng), cos(finHdng)])
        ls = np.linalg.norm(C-C2)
        if np.abs((A[1]+lam*(B[1]-A[1])+rho*cos(finHdng)-rho)/ls - sin(finHdng) )<1e-6 and np.abs((A[0]+lam*(B[0]-A[0])-rho*sin(finHdng) )/ls - cos(finHdng))<1e-6:
            lengthsLS = [rho*finHdng+ls, rho*finHdng, ls]
        else:
            finPos = [np.nan, np.nan]
    return finPos, lengthsLS

def RStoLine(lineSeg, finHdng, rho):

    lineSeg_rfl = du.LineReflectionXaxis(lineSeg)
    finPos, lengthsRS = LStoLine(lineSeg_rfl, -finHdng , rho)
    finPos = du.PtReflectionXaxis(finPos)
    return finPos, lengthsRS

def LocalMinLS(lineSeg, rho):
    A = lineSeg[0]; B = lineSeg[1]
    C = np.array([0,rho])
    crossprodCACB = np.cross(A-C, B-C)
    if crossprodCACB > 0:
        A = lineSeg[1]; B = lineSeg[0]
    gamma = np.arctan2(B[1]-A[1], B[0]-A[0])
    finHdng = np.mod(gamma+np.pi/2, 2*np.pi)
    finPos, minLenLS = LStoLine(lineSeg, finHdng, rho)
    return finPos, finHdng, minLenLS

def LocalMinRS(lineSeg, rho):

    lineSeg_rfl = du.LineReflectionXaxis(lineSeg)
    finPos, finHdng, minLenRS = LocalMinLS(lineSeg_rfl, rho)
    
    if np.isfinite(finPos[0]):
        finPos = du.PtReflectionXaxis(finPos)
        finHdng = np.mod(-finHdng, 2*np.pi)

    return finPos, finHdng, minLenRS

def SLtoLine(lineSeg, finHdng, rho):
    A = lineSeg[0]; B = lineSeg[1]
    finHdng = np.mod(finHdng, 2*np.pi)
   
    crossprodOAOB = np.cross(A, B)
    if crossprodOAOB > 0:
        A = lineSeg[1]; B = lineSeg[0]
    lengthsSL = [np.nan, np.nan, np.nan]
    finPos = [np.nan, np.nan]
    if (B[1]-A[1]) != 0:
        lam = (rho-rho*cos(finHdng)-A[1])/(B[1]-A[1])
        ls = A[0]+lam*(B[0]-A[0])-rho*sin(finHdng)
        if lam>=0 and lam <=1 and ls>0:
            finPos = (1-lam)*A + lam*B        
            lengthsSL = [ls+rho*finHdng, ls, rho*finHdng]
        
    else:
        if A[1]>0 and A[1]< 2*rho:
            print("To Do: SL to horizontal line")

    return finPos, lengthsSL

def SRtoLine(lineSeg, finHdng, rho):
    lineSeg_rfl = du.LineReflectionXaxis(lineSeg)
    finPos, lengthsSR = SLtoLine(lineSeg_rfl, -finHdng , rho)
    finPos = du.PtReflectionXaxis(finPos)

    return finPos, lengthsSR

def LocalMinSL(lineSeg, rho):
    A = lineSeg[0]; B = lineSeg[1]
    crossprodCACB = np.cross(A, B)
    if crossprodCACB > 0:
        A = lineSeg[1]; B = lineSeg[0]
    gamma = np.arctan2(B[1]-A[1], B[0]-A[0])
    finHdng = np.pi + 2*gamma
    finPos, minlengthsSL = SLtoLine(lineSeg, finHdng, rho)
    return finPos, finHdng, minlengthsSL

def LocalMinSL2(lineSeg, finHdng, rho):
# Local minimum for given final heading, and the initial heading can vary
#  This occurs when 'S' segment is perpendicular to linesegment
    A = lineSeg[0]; B = lineSeg[1]
    crossprodCACB = np.cross(A, B)
    if crossprodCACB > 0:
        A = lineSeg[1]; B = lineSeg[0]
    gamma = np.arctan2(B[1]-A[1], B[0]-A[0])

    Ar = utils.RotateVec(A, -gamma+3*np.pi/2)
    Br = utils.RotateVec(B, -gamma+3*np.pi/2)
    finHdng_rot = finHdng-(gamma-3*np.pi/2)
    finPos_rot, minlengthsSL = SLtoLine([Ar, Br], finHdng_rot, rho)

    finPos = utils.RotateVec(finPos_rot, gamma-3*np.pi/2)
    iniHdng = np.mod(gamma-3*np.pi/2, 2*np.pi)

    return finPos, iniHdng, minlengthsSL

def LocalMinSR(lineSeg, rho):

    lineSeg_rfl = du.LineReflectionXaxis(lineSeg)
    finPos, finHdng, minLenSR = LocalMinSL(lineSeg_rfl, rho)
    
    if np.isfinite(finPos[0]):
        finPos = du.PtReflectionXaxis(finPos)
        finHdng = np.mod(-finHdng, 2*np.pi)
    return finPos, finHdng, minLenSR

def LocalMinSR2(lineSeg, finHdng, rho):

    lineSeg_rfl = du.LineReflectionXaxis(lineSeg)
    finHdng_rfl = -finHdng
    finPos, iniHdng, minLenSR = LocalMinSL2(lineSeg_rfl, finHdng_rfl, rho)
    if np.isfinite(finPos[0]):
        finPos = du.PtReflectionXaxis(finPos)
        iniHdng = np.mod(-iniHdng, 2*np.pi)
    return finPos, iniHdng, minLenSR

def LocalMinLSL(lineSeg, finHdng, rho):

    A = lineSeg[0]; B = lineSeg[1]
    delX = B[0]-A[0]
    delY = B[1]-A[1]

    lam = (-(A[1]+rho*cos(finHdng)-rho)*delY-(A[0]-rho*sin(finHdng))*delX)/(delX**2+delY**2)
    if lam>=0 and lam<=1:
        finPt = A+lam*(B-A)
        pathLSL = dubins.path([0,0,0],[finPt[0], finPt[1],finHdng],rho, 0)
        lengthLSL = pathLSL.path_length()
        if np.abs(pathLSL.segment_length(0)-2*np.pi) <1e-8:
            lengthLSL = lengthLSL-2*pi
        if np.abs(pathLSL.segment_length(2)-2*np.pi) <1e-8:
            lengthLSL = lengthLSL-2*pi
    else:
        return [np.nan, np.nan], np.nan, None

    return finPt, lengthLSL, pathLSL

def LocalMinRSR(lineSeg, finHdng, rho):
    lineSeg_rfl = du.LineReflectionXaxis(lineSeg)
    finHdng_rfl = np.mod(-finHdng, 2*np.pi)
    finPt, lengthRSR, pathRSR = LocalMinLSL(lineSeg_rfl, finHdng_rfl, rho)
    finPt = du.PtReflectionXaxis(finPt)
    if isfinite(finPt[0]):
        pathRSR = dubins.path([0,0,0],[finPt[0], finPt[1],finHdng],rho, 3)
        lengthRSR = pathRSR.path_length()

    return finPt, lengthRSR, pathRSR

def LocalMinLSR(lineSeg, finHdng, rho):

    A = lineSeg[0]; B = lineSeg[1]
    C1 = np.array([0,rho])
    delX = B[0]-A[0]
    delY = B[1]-A[1]

    h = du.DistPtToLineSeg(C1, lineSeg)
    crossprodCACB = np.cross(A-C1, B-C1)
    if crossprodCACB<0:
        al = finHdng-np.arctan2(delY, delX)
    else:
        al = finHdng-np.arctan2(-delY, -delX)
    # print('alpha: ', al)

    t1 = A[0]+rho*sin(finHdng)
    t2 = A[1]-rho*cos(finHdng)-rho
    t3 = (h-rho*cos(al))**2 + 4*(rho**2)
    quadEqCoeffs = np.array([delX**2+delY**2, 2*delX*t1+2*delY*t2, t1**2 + t2**2-t3])
    lamRoots = np.roots(quadEqCoeffs)
    
    finPt =np.array([np.nan, np.nan])
    lengthLSR = np.nan
    pathLSR = None

    for lam in lamRoots:
        if np.imag(lam) == 0 and lam>=0 and lam<=1:
            finPos = A+lam*(B-A)
            pathLSR_cand = dubins.path([0,0,0],[finPos[0], finPos[1],finHdng],rho, 1)
            p1,p2 = du.DubinsInflexionPoints(pathLSR_cand)
            if np.abs(np.dot(p2-p1, B-A)) < 1e-6:
                finPt = finPos
                pathLSR = pathLSR_cand
                lengthLSR = pathLSR.path_length()
                if np.abs(pathLSR.segment_length(0)-2*np.pi) <1e-8:
                    lengthLSR = lengthLSR-2*pi
                if np.abs(pathLSR.segment_length(2)-2*np.pi) <1e-8:
                    lengthLSR = lengthLSR-2*pi

    return finPt, lengthLSR, pathLSR

def LocalMinRSL(lineSeg, finHdng, rho):
    lineSeg_rfl = du.LineReflectionXaxis(lineSeg)
    finHdng_rfl = np.mod(-finHdng, 2*np.pi)
    finPt, lengthRSL, pathRSL = LocalMinLSR(lineSeg_rfl, finHdng_rfl, rho)
    finPt = du.PtReflectionXaxis(finPt)
    if isfinite(finPt[0]):
        pathRSL = dubins.path([0,0,0],[finPt[0], finPt[1],finHdng],rho, 2)
        lengthRSL = pathRSL.path_length()
        if np.abs(pathRSL.segment_length(0)-2*np.pi) <1e-8:
            lengthRSL = lengthRSL-2*pi
        if np.abs(pathRSL.segment_length(2)-2*np.pi) <1e-8:
                lengthRSL = lengthRSL-2*pi

    return finPt, lengthRSL, pathRSL

def LRLFeasLimits(lineSeg, finHdng, rho):

    A = lineSeg[0]; B = lineSeg[1]
    delX = B[0]-A[0]
    delY = B[1]-A[1]
    
    t1 = A[0]-rho*sin(finHdng)
    t2 = A[1]+rho*cos(finHdng)-rho
    
    quadEqCoeffs = np.array([delX**2+delY**2, 2*delX*t1 + 2*delY*t2, t1**2 + t2**2-16*rho*rho])
    lamRoots = np.roots(quadEqCoeffs)
    finPtsList = []
    lengthsList = []
    for k in range(0,2):
        if np.imag(lamRoots[k]) == 0 and lamRoots[k]>=0 and lamRoots[k]<=1:
            finPt = A + lamRoots[k]*(B-A)
            phi1 = np.pi/2 + np.arctan2(A[1]+lamRoots[k]*delY+rho*cos(finHdng)-rho, A[0]+lamRoots[k]*delX-rho*sin(finHdng))
            phi3 = finHdng+np.pi-phi1
            finPtsList.append(finPt)
            lengths = [rho*np.mod(phi1, 2*pi), rho*np.pi, rho*np.mod(phi3, 2*pi)]
            pathlength = sum(lengths)
            lengthsList.append([pathlength, lengths[0], lengths[1], lengths[2]])
    return finPtsList, lengthsList

def RLRFeasLimits(lineSeg, finHdng, rho):
    lineSeg_rfl = du.LineReflectionXaxis(lineSeg)
    finHdng_rfl = np.mod(-finHdng, 2*np.pi)

    finPtsList_rfl, lengthsList = LRLFeasLimits(lineSeg_rfl, finHdng_rfl, rho)
    finPtsList = []
    for finPt in finPtsList_rfl:
        finPtsList.append(du.PtReflectionXaxis(finPt) )

    return finPtsList, lengthsList



if __name__ == "__main__":

    # LSL =0; LSR = 1; RSL = 2; RSR = 3; RLR = 4; LRL = 5; 
    iniConf = np.array([0,0,0])
    rho = 1
    finHdng = 1
    plotformat = namedtuple("plotformat","color linewidth linestyle marker")

    ############################# Test L to line #############################
    # C1 = np.array([0,rho])  
    # lineSeg = np.array([ (-1,2), (4,1.9)])
    # finPosList, finHdngList, finLengthsList = LtoLine(lineSeg, rho)
    # utils.PlotCircle(C1, rho,plotformat('b',2,'--',''))
    # utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'--',''))
    # for k in range(len(finPosList)):
    #     finPos = finPosList[k]
    #     finHdng = finHdngList[k]
    #     if np.isfinite(finPos[0]):
    #         du.PlotDubPathSegments(iniConf, 'L', [finLengthsList[k]], rho, plotformat('b',2,'-',''))
    #         utils.PlotArrow(finPos, finHdng, 1, plotformat('m',2,'-',''))
    #         print('Length L: ', finLengthsList[k])
    # plt.axis('equal')
    # plt.show()

    ############################# Test L to given heading #############################
    # C1 = np.array([0,-rho]) 
    # finHdng = 4
    # finPos, length = PathLtoFinHdng(finHdng, rho)
    # utils.PlotCircle(C1, rho,plotformat('b',2,'--',''))
    # du.PlotDubPathSegments(iniConf, 'L', [length], rho, plotformat('b',2,'-',''))
    # utils.PlotArrow(finPos, finHdng, 1, plotformat('m',2,'-',''))
    # print(f"{length=}")
    # plt.axis('equal')
    # plt.show()
    ############################# Test R to line #############################
    # lineSeg = np.array([ (-1,-2), (3,-1)])
    # C1 = np.array([0,-rho])
    # lineSeg_rfl = du.LineReflectionXaxis(lineSeg)
    # finPosList, finHdngList, finLengthsList = RtoLine(lineSeg, rho)
    # utils.PlotCircle(C1, rho,plotformat('b',2,'--',''))
    # utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'--',''))
    # for k in range(len(finPosList)):
    #     finPos = finPosList[k]
    #     finHdng = finHdngList[k]
    #     if np.isfinite(finPos[0]):
    #         du.PlotDubPathSegments(iniConf, 'R', [finLengthsList[k]], rho, plotformat('b',2,'-',''))
    #         utils.PlotArrow(finPos, finHdng, 1, plotformat('m',2,'-',''))
    #         print('Length R: ', finLengthsList[k])
    #         print('Final Hdng: ', finHdng)
    # plt.axis('equal')
    # plt.show()

    # ############################# Test LR to LineSeg #############################
    # lineSeg = np.array([ (-2.5,0), (-1,3)])
    # finHdng = 2
    # finPtsList, lengthsList = LRtoLine(lineSeg, finHdng, rho)
    
    # C1 = np.array([0,rho])
    
    # utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'--',''))
    # utils.PlotCircle(np.array([0,rho]), rho,plotformat('g',2,'--',''))

    # for k in range(len(lengthsList)):
    #     lengths = lengthsList[k]
    #     finPt = finPtsList[k]
    #     center2 = finPt + rho*np.array([sin(finHdng), -cos(finHdng)])
    #     du.PlotDubPathSegments(iniConf, 'LR', lengths[1:3], rho, plotformat('b',2,'-',''))
    #     utils.PlotCircle(center2, rho,plotformat('g',2,'--',''))
    #     utils.PlotArrow(finPt, finHdng, 1, plotformat('m',2,'--',''))
        
    # plt.axis('equal')
    # plt.show()

    ############################# Test LR Min #############################
    lineSeg = np.array([ (2,-1), (-2,3)])
    
    lineSeg = np.array([ (0.5,-2), (1,5)] )

    finPt, finHdng, lengths = LocalMinLR(lineSeg, rho)
    
    C1 = np.array([0,rho])
    
    if np.isfinite(lengths[0]):
        utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'--',''))
        utils.PlotCircle(np.array([0,rho]), rho,plotformat('g',2,'--',''))
        
        du.PlotDubPathSegments(iniConf, 'LR', lengths[1:3], rho, plotformat('b',2,'-',''))
        
        # center2 = finPt + rho*np.array([sin(finHdng), -cos(finHdng)])
        C2 = finPt + rho*np.array([cos(finHdng-pi/2), sin(finHdng-pi/2)])
        utils.PlotCircle(C2, rho,plotformat('g',2,'--',''))
        utils.PlotArrow(finPt, finHdng, 1, plotformat('m',2,'--',''))
        
        P1 = (C1+C2)*.5
        utils.PlotLineSeg(P1, finPt, plotformat('g',2,'-','x'))
        plt.axis('equal')
        plt.show()

    ############################# Test S Min #############################
    
    # # lineSeg = np.array([(2,1), (1,-2)] )
    # lineSeg = np.array([[-1.1480503 ,  2.7716386 ], [ 5.45042262, -5.31910643]])

    # finPt, finHdng, lengthL = LocalMinS(lineSeg, rho)
    # print(f"{finPt=}")
    # print(f"{finHdng=}")
    
    # C1 = np.array([0,rho])
    
    # if np.isfinite(lengthL):
    #     utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'--',''))
    #     utils.PlotLineSeg([0,0], finPt, plotformat('b',2,'-',''))
        
    #     utils.PlotArrow(finPt, finHdng, 1, plotformat('m',2,'--',''))

    #     plt.axis('equal')
    #     plt.show()

     ############################ Test LS Min #############################
    
    # lineSeg = np.array([(-2,3), (-3,-2)])
    # finPos, finHdng, minLenLS = LocalMinLS(lineSeg, rho)
    
    # C1 = np.array([0,rho])
    
    # utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'--',''))
    # utils.PlotCircle(C1, rho,plotformat('g',2,'--',''))

    # if np.isfinite(finPos[0]):
    #     du.PlotDubPathSegments(iniConf, 'LS', minLenLS[1:3], rho, plotformat('b',2,'-',''))
    #     utils.PlotArrow(finPos, finHdng, 1, plotformat('m',2,'--',''))

    # plt.axis('equal')
    # plt.show()  
    
     ############################ Test RS Min #############################
    # lineSeg = np.array([(2,-3), (3,2)])
    # finPos, finHdng, minLenRS = LocalMinRS(lineSeg, rho)
    
    # C1 = np.array([0,-rho])
    
    # utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'--',''))
    # utils.PlotCircle(C1, rho,plotformat('g',2,'--',''))

    # if np.isfinite(finPos[0]):
    #     du.PlotDubPathSegments(iniConf, 'RS', minLenRS[1:3], rho, plotformat('b',2,'-',''))
    #     utils.PlotArrow(finPos, finHdng, 1, plotformat('m',2,'--',''))

    # plt.axis('equal')
    # plt.show()   

     ############################ Test SL Min #############################
    
    # lineSeg = np.array([ (-1,2), (4,1.9)])

    # finPos, finHdng, minLenSL = LocalMinSL(lineSeg, rho)
    
    # utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'-',''))
    # if np.isfinite(finPos[0]):
    #     du.PlotDubPathSegments(iniConf, 'SL', minLenSL[1:3], rho, plotformat('b',2,'-',''))
    #     p1 = [minLenSL[1],0]
    #     utils.PlotLineSeg(p1, finPos, plotformat('c',2,'--',''))
    #     utils.PlotArrow(finPos, finHdng, 1, plotformat('m',2,'--',''))
    #     print('Min Length SL: ',minLenSL[0])
    #     print(f"{finHdng=}")

    # plt.axis('equal')
    # plt.show()  

    # lineSeg = np.array([(4,-4), (5,1)])
    # finHdng = 1
    # finPos, iniHdng, minLenSL = LocalMinSL2(lineSeg, finHdng, rho)

    # utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'-',''))
    # if np.isfinite(finPos[0]):
    #     du.PlotDubPathSegments([0,0,iniHdng], 'SL', minLenSL[1:3], rho, plotformat('b',2,'-',''))
    #     utils.PlotArrow(finPos, finHdng, 1, plotformat('m',2,'--',''))
    #     print('Min Length SL: ',minLenSL[0])
    #     print(f"{iniHdng=}")

    # plt.axis('equal')
    # plt.show()  

    ############################ Test SR Min #############################
    
    # lineSeg = np.array([ (-1,-2.5), (3,-1.5)])
    # finPos, finHdng, minLenSR = LocalMinSR(lineSeg, rho)
    
    # utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'-',''))
    # if np.isfinite(finPos[0]):
    #     du.PlotDubPathSegments(iniConf, 'SR', minLenSR[1:3], rho, plotformat('b',2,'-',''))
    #     p1 = [minLenSR[1],0]
    #     utils.PlotLineSeg(p1, finPos, plotformat('c',2,'--',''))
    #     utils.PlotArrow(finPos, finHdng, 1, plotformat('m',2,'--',''))
    #     print('Min Length SR: ',minLenSR[0])
    #     print(f"{finHdng=}")

    # plt.axis('equal')
    # plt.show() 


    # lineSeg = np.array([ (-2, 3), (-1,-4)])
    # finHdng = 4
    # finPos, iniHdng, minLenSR = LocalMinSR2(lineSeg, finHdng, rho)
    
    # utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'-',''))
    # if np.isfinite(finPos[0]):
    #     du.PlotDubPathSegments([0,0,iniHdng], 'SR', minLenSR[1:3], rho, plotformat('b',2,'-',''))
    #     utils.PlotArrow(finPos, finHdng, 1, plotformat('m',2,'--',''))
    #     print('Min Length SR: ',minLenSR[0])
    #     print(f"{iniHdng=}")

    # plt.axis('equal')
    # plt.show()

    ############################# Test RL #############################

    # lineSeg = np.array([ (2,-1), (-2,3)])
    # finPt, finHdngRL, lensLocalMinRL = LocalMinRL(lineSeg, rho)
    
    # C1 = np.array([0,-rho])
    # utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'--',''))
    # utils.PlotCircle(C1, rho,plotformat('g',2,'--',''))
    # du.PlotDubPathSegments(iniConf, 'RL', lensLocalMinRL[1:3], rho, plotformat('b',2,'',''))
    
    # plt.axis('equal')
    # plt.show()

    # du.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'--',''))
    # du.PlotCircle(np.array([0,rho]), rho,plotformat('g',2,'--',''))

    # for lengths in lengthsList:
    #     du.PlotDubPathSegments(iniConf, 'LR', lengths[1:3], rho, plotformat('b',2,'',''))
    # # du.PlotDubPathSegments(iniConf, 'LR', lengthsList[1][1:3], rho, plotformat('b',2,'',''))

    # for finPt in finPtsList:
    #     center2 = finPt + rho*np.array([sin(finHdng), -cos(finHdng)])
    #     du.PlotCircle(center2, rho,plotformat('g',2,'--',''))
    #     du.PlotArrow(finPt, finHdng, 1, plotformat('m',2,'-',''))
    #     C2 = finPt + rho*np.array([cos(finHdng-pi/2), sin(finHdng-pi/2)])
    #     P1 = (C1+C2)*.5
    #     du.PlotLineSeg(P1, finPt, plotformat('g',2,'-','x'))

    ############################# Test LSL #############################
    # lineSeg = np.array([(3, -2) , (2, 3)])
    # finHdng = 1
    # finPt, lengthLSL, pathLSL = LocalMinLSL(lineSeg, finHdng, rho)
    # if isfinite(finPt[0]):
    #     du.PlotDubinsPath(pathLSL,plotformat('b',2,'-',''))
    #     utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'-','x'))
    #     utils.PlotArrow(finPt, finHdng, 1, plotformat('m',2,'-',''))
    #     plt.axis('equal')
    #     plt.show()

    # ############################# Test RSR #############################
    # lineSeg = np.array([(1,4) , (4, 2)])
    # finHdng = -2
    # finPt, lengthRSR, pathRSR = LocalMinRSR(lineSeg, finHdng, rho)
    # print("finPt: ",finPt)
    # if isfinite(finPt[0]):
    #     du.PlotDubinsPath(pathRSR, plotformat('b',2,'-',''))
    #     utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'-','x'))
    #     utils.PlotArrow(finPt, finHdng, 1, plotformat('m',2,'-',''))
    #     plt.axis('equal')
    #     plt.show()
    ############################# Test LSR #############################
    # lineSeg = np.array([(4,4) , (4, -2)])
    # finPt, lengthLSR, pathLSR = LocalMinLSR(lineSeg, finHdng, rho)
    # utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'-','x'))
    # utils.PlotCircle(np.array([0,rho]), rho,plotformat('g',2,'--',''))
    
    # if isfinite(finPt[0]):
    #     du.PlotDubinsPath(pathLSR,plotformat('b',2,'-',''))
    #     utils.PlotArrow(finPt, finHdng, 1, plotformat('m',2,'--',''))

    # plt.axis('equal')
    # plt.show()

    ############################# Test RSL #############################
    # lineSeg = np.array([(4,4) , (6, 0)])
    # finPt, lengthRSL, pathRSL = LocalMinRSL(lineSeg, finHdng, rho)

    # utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'-','x'))
    # utils.PlotCircle(np.array([0,rho]), rho,plotformat('g',2,'--',''))
    
    # if isfinite(finPt[0]):
    #     du.PlotDubinsPath(pathRSL,plotformat('b',2,'-',''))
    #     utils.PlotArrow(finPt, finHdng, 1, plotformat('m',2,'--',''))

    # plt.axis('equal')
    # plt.show()

    ############################# Test LRL feasibility to LineSeg #############################
    # lineSeg = np.array([ (-4,-4), (-1,6)])
    # finHdng = 4
    # finPtsList, lengthsList = LRLFeasLimits(lineSeg, finHdng, rho)
    
    # C1 = np.array([0,rho])
    
    # utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'--',''))
    # utils.PlotCircle(C1, rho,plotformat('g',2,'--',''))

    # for k in range(len(finPtsList)):        
    #     finPt = finPtsList[k]
    #     C2 = finPt + rho*np.array([-sin(finHdng), cos(finHdng)])
    #     utils.PlotCircle(C2, rho,plotformat('g',2,'--',''))

    #     # pathLRL = dubins.path([0,0,0],[finPt[0], finPt[1],finHdng+.001],rho, 5)
    #     # segmentLengths = [pathLRL.segment_length(0), pathLRL.segment_length(1), pathLRL.segment_length(2)]
    #     # print(f"{segmentLengths=}")
    #     print(f"{lengthsList[k]=}")

    #     # du.PlotDubinsPath(pathLRL,plotformat('b',2,'-',''))
    #     du.PlotDubPathSegments(iniConf, 'LRL', lengthsList[k], rho, plotformat('b',2,'-',''))
    #     utils.PlotArrow(finPt, finHdng, 1, plotformat('m',2,'--',''))
        
    # plt.axis('equal')
    # plt.show()

    ############################# Test RLR feasibility to LineSeg #############################
    # lineSeg = np.array([ (-3,3), (5,1)])
    # finHdng = 1
    # finPtsList = RLRFeasLimits(lineSeg, finHdng, rho)
    
    # C1 = np.array([0,-rho])
    
    # utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'--',''))
    # utils.PlotCircle(C1, rho,plotformat('g',2,'--',''))

    # for k in range(len(finPtsList)):        
    #     finPt = finPtsList[k]
    #     C2 = finPt + rho*np.array([sin(finHdng), -cos(finHdng)])
    #     utils.PlotCircle(C2, rho,plotformat('g',2,'--',''))

    #     pathRLR = dubins.path([0,0,0],[finPt[0], finPt[1],finHdng-.000001],rho, 4)
    #     midSegLength = pathRLR.segment_length(1)
    #     print(f"{midSegLength=}")
    #     du.PlotDubinsPath(pathRLR, plotformat('b',2,'-',''))
    #     utils.PlotArrow(finPt, finHdng, 1, plotformat('m',2,'--',''))
        
    # plt.axis('equal')
    # plt.show()

# ############################# Test LS to line #############################
    # lineSeg = np.array([ (4,3), (-1,-3)])
    # lineSeg = np.array([ (2,-3), (3,4) ])

    # finHdng = 1
    # finPos, lengthsLS = LStoLine(lineSeg, finHdng, rho)
    # C1 = np.array([0,rho])
    
    # utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'--',''))
    # utils.PlotCircle(C1, rho,plotformat('g',2,'--',''))

    # if isfinite(finPos[0]):
    #     du.PlotDubPathSegments([0,0,0],'LS',lengthsLS[1:3],rho,plotformat('b',2,'-',''))
    #     utils.PlotArrow(finPos, finHdng, 1, plotformat('m',2,'--',''))
    #     print(f"{lengthsLS=}")
        
        
    # plt.axis('equal')
    # plt.show()

# ############################# Test RS to line #############################
#     # lineSeg = np.array([ (4,3), (-1,-3)])
#     lineSeg = np.array([ (-2,3), (4,2) ])

#     finHdng = 1.8
#     finPos, lengthsRS = RStoLine(lineSeg, finHdng, rho)
#     C1 = np.array([0,-rho])
    
#     utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'--',''))
#     utils.PlotCircle(C1, rho,plotformat('g',2,'--',''))

#     if isfinite(finPos[0]):
#         du.PlotDubPathSegments([0,0,0],'RS',lengthsRS[1:3],rho,plotformat('b',2,'-',''))
#         utils.PlotArrow(finPos, finHdng, 1, plotformat('m',2,'--',''))
#         print(f"{lengthsRS=}")
        
        
#     plt.axis('equal')
#     plt.show()

############################# Test SL to line #############################
    # # lineSeg = np.array([ (4,3), (-1,-3)])
    # lineSeg = np.array([ (0,2), (4,3) ])

    # finHdng = 2
    # finPos, lengthsSL = SLtoLine(lineSeg, finHdng, rho)
    
    # utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'--',''))
    
    # if isfinite(finPos[0]):
    #     du.PlotDubPathSegments([0,0,0],'SL',lengthsSL[1:3],rho,plotformat('b',2,'-',''))
    #     utils.PlotArrow(finPos, finHdng, 1, plotformat('m',2,'--',''))
    #     print(f"{lengthsSL=}")
        
        
    # plt.axis('equal')
    # plt.show()

############################# Test SR to line #############################
    # lineSeg = np.array([ (4,3), (-1,-3)])
    # lineSeg = np.array([ (0,-3), (2.19,1) ])

    # finHdng = 4.5
    # finPos, lengthsSR = SRtoLine(lineSeg, finHdng, rho)
    
    # utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'--',''))
    
    # if isfinite(finPos[0]):
    #     du.PlotDubPathSegments([0,0,0],'SR',lengthsSR[1:3],rho,plotformat('b',2,'-',''))
    #     utils.PlotArrow(finPos, finHdng, 1, plotformat('m',2,'--',''))
    #     print(f"{lengthsSR=}")
        
        
    # plt.axis('equal')
    # plt.show()