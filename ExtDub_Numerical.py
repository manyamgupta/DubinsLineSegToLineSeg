# Shortest Dubins LSR path to a line segment
# Author: Satya Gupta Manyam
import numpy as np
from numpy import pi,cos,sin
import matplotlib.pyplot as plt
import dubins 
import dubutils as du
import utils
import DubToLineSeg as dls
from collections import namedtuple


if __name__ == "__main__":

    LSL =0; LSR = 1; RSL = 2; RSR = 3; RLR = 4; LRL = 5; 
    pathType = LRL
    plotformat = namedtuple("plotformat","color linewidth linestyle marker")
    iniConf = np.array([0,0,0])
    rho = 1
    finHdng = 1
    
    ndisc = 10000
    lamVec = np.linspace(0,1,ndisc)
    lenVec = np.zeros(ndisc)
    finHdngVec = np.linspace(0,2*pi,ndisc)
    # utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'-','x') )
    # utils.PlotCircle( np.array([0,rho]), rho, plotformat('g',2,'--','') )
    # plt.axis('equal')
    # for k in range(ndisc):
    #     lam = lamVec[k]
    #     finPos = (1-lam)*lineSeg[0]+lam*lineSeg[1]
    #     pathDub = dubins.path([0,0,0],[finPos[0], finPos[1],finHdng],rho, pathType)
    #     if pathDub == None:
    #         lenVec[k] = np.nan
    #     else:
    #         lenVec[k] = pathDub.path_length()
    #         # du.PlotDubinsPath(pathDub,plotformat('b',2,'-',''))
    # plt.figure()
    # plt.plot(lamVec, lenVec)
    
    # lenVecFinite = lenVec[np.isfinite(lenVec)]
    # lamVecFinite = lamVec[np.isfinite(lenVec)]
    # minInd = np.argmin(lenVecFinite)
    # minLam = lamVecFinite[minInd]
    # # minLam = lamVecFinite[-1]
    
    # finPos = (1-minLam)*lineSeg[0]+minLam*lineSeg[1]
    # plt.figure()
    # utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'-','x') )
    # utils.PlotCircle( np.array([0,rho]), rho, plotformat('g',2,'--','') )
    # pathDub = dubins.path([0,0,0],[finPos[0], finPos[1],finHdng],rho, pathType)        
    # du.PlotDubinsPath(pathDub,plotformat('b',2,'-',''))
    # utils.PlotArrow(finPos, finHdng, 1, plotformat('m',2,'--','') )
    
    lineSeg = np.array([ (-1,-2.5), (3,-1.5)])
    for k in range(ndisc):
        lam = lamVec[k]
        finHdng = finHdngVec[k]
        finPos, lengths = dls.SRtoLine(lineSeg, finHdng, rho)
        lenVec[k] = lengths[0]

    plt.figure()
    plt.plot(finHdngVec, lenVec)
    
    lenVecFinite = lenVec[np.isfinite(lenVec)]
    finHdngVecFinite = finHdngVec[np.isfinite(lenVec)]
    minInd = np.argmin(lenVecFinite)
    minFinHdng = finHdngVecFinite[minInd]
    minFinPos, minLengths = dls.SRtoLine(lineSeg, minFinHdng, rho)
    print(f"{minFinHdng=}")
    print(f"{minLengths=}")
    
    plt.figure()
    utils.PlotLineSeg(lineSeg[0], lineSeg[1], plotformat('c',2,'-','x') )
    # utils.PlotCircle( np.array([0,-rho]), rho, plotformat('g',2,'--','') )
    utils.PlotArrow(minFinPos, minFinHdng, 1, plotformat('m',2,'--','') )
    du.PlotDubPathSegments([0,0,0],'SR',minLengths[1:3],rho,plotformat('b',2,'-',''))

    plt.axis('equal')
    plt.show()


