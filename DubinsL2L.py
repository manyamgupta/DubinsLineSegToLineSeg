# Author: Satyanarayana Gupta Manyam
# Udpate Date: Feb 8th, 2023
# This file contains the Line2LineDubins Class and the functions necessary to compute the
# minimum Dubins Line to Line path
# The inputs to the class are the line1, interval1, line2, interval2, rho
# The function 'MinDub_L2L' computes the minimum path and its length

import numpy as np
from numpy import pi,cos,sin
import matplotlib.pyplot as plt
import dubins
import utils
import dubutils as du
from types import SimpleNamespace
from dataclasses import dataclass
from timeit import default_timer as timer
import copy

@dataclass
class Config:
    pos_x: float # position: x-coordinate
    pos_y: float # position: y-coordinate
    heading: float # heading

@dataclass
class LineSegment:
    point1: tuple # tuple of x,y coordinates of end point
    point2: tuple # tuple of x,y coordinates of end point
    
    
@dataclass
class AngIntrvl:
    # interval goes from lower bound to upper bound in counter-clockwise direction
    int_lb: float # lower bound of the interval
    int_ub: float # upper bound of the interval   
        
@dataclass
class CandidatePath:
    pathType: str 
    iniPos: tuple # position on line 1
    iniHead: float # initial heading
    finalPos: tuple # position on line 2
    finalHead: float # final heading
    segLengths: tuple
    
class Line2LineDubins:
    
    def __init__(self, line1, int1, line2, int2, rho):   
        
        self.line1 = copy.copy(line1)
        self.int1 = copy.copy(int1)        
        self.line2 = copy.copy(line2)
        self.int2 = copy.copy(int2)        
        self.rho = rho 
        
        
        # translated and rotates lines
        # translates line 1 such that line2_point1_int_lb moved to [0,0,0]
        self.line1_trans1 = self.RotateTransLine(Config(self.line2.point1[0], self.line2.point1[1], self.int2[0]+np.pi), self.line1)
        # translates line 1 such that line2_point1_int_ub moved to [0,0,0]
        self.line1_trans2 = self.RotateTransLine(Config(self.line2.point1[0], self.line2.point1[1], self.int2[1]+np.pi), self.line1)
        # translates line 1 such that line2_point2_int_lb moved to [0,0,0]
        self.line1_trans3 = self.RotateTransLine(Config(self.line2.point2[0], self.line2.point2[1], self.int2[0]+np.pi), self.line1)
        # translates line 1 such that line2_point2_int_ub moved to [0,0,0]
        self.line1_trans4 = self.RotateTransLine(Config(self.line2.point2[0], self.line2.point2[1], self.int2[1]+np.pi), self.line1)
        
        # translates line 2 such that line1_point1_int_lb moved to [0,0,0]
        self.line2_trans1 = self.RotateTransLine(Config(self.line1.point1[0], self.line1.point1[1], self.int1[0]), self.line2)
        # translates line 2 such that line1_point1_int_ub moved to [0,0,0]
        self.line2_trans2 = self.RotateTransLine(Config(self.line1.point1[0], self.line1.point1[1], self.int1[1]), self.line2)
        # translates line 2 such that line1_point2_int_lb moved to [0,0,0]
        self.line2_trans3 = self.RotateTransLine(Config(self.line1.point2[0], self.line1.point2[1], self.int1[0]), self.line2)
        # translates line 2 such that line1_point2_int_ub moved to [0,0,0]
        self.line2_trans4 = self.RotateTransLine(Config(self.line1.point2[0], self.line1.point2[1], self.int1[1]), self.line2)
                
        # empty lists for lengths and candidate paths, these are populated when A2AMinDubins() is called
        self.lengthsVec = []
        self.candPathsList = []
        
        return
    
    def RotateTransLine(self, startConfig, line):
        # Rotates and translates the line such that the start config aligns with [0,0,0]
        
        pt1_trans = (line.point1[0]-startConfig.pos_x, line.point1[1]-startConfig.pos_y)
        pt2_trans = (line.point2[0]-startConfig.pos_x, line.point2[1]-startConfig.pos_y)
        pt1_transrot = utils.RotateVec(pt1_trans, -startConfig.heading )
        pt2_transrot = utils.RotateVec(pt2_trans, -startConfig.heading )
        
        return LineSegment(pt1_transrot, pt2_transrot)
        
    def MinDub_L2L(self):
        # Computes the list of candidate paths for minimum line to line dubins paths        
        # returns the length of the minimum path, minimum path, and candidate paths
        # theta_1_l = self.int1[0]
        # theta_1_u = self.int1[1]
        # theta_2_l = self.int2[0]
        # theta_2_u = self.int2[1]
        
        ######################## Paths b/w Boundary points, interval boundaries ########################
        self.AddBoundaryPaths()
                
        ######################## LSL Paths ########################
        self.AddLSLPaths()

        ######################## RSR Paths ########################
        self.AddRSRPaths()

        ######################## LSR Paths ########################
        self.AddLSRPaths()

        ######################## RSL Paths ########################
        self.AddRSLPaths()

        ######################## LS Paths ########################
        self.AddLSPaths()
        
        ######################## RS Paths ########################           
        self.AddRSPaths()
        
        ######################## SR Paths ########################                   
        self.AddSRPaths()
        
        ######################## SL Paths ########################                   
        self.AddSLPaths()
                
        ######################## LR Paths ########################                   
        self.AddLRPaths()

        ######################## RL Paths ########################                   
        self.AddRLPaths()
        
        ######################## L Paths ########################                   
        self.AddLPaths()

        ######################## R Paths ########################                   
        self.AddRPaths()        
        
        ######################## S Paths ########################                   
        self.AddSPaths() 
        
        lengthsVec = np.array(self.lengthsVec)
        if lengthsVec.size >0:
            minLength = np.min(lengthsVec)
            minInd = np.argmin(lengthsVec)
            minPath = self.candPathsList[minInd]  
        else:
            return None, None
        return minLength, minPath
    
    def AddCandToList(self, pathType, pos_start, theta_start, finPos, finHdng, candLen, candSegLens):
        
        # finPos, finHdng are in the transformed/rotated frame
        # pos_start and theta_start are the actual start position and start heading
        # This functions rotate the final position and heading and adds to the list of candidate paths
        finPos_rotback = utils.RotateVec(finPos, theta_start)
        finPos_transrot = (finPos_rotback[0]+pos_start[0], finPos_rotback[1]+pos_start[1])
        cp = CandidatePath(pathType, pos_start, theta_start, finPos_transrot, finHdng+theta_start, candSegLens)
        self.candPathsList.append(cp)
        self.lengthsVec.append(candLen)  
        
        return cp

    def AddCandToList2(self, pathType, pos_start, theta_start, finPos, finHdng, candLen, candSegLens):
        
        # finPos, finHdng are in the transformed/rotated frame, and belong to start config of the fwd path
        # pos_start and theta_start-pi are the actual final position and heading
        # This function rotates and reverses the paths and add to the candidate paths
        finPos_rotback = utils.RotateVec(finPos, theta_start)
        finPos_transrot = (finPos_rotback[0]+pos_start[0], finPos_rotback[1]+pos_start[1])
        cp = CandidatePath(pathType, finPos_transrot, finHdng+theta_start-np.pi, pos_start, theta_start-np.pi, candSegLens[::-1])
        self.candPathsList.append(cp)
        self.lengthsVec.append(candLen)  
        
        return cp
    def AddBoundaryPaths(self):
        
        for p1 in (self.line1.point1, self.line1.point2):
            for t1 in self.int1:
                for p2 in (self.line2.point1, self.line2.point2):
                    for t2 in self.int2:
                        dPath = dubins.shortest_path([p1[0], p1[1], t1],[p2[0], p2[1], t2], self.rho)
                        candSegLens = (dPath.segment_length(0), dPath.segment_length(1), dPath.segment_length(2))
                        pathType = du.PathNumtoType(dPath.path_type())
                        cp = CandidatePath(pathType, p1, t1, p2, t2, candSegLens)
                        self.candPathsList.append(cp)
                        self.lengthsVec.append(dPath.path_length())  
        return
    
    def AddLSLPaths(self):
        candLen = None
        cp = None        
        theta_1_u = self.int1[1]
        theta_2_l = self.int2[0]
        
        # LSL from line1_pt1, theta_1_u to theta_2_l
        finHdng = theta_2_l-theta_1_u
        finPos, candLen, candSegLens = self.LocalMinLSL(self.line2_trans2, finHdng)
        if candLen:
            cp = self.AddCandToList('LSL', self.line1.point1, theta_1_u, finPos, finHdng, candLen, candSegLens)
        
        # LSL from line1_pt2, theta_1_u to theta_2_l
        # final heading same as the one before 
        finPos, candLen, candSegLens = self.LocalMinLSL(self.line2_trans4, finHdng)
        if candLen:
            cp = self.AddCandToList('LSL', self.line1.point2, theta_1_u, finPos, finHdng, candLen, candSegLens)
            # self.PlotScenario('LSL', (0,0,0), candSegLens, self.line2_trans4)      

        # LSL from theta_1_u to line2_pt1, theta_2_l
        # This is calculated using RSR from final to initial
        line1_trans = self.RotateTransLine(Config(self.line2.point1[0], self.line2.point1[1], theta_2_l+np.pi), self.line1)
        finPos, candLen, candSegLens = self.LocalMinRSR(line1_trans, theta_1_u-theta_2_l)
        # self.PlotScenario('RSR', (0,0,0), candSegLens, line1_trans)   
        # plt.show()   
        if candLen:
            cp = self.AddCandToList2('LSL', self.line2.point1, theta_2_l+np.pi, finPos, theta_1_u-theta_2_l, candLen, candSegLens)

        # LSL from theta_1_u to line2_pt2, theta_2_l
        # This is calculated using RSR from final to initial
        line1_trans = self.RotateTransLine(Config(self.line2.point2[0], self.line2.point2[1], theta_2_l+np.pi), self.line1)
        finPos, candLen, candSegLens = self.LocalMinRSR(line1_trans, theta_1_u-theta_2_l)  
        if candLen:
            cp = self.AddCandToList2('LSL', self.line2.point2, theta_2_l+np.pi, finPos, theta_1_u-theta_2_l, candLen, candSegLens)
    
        return candLen, cp
    
    def AddRSRPaths(self):
        
        theta_1_l = self.int1[0]
        theta_2_u = self.int2[1]
        
        # RSR local min from line1_pt1, theta_1_l to theta_2_u
        finHdng = theta_2_u-theta_1_l
        candLen = None
        cp = None
        finPos, candLen, candSegLens = self.LocalMinRSR(self.line2_trans1, finHdng)
        if candLen:
            cp = self.AddCandToList('RSR', self.line1.point1, theta_1_l, finPos, finHdng, candLen, candSegLens)

        # RSR local min from line1_pt2, theta_1_l to theta_2_u
        finHdng = theta_2_u-theta_1_l
        finPos, candLen, candSegLens = self.LocalMinRSR(self.line2_trans3, finHdng)
        if candLen:
            cp = self.AddCandToList('RSR', self.line1.point2, theta_1_l, finPos, finHdng, candLen, candSegLens)

        # RSR local min from theta_1_l to line2_pt1, theta_2_u
        # This is calculated using LSL from final to initial
        line1_trans = self.RotateTransLine(Config(self.line2.point1[0], self.line2.point1[1], theta_2_u+np.pi), self.line1)
        finPos, candLen, candSegLens = self.LocalMinLSL(line1_trans, theta_1_l-theta_2_u)
        if candLen:
            cp = self.AddCandToList2('RSR', self.line2.point1, theta_2_u+np.pi, finPos, theta_1_l-theta_2_u, candLen, candSegLens)
            # self.PlotScenario('LSL', (0,0,0), candSegLens, line1_trans)   
            # plt.show()    

        # RSR local min from theta_1_l to line2_pt2, theta_2_u
        # This is calculated using LSL from final to initial
        line1_trans = self.RotateTransLine(Config(self.line2.point2[0], self.line2.point2[1], theta_2_u+np.pi), self.line1)
        finPos, candLen, candSegLens = self.LocalMinLSL(line1_trans, theta_1_l-theta_2_u)
        if candLen:
            cp = self.AddCandToList2('RSR', self.line2.point2, theta_2_u+np.pi, finPos, theta_1_l-theta_2_u, candLen, candSegLens)

        return candLen, cp
    
    def AddLSRPaths(self):
        
        candLen = None
        cp = None        
        theta_1_u = self.int1[1]        
        theta_2_u = self.int2[1]
        
        # LSR from line1_pt1, theta_1_u to theta_2_u
        finHdng = theta_2_u-theta_1_u
        finPos, candLen, candSegLens = self.LocalMinLSR(self.line2_trans2, finHdng)
        if candLen:
            cp = self.AddCandToList('LSR', self.line1.point1, theta_1_u, finPos, finHdng, candLen, candSegLens)

        # LSR from line1_pt2, theta_1_u to theta_2_u
        finHdng = theta_2_u-theta_1_u
        finPos, candLen, candSegLens = self.LocalMinLSR(self.line2_trans4, finHdng)
        if candLen:
            cp = self.AddCandToList('LSR', self.line1.point2, theta_1_u, finPos, finHdng, candLen, candSegLens)

        # LSR from theta_1_u to line2_pt1, theta_2_u
        # This is calculated using LSR from final to initial        
        finPos, candLen, candSegLens = self.LocalMinLSR(self.line1_trans2, theta_1_u-theta_2_u) 
        if candLen:
            cp = self.AddCandToList2('LSR', self.line2.point1, theta_2_u+np.pi, finPos, theta_1_u-theta_2_u, candLen, candSegLens)

        # LSR from theta_1_u to line2_pt2, theta_2_u
        # This is calculated using LSR from final to initial        
        finPos, candLen, candSegLens = self.LocalMinLSR(self.line1_trans4, theta_1_u-theta_2_u) 
        if candLen:
            cp = self.AddCandToList2('LSR', self.line2.point2, theta_2_u+np.pi, finPos, theta_1_u-theta_2_u, candLen, candSegLens)

        return candLen, cp
    
    def AddRSLPaths(self):
        candLen = None
        cp = None        
        theta_1_l = self.int1[0]        
        theta_2_l = self.int2[0]
        
        # RSL from line1_pt1, theta_1_l to theta_2_l
        finHdng = theta_2_l-theta_1_l
        finPos, candLen, candSegLens = self.LocalMinRSL(self.line2_trans1, finHdng)
        if candLen:
            cp = self.AddCandToList('RSL', self.line1.point1, theta_1_l, finPos, finHdng, candLen, candSegLens)

        # RSL from line1_pt2, theta_1_l to theta_2_l
        finHdng = theta_2_l-theta_1_l
        finPos, candLen, candSegLens = self.LocalMinRSL(self.line2_trans3, finHdng)
        if candLen:
            cp = self.AddCandToList('RSL', self.line1.point2, theta_1_l, finPos, finHdng, candLen, candSegLens)

        # RSL from theta_1_l to line2_pt1, theta_2_l
        # This is calculated using RSL from final to initial        
        finPos, candLen, candSegLens = self.LocalMinRSL(self.line1_trans1, theta_1_l-theta_2_l) 
        if candLen:
            cp = self.AddCandToList2('RSL', self.line2.point1, theta_2_l+np.pi, finPos, theta_1_l-theta_2_l, candLen, candSegLens)

        # RSL from theta_1_l to line2_pt2, theta_2_l
        # This is calculated using RSL from final to initial        
        finPos, candLen, candSegLens = self.LocalMinRSL(self.line1_trans3, theta_1_l-theta_2_l) 
        if candLen:
            cp = self.AddCandToList2('RSL', self.line2.point2, theta_2_l+np.pi, finPos, theta_1_l-theta_2_l, candLen, candSegLens)

        return candLen, cp

    def AddSRPaths(self):
        L2LDubRev = Line2LineDubins(self.line2, (self.int2[0]+np.pi,self.int2[1]+np.pi), self.line1, (self.int1[0]+np.pi,self.int1[1]+np.pi), self.rho) 
        L2LDubRev.AddLSPaths()
        for cp_r in L2LDubRev.candPathsList:            
            cp = CandidatePath('SR', cp_r.finalPos , np.mod(cp_r.finalHead+np.pi, 2*np.pi), cp_r.iniPos, np.mod(cp_r.iniHead+np.pi, 2*np.pi), cp_r.segLengths[::-1])
            self.candPathsList.append(cp)
            self.lengthsVec.append(cp_r.segLengths[0]+cp_r.segLengths[1])  
            # self.PlotScenario('SR', (cp.iniPos[0], cp.iniPos[1], cp.iniHead), cp.segLengths, self.line2)   
            # plt.show()
            
        return

    def AddSLPaths(self):
        L2LDubRev = Line2LineDubins(self.line2, (self.int2[0]+np.pi,self.int2[1]+np.pi), self.line1, (self.int1[0]+np.pi,self.int1[1]+np.pi), self.rho) 
        L2LDubRev.AddRSPaths()
        for cp_r in L2LDubRev.candPathsList:            
            cp = CandidatePath('SL', cp_r.finalPos , np.mod(cp_r.finalHead+np.pi, 2*np.pi), cp_r.iniPos, np.mod(cp_r.iniHead+np.pi, 2*np.pi), cp_r.segLengths[::-1])
            self.candPathsList.append(cp)
            self.lengthsVec.append(cp_r.segLengths[0]+cp_r.segLengths[1])  
            # self.PlotScenario('SL', (cp.iniPos[0], cp.iniPos[1], cp.iniHead), cp.segLengths, self.line2)   
            # plt.show()
            
        return
        
    def AddLSPaths(self):
        
        theta_1_u = self.int1[1]   
        
        # LS path from line1 pt1, theta_1_u to line2 pt1
        finHead, candLen, candSegLens = self.LSpt2ptPaths(self.line1.point1, theta_1_u, self.line2.point1)
        self.VerAddTwoSegPath('LS', self.line1.point1, theta_1_u, self.line2.point1, finHead, candLen, candSegLens)        
        # LS path from line1 pt1, theta_1_u to line2 pt2             
        finHead, candLen, candSegLens = self.LSpt2ptPaths(self.line1.point1, theta_1_u, self.line2.point2)
        self.VerAddTwoSegPath('LS', self.line1.point1, theta_1_u, self.line2.point2, finHead, candLen, candSegLens)
        # LS path from line1 pt2, theta_1_u to line2 pt1
        finHead, candLen, candSegLens = self.LSpt2ptPaths(self.line1.point2, theta_1_u, self.line2.point1)
        self.VerAddTwoSegPath('LS', self.line1.point2, theta_1_u, self.line2.point1, finHead, candLen, candSegLens)
        # LS path from line1 pt2, theta_1_u to line2 pt2
        finHead, candLen, candSegLens = self.LSpt2ptPaths(self.line1.point2, theta_1_u, self.line2.point2)
        self.VerAddTwoSegPath('LS', self.line1.point2, theta_1_u, self.line2.point2, finHead, candLen, candSegLens)

        # LS path from line1 pt1, theta_1_u to line2 theta_2_l
        finPos, candLen, candSegLens = self.LSpt2LineT2(self.line1.point1, theta_1_u, self.line2, self.int2[0])
        self.VerAddTwoSegPath('LS', self.line1.point1, theta_1_u, finPos, self.int2[0], candLen, candSegLens)        
        # LS path from line1 pt1, theta_1_u to line2 theta_2_u
        finPos, candLen, candSegLens = self.LSpt2LineT2(self.line1.point1, theta_1_u, self.line2, self.int2[1])
        self.VerAddTwoSegPath('LS', self.line1.point1, theta_1_u, finPos, self.int2[1], candLen, candSegLens)
        # LS path from line1 pt2, theta_1_u to line2 theta_2_l
        finPos, candLen, candSegLens = self.LSpt2LineT2(self.line1.point2, theta_1_u, self.line2, self.int2[0])
        self.VerAddTwoSegPath('LS', self.line1.point2, theta_1_u, finPos, self.int2[0], candLen, candSegLens)
        # LS path from line1 pt2, theta_1_u to line2 theta_2_u
        finPos, candLen, candSegLens = self.LSpt2LineT2(self.line1.point2, theta_1_u, self.line2, self.int2[1])
        self.VerAddTwoSegPath('LS', self.line1.point2, theta_1_u, finPos, self.int2[1], candLen, candSegLens)

        # local min LS (wrt lam_2) path from line1 pt1, theta_1_u to line2 
        finHdng= self.LocalMinLS(self.line2_trans2)
        finHdng = finHdng+theta_1_u
        if utils.InInt(self.int2[0], self.int2[1], finHdng ):
            finPos, candLen, candSegLens = self.LSpt2LineT2(self.line1.point1, theta_1_u, self.line2, finHdng)
            self.VerAddTwoSegPath('LS', self.line1.point1, theta_1_u, finPos, finHdng, candLen, candSegLens)
        
        # local min LS (wrt lam_2) path from line1 pt2, theta_1_u to line2 
        finHdng= self.LocalMinLS(self.line2_trans4)
        finHdng = finHdng+theta_1_u
        if utils.InInt(self.int2[0], self.int2[1], finHdng ):
            finPos, candLen, candSegLens = self.LSpt2LineT2(self.line1.point2, theta_1_u, self.line2, finHdng)
            self.VerAddTwoSegPath('LS', self.line1.point2, theta_1_u, finPos, finHdng, candLen, candSegLens)
        
        # local min LS (wrt lam_1) path from (line1, theta_1_u) to line2 pt1 or pt2
        finHdngs = self.LocalMinLS2(self.line1)
        for finHdng in finHdngs:
            if utils.InInt(self.int2[0], self.int2[1], finHdng ):
                startHdng = finHdng+np.pi
                for startPt in (self.line2.point1, self.line2.point2):
                    line1_trans= self.RotateTransLine(Config(startPt[0], startPt[1], startHdng), self.line1)
                    finPos, candLen, candSegLens = self.SRtoLine(line1_trans, theta_1_u-startHdng+np.pi)
                    if candLen:
                        finPos_r = utils.RotateVec(finPos, startHdng )
                        finPos_tr = (finPos_r[0]+startPt[0], finPos_r[1]+startPt[1])                
                        self.VerAddTwoSegPath('LS', finPos_tr, theta_1_u, startPt, finHdng, candLen, candSegLens[::-1])
        
        return 
    
    def LSpt2LineT2(self, pt1, t1, targLine, t2):
        
        targLine_trans = self.RotateTransLine(Config(pt1[0], pt1[1], t1), targLine)
        finalHead = t2-t1       
        finPos, lengthLS, segLengths = self.LStoLineT2(targLine_trans, finalHead)
        
        if lengthLS:
            finPos_r = utils.RotateVec(finPos, t1 )
            finPos_tr = (finPos_r[0]+pt1[0], finPos_r[1]+pt1[1])
            
        else:
            return (None, None), None, (None, None)
       
        return finPos_tr, lengthLS, segLengths
    
    def LSpt2ptPaths(self, pt1, t1, pt2):
             
        # LS from pt1, t1 to pt2
        finHead = None
        pt2_tans = (pt2[0]-pt1[0], pt2[1]-pt1[1])
        pt2_transrot = utils.RotateVec(pt2_tans, -t1 )

        lenLS, segLens = du.PathLS(pt2_transrot[0], pt2_transrot[1], self.rho)        
        if lenLS:
            finHead = segLens[0]/self.rho+t1
       
        return finHead, lenLS, segLens
    
    def VerAddTwoSegPath(self, pathType, p1, t1, p2, t2, candLen, candSegLens):
        if candLen != None:
            if utils.InInt(self.int1[0], self.int1[1], t1 ) and  utils.InInt(self.int2[0], self.int2[1], t2 ):
                cp = CandidatePath(pathType, p1, t1, p2, t2, candSegLens)
                self.candPathsList.append(cp)
                self.lengthsVec.append(candLen)          
        
        return
    
    def AddRSPaths(self):
        
        theta_1_l = self.int1[0]   
        
        # RS path from line1 pt1, theta_1_l to line2 pt1
        finHead, candLen, candSegLens = self.RSpt2ptPaths(self.line1.point1, theta_1_l, self.line2.point1)
        self.VerAddTwoSegPath('RS', self.line1.point1, theta_1_l, self.line2.point1, finHead, candLen, candSegLens)        
        # RS path from line1 pt1, theta_1_l to line2 pt2
        finHead, candLen, candSegLens = self.RSpt2ptPaths(self.line1.point1, theta_1_l, self.line2.point2)
        self.VerAddTwoSegPath('RS', self.line1.point1, theta_1_l, self.line2.point2, finHead, candLen, candSegLens)
        # RS path from line1 pt2, theta_1_l to line2 pt1
        finHead, candLen, candSegLens = self.RSpt2ptPaths(self.line1.point2, theta_1_l, self.line2.point1)
        self.VerAddTwoSegPath('RS', self.line1.point2, theta_1_l, self.line2.point1, finHead, candLen, candSegLens)
        # RS path from line1 pt2, theta_1_l to line2 pt2
        finHead, candLen, candSegLens = self.RSpt2ptPaths(self.line1.point2, theta_1_l, self.line2.point2)
        self.VerAddTwoSegPath('RS', self.line1.point2, theta_1_l, self.line2.point2, finHead, candLen, candSegLens)

        # RS path from line1 pt1, theta_1_l to line2 theta_2_l
        finPos, candLen, candSegLens = self.RSpt2LineT2(self.line1.point1, theta_1_l, self.line2, self.int2[0])
        self.VerAddTwoSegPath('RS', self.line1.point1, theta_1_l, finPos, self.int2[0], candLen, candSegLens)
        # RS path from line1 pt1, theta_1_l to line2 theta_2_u
        finPos, candLen, candSegLens = self.RSpt2LineT2(self.line1.point1, theta_1_l, self.line2, self.int2[1])
        self.VerAddTwoSegPath('RS', self.line1.point1, theta_1_l, finPos, self.int2[1], candLen, candSegLens)
        # RS path from line1 pt2, theta_1_l to line2 theta_2_l
        finPos, candLen, candSegLens = self.RSpt2LineT2(self.line1.point2, theta_1_l, self.line2, self.int2[0])
        self.VerAddTwoSegPath('RS', self.line1.point2, theta_1_l, finPos, self.int2[0], candLen, candSegLens)
        # RS path from line1 pt2, theta_1_l to line2 theta_2_u
        finPos, candLen, candSegLens = self.RSpt2LineT2(self.line1.point2, theta_1_l, self.line2, self.int2[1])
        self.VerAddTwoSegPath('RS', self.line1.point2, theta_1_l, finPos, self.int2[1], candLen, candSegLens)

        # local min RS (wrt lam_2) path from line1 pt1, theta_1_l to line2 
        finHdng= self.LocalMinRS(self.line2_trans1)
        finHdng = finHdng+theta_1_l
        if utils.InInt(self.int2[0], self.int2[1], finHdng ):
            finPos, candLen, candSegLens = self.RSpt2LineT2(self.line1.point1, theta_1_l, self.line2, finHdng)
            self.VerAddTwoSegPath('RS', self.line1.point1, theta_1_l, finPos, finHdng, candLen, candSegLens)
        # local min RS (wrt lam_2) path from line1 pt2, theta_1_l to line2 
        finHdng= self.LocalMinRS(self.line2_trans3)
        finHdng = finHdng+theta_1_l
        if utils.InInt(self.int2[0], self.int2[1], finHdng ):
            finPos, candLen, candSegLens = self.RSpt2LineT2(self.line1.point2, theta_1_l, self.line2, finHdng)
            self.VerAddTwoSegPath('RS', self.line1.point2, theta_1_l, finPos, finHdng, candLen, candSegLens)

        # local min RS (wrt lam_1) path from (line1, theta_1_l) to line2 pt1 or pt2
        finHdngs = self.LocalMinLS2(self.line1) # this function works for both LS and RS, since it outputs heading orthogonal to line 1
        for finHdng in finHdngs:
            if utils.InInt(self.int2[0], self.int2[1], finHdng ):
                startHdng = finHdng+np.pi
                for startPt in (self.line2.point1, self.line2.point2):
                    line1_trans= self.RotateTransLine(Config(startPt[0], startPt[1], startHdng), self.line1)
                    # self.PlotScenario('SR', (0,0,0), (1,0), line1_trans)
                    plt.show()
                    finPos, candLen, candSegLens = self.SLtoLine(line1_trans, theta_1_l-startHdng+np.pi)
                    if candLen:
                        finPos_r = utils.RotateVec(finPos, startHdng )
                        finPos_tr = (finPos_r[0]+startPt[0], finPos_r[1]+startPt[1])                
                        self.VerAddTwoSegPath('RS', finPos_tr, theta_1_l, startPt, finHdng, candLen, candSegLens[::-1])
        
        return
    
    def RSpt2ptPaths(self, pt1, t1, pt2):
             
        # RS from pt1, t1 to pt2
        finHead = None
        pt2_trans = (pt2[0]-pt1[0], pt2[1]-pt1[1])
        pt2_transrot = utils.RotateVec(pt2_trans, -t1 )

        lenRS, segLens = du.PathRS(pt2_transrot[0], pt2_transrot[1], self.rho)        
        if lenRS:
            finHead = np.mod(-segLens[0]/self.rho+t1, 2*np.pi)
       
        return finHead, lenRS, segLens
    
    def RSpt2LineT2(self, pt1, t1, targLine, t2):
        
        targLine_trans = self.RotateTransLine(Config(pt1[0], pt1[1], t1), targLine)
        finalHead = t2-t1       
        finPos, lengthRS, segLengths = self.RStoLineT2(targLine_trans, finalHead)
        
        if lengthRS:
            finPos_r = utils.RotateVec(finPos, t1 )
            finPos_tr = (finPos_r[0]+pt1[0], finPos_r[1]+pt1[1])            
        else:
            return (None, None), None, (None, None)
       
        return finPos_tr, lengthRS, segLengths
    
    def RStoLineT2(self, lineSeg, finHdng):
        
        lineSeg_rfl = self.LineReflectionXaxis2(lineSeg)
        finPos, lengthsRS, segLengths= self.LStoLineT2(lineSeg_rfl, -finHdng)
        if lengthsRS:
            finPos = du.PtReflectionXaxis(finPos)
        return finPos, lengthsRS, segLengths

    def LocalMinRS(self, lineSeg):

        lineSeg_rfl = self.LineReflectionXaxis2(lineSeg)
        finHdng = self.LocalMinLS(lineSeg_rfl)
        finHdng = np.mod(-finHdng, 2*np.pi)
        return finHdng
    
    def PlotScenario(self, pathType, startConf, segLengths, lineSeg):
        lenLineSeg = np.sqrt((lineSeg.point1[0]-lineSeg.point2[0])**2+(lineSeg.point1[1]-lineSeg.point2[1])**2)
        arrowLen = .2*lenLineSeg
        intrfmt = SimpleNamespace(color='orange', linewidth=1, linestyle='--', marker='x', arrLen=arrowLen)                
        arrfmt = SimpleNamespace(color='cyan', linewidth=1, linestyle='--', marker='x', arrLen=arrowLen)                
        linesegfmt0 = SimpleNamespace(color='g', linewidth=1.5, linestyle='-', marker='x',endPoints=True)               
        utils.PlotLineSeg(self.line1.point1, self.line1.point2, linesegfmt0)        
        utils.PlotLineSeg(self.line2.point1, self.line2.point2, linesegfmt0) 
        utils.PlotInterval(self.line1.point1, self.int1, intrfmt)  
        utils.PlotInterval(self.line2.point1, self.int2, intrfmt)                   
        linesegfmt = SimpleNamespace(color='m', linewidth=1.5, linestyle='-', marker='x',endPoints=True)        
        utils.PlotLineSeg(lineSeg.point1, lineSeg.point2, linesegfmt)        
        pathfmt = SimpleNamespace(color='b', linewidth=2, linestyle='-', marker='x', endPoints=False)    
        finalConf = du.PlotDubPathSegments(startConf, pathType, segLengths, self.rho, pathfmt)
        
        utils.PlotArrow((finalConf[0], finalConf[1]), finalConf[2], arrowLen, arrfmt) 
        plt.axis('equal')                
        return
    
    def LocalMinLSL(self, lineSeg, finHdng):

        A = lineSeg.point1; B = lineSeg.point2
        delX = B[0]-A[0]
        delY = B[1]-A[1]
        rho = self.rho
        lam = (-(A[1]+rho*cos(finHdng)-rho)*delY-(A[0]-rho*sin(finHdng))*delX)/(delX**2+delY**2)
        if lam>=0 and lam<=1:
            finPt = (A[0]+lam*(B[0]-A[0]), A[1]+lam*(B[1]-A[1]) )
            
            pathLSL = dubins.path([0,0,0],[finPt[0], finPt[1],finHdng],rho, 0)   
            length_seg1 = pathLSL.segment_length(0)
            length_seg3 = pathLSL.segment_length(2)
            
            if np.abs(length_seg1-2*np.pi*rho) <1e-8:
                length_seg1 = length_seg1-2*np.pi*rho
            if np.abs(length_seg3/rho-2*np.pi*rho) <1e-8:
                length_seg3 = length_seg3-2*np.pi*rho
            # finHdng = np.mod(length_seg1/rho + length_seg3/rho, 2*np.pi)
            
            lengthLSL = length_seg1+length_seg3+pathLSL.segment_length(1)
        else:
            return None, None, None

        return finPt, lengthLSL, (length_seg1, pathLSL.segment_length(1), length_seg3)
    
    def LocalMinRSR(self, lineSeg, finHdng):
        lineSeg_rfl = self.LineReflectionXaxis2(lineSeg)
        finHdng_rfl = np.mod(-finHdng, 2*np.pi)
        finPt, lengthRSR, segLengths = self.LocalMinLSL(lineSeg_rfl, finHdng_rfl)
        if lengthRSR:
            finPt = du.PtReflectionXaxis(finPt)
        
        return finPt, lengthRSR, segLengths
    
    def LocalMinLSR(self, lineSeg, finHdng):

        rho = self.rho
        
        A = np.array(lineSeg.point1); B = np.array(lineSeg.point2)
        C1 = np.array([0,rho])
        delX = B[0]-A[0]
        delY = B[1]-A[1]

        h = du.DistPtToLineSeg2(C1, lineSeg)
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
        
        for lam in lamRoots:
            if np.imag(lam) == 0 and lam>=0 and lam<=1:
                finPos = A+lam*(B-A)
                pathLSR_cand = dubins.path([0,0,0],[finPos[0], finPos[1],finHdng],rho, 1)
                p1,p2 = du.DubinsInflexionPoints(pathLSR_cand)
                if np.abs(np.dot(p2-p1, B-A)) < 1e-6:
                    finPt = finPos                    
                    lenSeg1 = pathLSR_cand.segment_length(0)
                    lenSeg2 = pathLSR_cand.segment_length(1)
                    lenSeg3 = pathLSR_cand.segment_length(2)
                    
                    if np.abs(lenSeg1-2*rho*np.pi) <1e-8:
                        lenSeg1 = lenSeg1-2*pi*rho
                    if np.abs(lenSeg3-2*np.pi*rho) <1e-8:
                        lenSeg3 = lenSeg3-2*pi*rho
                    
                    return finPt, lenSeg1+lenSeg2+lenSeg3 , (lenSeg1, lenSeg2, lenSeg3)
                    

        return None, None, None
    
    def LocalMinRSL(self, lineSeg, finHdng):
        
        lineSeg_rfl = self.LineReflectionXaxis2(lineSeg)
        finHdng_rfl = np.mod(-finHdng, 2*np.pi)
        finPt, lengthRSL, segLengths = self.LocalMinLSR(lineSeg_rfl, finHdng_rfl)
        
        if lengthRSL:
            finPt = du.PtReflectionXaxis(finPt)

        return finPt, lengthRSL, segLengths

    def LineReflectionXaxis2(self, lineSeg):

        pt1 = (lineSeg.point1[0], -lineSeg.point1[1])
        pt2 = (lineSeg.point2[0], -lineSeg.point2[1])
        
        return LineSegment(pt1, pt2)
    
    def LStoLineT2(self, lineSeg, finHdng):
        # LS path from (0,0,0) to a line segment and given final heading
        rho = self.rho
        A = np.array(lineSeg.point1); B = np.array(lineSeg.point2)        
        finHdng = np.mod(finHdng, 2*np.pi)
        C = np.array([0,rho])
        crossprodCACB = np.cross(A-C, B-C)
        if crossprodCACB > 0:
            B = np.array(lineSeg.point1); A = np.array(lineSeg.point2)        

        lengthsLS = None
        segLengths = (None, None)
        finPos = (None, None)

        lam_den = (B[0]-A[0])*sin(finHdng)-(B[1]-A[1])*cos(finHdng) 
        if np.abs(lam_den) <= 1e-10:
            return finPos, lengthsLS, segLengths
        lam = (A[1]*cos(finHdng)-A[0]*sin(finHdng)+rho-rho*cos(finHdng))/lam_den
        # if utils.InInt(finHdngRange[0], finHdngRange[1], finHdng):
        if lam>=0 and lam<=1:
            finPos = (1-lam)*A + lam*B
            C2 = finPos + rho*np.array([-sin(finHdng), cos(finHdng)])
            ls = np.linalg.norm(C-C2)
            if np.abs((A[1]+lam*(B[1]-A[1])+rho*cos(finHdng)-rho)/ls - sin(finHdng) )<1e-6 and np.abs((A[0]+lam*(B[0]-A[0])-rho*sin(finHdng) )/ls - cos(finHdng))<1e-6:
                lengthsLS = rho*finHdng+ls
                segLengths = (rho*finHdng, ls)
                
            else:
                finPos = (None, None)
        return finPos, lengthsLS, segLengths
    
    def LocalMinLS(self, lineSeg):
        rho = self.rho
        A = np.array(lineSeg.point1); B = np.array(lineSeg.point2)                
        C = np.array([0,rho])
        crossprodCACB = np.cross(A-C, B-C)
        if crossprodCACB > 0:
            B = np.array(lineSeg.point1); A = np.array(lineSeg.point2)        

        gamma = np.arctan2(B[1]-A[1], B[0]-A[0])
        finHdng = np.mod(gamma+np.pi/2, 2*np.pi)
        # finPos, lengthsLS, segLengths = self.LStoLineT2(lineSeg, finHdng)
        return finHdng
    
    def LocalMinLS2(self, lineSeg):
        # local minimum of LS wrt to lam_1, segment S is perpendicular to first line seg
        # This works for RS as well
        psi = np.arctan2(lineSeg.point2[1]-lineSeg.point1[1], lineSeg.point2[0]-lineSeg.point1[0])
        finHdngs = (np.mod(psi+np.pi/2, 2*np.pi), np.mod(psi-np.pi/2, 2*np.pi))
        return finHdngs
    
    def SLtoLine(self, lineSeg, finHdng):
        rho = self.rho
        A = np.array(lineSeg.point1); B = np.array(lineSeg.point2)                
        finHdng = np.mod(finHdng, 2*np.pi)
        C = np.array([0,rho])
        crossprodCACB = np.cross(A-C, B-C)
        if crossprodCACB > 0:
            B = np.array(lineSeg.point1); A = np.array(lineSeg.point2) 
            
        lam = (rho-rho*cos(finHdng)-A[1])/(B[1]-A[1])
        ls = A[0]+lam*(B[0]-A[0])-rho*sin(finHdng)
        if lam>=0 and lam <=1 and ls>0:
            finPos = (1-lam)*A + lam*B        
            lengthSL = ls+rho*finHdng
            segLengths = (ls, rho*finHdng)
        else:
            lengthSL = None
            segLengths = (None, None)            
            finPos = (None, None)
        return finPos, lengthSL, segLengths

    def SRtoLine(self, lineSeg, finHdng):
        rho = self.rho        
        lineSeg_rfl = self.LineReflectionXaxis2(lineSeg)
        finPos, lengthSR, segLengths = self.SLtoLine(lineSeg_rfl, -finHdng)
        if lengthSR:
            finPos = du.PtReflectionXaxis(finPos)

        return finPos, lengthSR, segLengths

    def LocalMinSL(self, lineSeg):
        rho = self.rho        
        A = np.array(lineSeg.point1); B = np.array(lineSeg.point2)                

        crossprodCACB = np.cross(A, B)
        if crossprodCACB > 0:
            B = np.array(lineSeg.point1); A = np.array(lineSeg.point2)        
        gamma = np.arctan2(B[1]-A[1], B[0]-A[0])
        finHdng = np.pi + 2*gamma
        finPos, minlengthSL, segLengths = self.SLtoLine(lineSeg, finHdng, rho)
        return finPos, finHdng, minlengthSL, segLengths

    def LocalMinSR(self, lineSeg):
        rho = self.rho
        lineSeg_rfl = self.LineReflectionXaxis2(lineSeg)
        finPos, finHdng, minLenSR = self.LocalMinSL(lineSeg_rfl, rho)
        
        finPos = du.PtReflectionXaxis(finPos)
        finHdng = np.mod(-finHdng, 2*np.pi)
        return finPos, finHdng, minLenSR
    
    def AddLRPaths(self):
        
        theta_1_u = self.int1[1]   
        theta_2_u = self.int2[1]   
        
        # LR paths from line1_pt1 or line1_pt2 to line2_pt1 or line2_pt2
        listPathsLR = []
        for iniPt in (self.line1.point1, self.line1.point2):
            for finPt in (self.line2.point1, self.line2.point2):
                LR_list = self.LRpt2ptPaths(iniPt, theta_1_u, finPt)
                for pathLR in LR_list:
                    listPathsLR.append((pathLR, iniPt, finPt))
   
        for tup_path_finpt in listPathsLR:  
            
            pathLengths = tup_path_finpt[0]
            iniPt = tup_path_finpt[1]
            finPt = tup_path_finpt[2]
            
            # if pathLengths[2][1]/self.rho > np.pi:       
            self.VerAddTwoSegPath('LR', iniPt, theta_1_u, finPt, pathLengths[0], pathLengths[1], pathLengths[2])        
            # self.PlotScenario('LR', (iniPt[0], iniPt[1], theta_1_u), pathLengths[2], self.line2)      
            # plt.show()

        # LR paths from line1 pt1, theta_1_u to line2, theta_2_l, theta_2_u
        for startPt in (self.line1.point1, self.line1.point2):
            for startHead in self.int1:
                line2_trans = self.RotateTransLine(Config(startPt[0], startPt[1], startHead), self.line2)
                for finalHead in (self.int2[0], self.int2[1]):                
                    list_tup_pathLR = self.LRtoLineT2(line2_trans, finalHead-startHead)
                    for tup_LR in list_tup_pathLR:
                        # if tup_LR[2][1]/self.rho > np.pi:
                        finPt = tup_LR[0]
                        finPos_r = utils.RotateVec(finPt, startHead )
                        finPos_tr = (finPos_r[0]+startPt[0], finPos_r[1]+startPt[1])                    
                        self.VerAddTwoSegPath('LR', startPt, startHead, finPos_tr, finalHead, tup_LR[1], tup_LR[2])        
                        # self.PlotScenario('LR', (startPt[0], startPt[1], startHead), tup_LR[2], self.line2)      
                        # plt.show()

        # LR paths from line1_pt1, line1_pt2 to line2_pt1, line2_pt2
        # Final heading is theta_2_u, and initial heading is free
        # This is done by finding reverse path from final point to start point, with finalHead+pi as initial heading

        listPathsLR = []
        for iniPt in [self.line2.point1, self.line2.point2]:
            for finPt in [self.line1.point1, self.line1.point2]:
                LR_list = self.LRpt2ptPaths(iniPt, theta_2_u+np.pi, finPt)
                for pathLR in LR_list:
                    listPathsLR.append((pathLR, iniPt, finPt))    
        for tup_path_finpt in listPathsLR:  
            
            pathLengths = tup_path_finpt[0]
            iniPt = tup_path_finpt[1]
            finPt = tup_path_finpt[2]
            
            # if pathLengths[2][1]/self.rho > np.pi:       
            self.VerAddTwoSegPath('LR', finPt, pathLengths[0]+np.pi, iniPt, theta_2_u , pathLengths[1], pathLengths[2][::-1])        
            # self.PlotScenario('LR', (finPt[0], finPt[1], pathLengths[0]+np.pi), pathLengths[2][::-1], self.line2)      
            # plt.show()

        # LR paths from line1, theta_1_e to line2_pt1 or line2_pt2, theta_2_u
        for startPt in (self.line2.point1, self.line2.point2):
            line1_trans = self.RotateTransLine(Config(startPt[0], startPt[1], theta_2_u+np.pi), self.line1)
            for finalHead in (self.int1[0]+np.pi, self.int1[1]+np.pi):                
                list_tup_pathLR = self.LRtoLineT2(line1_trans, finalHead-theta_2_u-np.pi)
                for tup_LR in list_tup_pathLR:
                    # if tup_LR[2][1]/self.rho > np.pi:
                    finPt = tup_LR[0]
                    finPos_r = utils.RotateVec(finPt, theta_2_u+np.pi )
                    finPos_tr = (finPos_r[0]+startPt[0], finPos_r[1]+startPt[1])                    
                    self.VerAddTwoSegPath('LR', finPos_tr, finalHead-np.pi, startPt, theta_2_u, tup_LR[1], tup_LR[2][::-1])        
                    # self.PlotScenario('LR', (finPos_tr[0], finPos_tr[1], finalHead-np.pi), tup_LR[2][::-1], self.line2)      
                    # plt.show()
                
        return
    
    def LRpt2ptPaths(self, pt1, t1, pt2):
             
        # LR from pt1, t1 to pt2
        finHead = None
        pt2_tans = (pt2[0]-pt1[0], pt2[1]-pt1[1])
        pt2_transrot = utils.RotateVec(pt2_tans, -t1 )

        lengthsLR1, lengthsLR2 = du.PathLR(pt2_transrot[0], pt2_transrot[1], self.rho)        
        
        lengthsList = []
        for segLengths in (lengthsLR1, lengthsLR2 ):
            if segLengths[0]:
                finHead = segLengths[0]/self.rho - segLengths[1]/self.rho + t1
                lengthsList.append((finHead, segLengths[0]+segLengths[1], segLengths))
       
        return lengthsList
    
    def LRtoLineT2(self, lineSeg, finHdng):
    # Computes path from (0,0,0) to any point on the lineSeg (AB), with a given final heading (finHdng)
        rho = self.rho
        A = np.array(lineSeg.point1); B = np.array(lineSeg.point2)                
        h = du.DistPtToLineSeg2(np.array([0,rho]), lineSeg)
        if h>3*rho:
            return []
        else:
            delX = B[0]-A[0]
            delY = B[1]-A[1]
            T1 = A[0]+rho*sin(finHdng)
            T2 = A[1]-rho*cos(finHdng)-rho
            quadEqnCoeffs = np.array([delX**2 + delY**2, 2*T1*delX+2*T2*delY, T1**2+T2**2-4*rho*rho])
            lams = np.roots(quadEqnCoeffs)

            returnList = []
            for k in range(0,2):
                if np.imag(lams[k]) == 0 and lams[k]>=0 and lams[k]<=1:
                    finPt = A + lams[k]*(B-A)
                    lenLR, segLengths = self.LengthOfLR(np.array([0,0,0]), np.array([finPt[0], finPt[1], finHdng]))
                    if lenLR:
                        returnList.append((finPt, lenLR, segLengths))
            
        return returnList
    
    def LengthOfLR(self, iniConf, finConf):
        
        rho = self.rho
        C2 = finConf[0:2] + rho*np.array([sin(finConf[2]), -cos(finConf[2])])
        C1 = iniConf[0:2] + rho*np.array([-sin(iniConf[2]), cos(iniConf[2])])

        if np.abs(np.linalg.norm(C1-C2)-2*rho)<0.001:
            gamma = np.arctan2(C2[1]-C1[1], C2[0]-C1[0])
            phi1 = np.mod(pi/2 + gamma, 2*pi)
            phi2 = np.mod(phi1-finConf[2], 2*pi)
            lenLR = rho*(phi1+phi2)
            
            return lenLR, (rho*phi1, rho*phi2)
        else:
            return None, (None, None)
    
    def AddRLPaths(self):
        
        line1_rfl = self.LineReflectionXaxis2(self.line1)
        line2_rfl = self.LineReflectionXaxis2(self.line2)
        int1_rfl = (-self.int1[1], -self.int1[0])
        int2_rfl = (-self.int2[1], -self.int2[0])
        
        
        L2LDubRfl = Line2LineDubins(line1_rfl, int1_rfl, line2_rfl, int2_rfl, self.rho) 
        L2LDubRfl.AddLRPaths()
        for cp_r in L2LDubRfl.candPathsList:            
            cp = CandidatePath('RL', du.PtReflectionXaxis(cp_r.iniPos) , np.mod(-cp_r.iniHead, 2*np.pi), du.PtReflectionXaxis(cp_r.finalPos), np.mod(-cp_r.finalHead, 2*np.pi), cp_r.segLengths)
            self.candPathsList.append(cp)
            self.lengthsVec.append(cp_r.segLengths[0]+cp_r.segLengths[1])  
            # self.PlotScenario('RL', (cp_r.iniPos[0], -cp_r.iniPos[1], -cp_r.iniHead), cp_r.segLengths, self.line2)      
            # plt.show()
        return
    
    def AddLPaths(self):
        
        # given start points be line1 end points
        # final points to be line2 end points
        for startPt in (self.line1.point1, self.line1.point2):
            for finalPt in (self.line2.point1, self.line2.point2):
                # iniHead, finalHead, lenL = self.PathL(startPt, finalPt)
                pathL_cands = self.PathL(startPt, finalPt)                
                for pathL_params in pathL_cands:
                    iniHead = pathL_params[0]
                    finalHead = pathL_params[1]
                    lenL = pathL_params[2]
                    if lenL:
                        # self.PlotScenario('LR', (startPt[0], startPt[1], iniHead), (lenL,0), self.line2) 
                        self.VerAddTwoSegPath('L', startPt, iniHead , finalPt, finalHead, lenL, (lenL,0))
        
        # given start position be line1 end points
        # initial heading to be interval 1 boundaries
        for startPt in [self.line1.point1, self.line1.point2]:
            for iniHead in self.int1:
                
                line2_trans = self.RotateTransLine(Config(startPt[0], startPt[1], iniHead), self.line2)
                returnList = self.LtoLine(line2_trans)
        
                for pathL_vals in returnList:
                    finPos_r = utils.RotateVec(pathL_vals[0], iniHead)
                    finPos_tr = (finPos_r[0]+startPt[0], finPos_r[1]+startPt[1])  
                    finHead = pathL_vals[1]+iniHead
                    self.VerAddTwoSegPath('L', startPt, iniHead , finPos_tr, finHead, pathL_vals[2], (pathL_vals[2],0))
                    # self.PlotScenario('LR', (startPt[0], startPt[1], iniHead), (pathL_vals[2],0), self.line2)      
               
        # given final position be line2 end points
        # initial heading to be int1 boundaries                
        for endPt in [self.line2.point1, self.line2.point2]:        
            for t1 in self.int1:
                
                p1_rev = endPt
                t2_rev = t1+np.pi
                pathsList = self.PathRtoLineT2(p1_rev, self.line1, t2_rev)
                
                for pathR_params in pathsList:
                # t1_rev, p2_rev, lenL = self.PathRtoLineT2(p1_rev, self.line1, t2_rev)
                # if lenL:
                    t1_rev = pathR_params[0]
                    p2_rev = pathR_params[1]
                    lenL = pathR_params[2]
                
                    t2 = np.mod(t1_rev+np.pi, 2*np.pi)
                    p1 = p2_rev
                    self.VerAddTwoSegPath('L', (p1[0], p1[1]), t1, endPt, t2, lenL, (lenL,0))                    
                    # self.PlotScenario('LR', (p1[0], p1[1], t1), (lenL,0), self.line2)      

        # given start position be line1 end points
        # final heading to be int2 boundaries                
        for startPt in [self.line1.point1, self.line1.point2]:        
            for t2 in self.int2:
                pathsL_list = self.PathLtoLineT2(startPt, self.line2, t2)
                
                # iniHead, p2, lenL = self.PathLtoLineT2(startPt, self.line2, t2)
                for pathL_params in pathsL_list:
                    iniHead = pathL_params[0]
                    p2 = pathL_params[1]
                    lenL = pathL_params[2]
                # if lenL:
                    self.VerAddTwoSegPath('L', startPt, iniHead , p2, t2, lenL, (lenL,0))                    
                    # self.PlotScenario('LR', (startPt[0], startPt[1], iniHead), (lenL,0), self.line2)      

        # given end position to be line2 end points
        # final heading to be int2 boundaries                
        for endPt in [self.line2.point1, self.line2.point2]:        
            for t2 in self.int2:
                
                p1_rev = endPt
                t1_rev = t2+np.pi                
                line1_trans = self.RotateTransLine(Config(p1_rev[0], p1_rev[1], t1_rev), self.line1)                
                returnList = self.RtoLine(line1_trans)
                
                for path_cand in returnList:
                    p2_rev = path_cand[0]
                    t2_rev = path_cand[1]
                    lenL = path_cand[2]
                    
                    p1_r = utils.RotateVec(p2_rev, t1_rev )
                    startPos = (p1_r[0]+p1_rev[0], p1_r[1]+p1_rev[1])  
                    t1 = np.mod(t2_rev+t1_rev+np.pi, 2*np.pi)
                    if lenL:
                        self.VerAddTwoSegPath('L', (startPos[0], startPos[1]), t1 ,(endPt[0], endPt[1]), t2, lenL, (lenL,0))                    
                        # self.PlotScenario('LR', (startPos[0], startPos[1], t1), (lenL,0), self.line2)      
      
                                 
        # given initial heading to be int1 boundaries
        # final heading be int2 boundaries        
        for iniHead in self.int1:
            for finHead in self.int2:
                p1, p2, lenL = self.PathLT1toT2(self.line1, self.line2, iniHead, finHead)
                # self.PlotScenario('LR', (self.line1.point1[0], self.line1.point1[1], iniHead), (0,0), self.line2)      
                # self.PlotScenario('LR', (self.line2.point1[0], self.line2.point1[1], finHead), (0,0), self.line2)      
                
                if lenL:
                    self.VerAddTwoSegPath('L', (p1[0], p1[1]), iniHead , (p2[0], p2[1]), finHead, lenL, (lenL,0))                    
                    # self.PlotScenario('LR', (p1[0], p1[1], iniHead), (lenL,0), self.line2)      
                    
        
        return
    
    def PathLT1toT2(self, line1, line2, t1, t2):
        
        iniPos = None
        finalPos = None
        lenL = None
        alpha = np.mod(t2-t1, 2*np.pi) # arc angle
        startConf = [line1.point1[0], line1.point1[1], t1]
        finalConf = du.MoveAlongSeg(startConf, alpha*self.rho, 'L', self.rho)
        pf = np.array([finalConf[0], finalConf[1]])
        v1 = np.array(line1.point2)-np.array(line1.point1)        
        v2 = np.array(line2.point2)-np.array(line2.point1)
        A = np.array(line1.point1)        
        C = np.array(line2.point1)
        vmat = np.transpose(np.array([v1, -v2]))
        if np.linalg.det(vmat) != 0:
            lams = np.linalg.solve(vmat, C-pf)
            if lams[0]>=0 and lams[0]<=1:
                if lams[1]>=0 and lams[1]<=1:
                    iniPos = A+lams[0]*v1
                    finalPos = C+lams[1]*v2
                    lenL = alpha*self.rho
            
        return iniPos, finalPos, lenL
    
    def LtoLine(self, lineSeg):    
        rho=self.rho        
        center = np.array([0, rho])
        h = du.DistPtToLineSeg2(center, lineSeg)
        returnList = []
        if h<=rho:
            finPosList = du.IntersectionLineCircle2(lineSeg, [0,rho], rho)
            
            for finPos in finPosList:            
                # if du.CheckPtLiesOnLineSeg(finPos, lineSeg):
                finHdng = np.arctan2(finPos[1]-rho, finPos[0])+pi/2
                returnList.append(( finPos, finHdng, rho*np.mod(finHdng, 2*np.pi) ))
                # finHdngList.append(finHdng)
                # finLengthsList.append(rho*np.mod(finHdng, 2*np.pi))
                    
        return returnList

    def RtoLine(self, lineSeg):

        # assumes initial configuration is [0,0,0]
        lineSeg_rfl = self.LineReflectionXaxis2(lineSeg)
        returnListL = self.LtoLine(lineSeg_rfl)
        
        # finPosList_rfl, finHdngList_rfl, finLengthsList = self.LtoLine(lineSeg_rfl)
        # finPosList = []
        # finHdngList = []
        returnList = []
        for revPathL in returnListL:
            finPos = du.PtReflectionXaxis(revPathL[0])
            finHdng =  np.mod(-revPathL[1], 2*np.pi )            
            returnList.append( (finPos, finHdng, revPathL[2]) )
            
        return returnList

    def PathL(self, pt1, pt2):    
        rho=self.rho        
        d = np.sqrt( (pt1[0]-pt2[0])**2 + (pt1[1]-pt2[1])**2)
        if d>2*rho:
            return [(None, None, None)]
        psi = np.arctan2(pt2[1]-pt1[1], pt2[0]-pt1[0])
        alpha = np.arccos(d/2/rho)
        returnList = [] #each entry is a tuple of t1, t2, length
        
        t1_v = [psi-alpha-np.pi/2, psi+alpha-np.pi/2] 
        phi_v = [np.pi+2*alpha, np.pi-2*alpha]
        for indx, t1 in enumerate(t1_v):
            phi = phi_v[indx]
            t2 = np.mod(t1+phi, 2*np.pi)
            returnList.append((t1, t2, rho*np.mod(phi, 2*np.pi)))
            
        return returnList
    
    def PathR(self, pt1, pt2):    
        rho=self.rho        
        d = np.sqrt( (pt1[0]-pt2[0])**2 + (pt1[1]-pt2[1])**2)
        if d>2*rho:
            return [(None, None, None)]
        psi = np.arctan2(pt2[1]-pt1[1], pt2[0]-pt1[0])
        alpha = np.arccos(d/2/rho)
        
        returnList = [] #each entry is a tuple of t1, t2, length
        t1_v = [psi-alpha+np.pi/2, psi+alpha+np.pi/2]
        phi_v = [np.pi - 2*alpha, np.pi + 2*alpha]
        for indx, t1 in enumerate(t1_v):                  
            phi = phi_v[indx]
            t2 = np.mod(t1-phi, 2*np.pi)
            returnList.append((t1, t2, rho*np.mod(phi, 2*np.pi)))
        return returnList
    
    def PathLtoLineT2(self, p1, lineSeg, t2):
        # initial point, final line, and final heading are given
        # initial heading is free
        t2 = np.mod(t2, 2*np.pi)
        rho=self.rho    
        A = np.array(lineSeg.point1); B = np.array(lineSeg.point2)                            
        xa = lineSeg.point1[0]; ya = lineSeg.point1[1]        
        vx = lineSeg.point2[0]-xa; vy = lineSeg.point2[1]-ya
        
        xbar = xa-rho*np.sin(t2)-p1[0]
        ybar = ya+rho*np.cos(t2)-p1[1]
        
        quadEqCoeffs = np.array([vx**2+vy**2, 2*vx*xbar+2*vy*ybar, xbar**2+ybar**2-rho**2])
        lamRoots = np.roots(quadEqCoeffs)
        returnList = []        
        for lam in lamRoots:
            if np.imag(lam) == 0 and lam>=0 and lam<=1:
                p2 = A + lam*(B-A)
                # iniHead, finalHead, lenL = self.PathL(p1, p2)
                
                pathL_cands = self.PathL(p1, p2)
                
                for pathL_params in pathL_cands:
                    iniHead = pathL_params[0]
                    finalHead = pathL_params[1]
                    lenL = pathL_params[2]
                    if np.abs(finalHead-t2) < 0.001:
                        returnList.append((iniHead, p2, lenL))                        
                        # return iniHead, p2, lenL
        
        return returnList
    
    def PathRtoLineT2(self, p1, lineSeg, t2):
        # initial point, final line, and final heading are given
        # initial heading is free
        # finds path R
        t2 = np.mod(t2, 2*np.pi)
        rho=self.rho    
        A = np.array(lineSeg.point1); B = np.array(lineSeg.point2)                            
        xa = lineSeg.point1[0]; ya = lineSeg.point1[1]        
        vx = lineSeg.point2[0]-xa; vy = lineSeg.point2[1]-ya
        
        xbar = xa+rho*np.sin(t2)-p1[0]
        ybar = ya-rho*np.cos(t2)-p1[1]
        
        quadEqCoeffs = np.array([vx**2+vy**2, 2*vx*xbar+2*vy*ybar, xbar**2+ybar**2-rho**2])
        lamRoots = np.roots(quadEqCoeffs)
        
        returnList = []
        for lam in lamRoots:
            if np.imag(lam) == 0 and lam>=0 and lam<=1:
                p2 = A + lam*(B-A)
                # iniHead, finalHead, lenR = self.PathR(p1, p2)
                pathL_cands = self.PathR(p1, p2)
                
                for pathL_params in pathL_cands:
                    iniHead = pathL_params[0]
                    finalHead = pathL_params[1]
                    lenR = pathL_params[2]
                    if np.abs(finalHead-t2) < 0.001:
                        returnList.append((iniHead, p2, lenR))
                        # return iniHead, p2, lenR
        
        return returnList
    
    def AddRPaths(self):
        line1_rfl = self.LineReflectionXaxis2(self.line1)
        line2_rfl = self.LineReflectionXaxis2(self.line2)
        int1_rfl = (-self.int1[1], -self.int1[0])
        int2_rfl = (-self.int2[1], -self.int2[0])
        
        
        L2LDubRfl = Line2LineDubins(line1_rfl, int1_rfl, line2_rfl, int2_rfl, self.rho) 
        L2LDubRfl.AddLPaths()
        for cp_r in L2LDubRfl.candPathsList:            
            # cp = CandidatePath('R', du.PtReflectionXaxis(cp_r.iniPos) , np.mod(-cp_r.iniHead, 2*np.pi), du.PtReflectionXaxis(cp_r.finalPos), np.mod(-cp_r.finalHead, 2*np.pi), cp_r.segLengths)
            self.VerAddTwoSegPath('R', du.PtReflectionXaxis(cp_r.iniPos), np.mod(-cp_r.iniHead, 2*np.pi) , du.PtReflectionXaxis(cp_r.finalPos), np.mod(-cp_r.finalHead, 2*np.pi), cp_r.segLengths[0], cp_r.segLengths)                    
            
            # self.candPathsList.append(cp)
            # self.lengthsVec.append(cp_r.segLengths[0])  
            # self.PlotScenario('R', (cp_r.iniPos[0], -cp_r.iniPos[1], -cp_r.iniHead), cp_r.segLengths, self.line2)      
            # plt.show()
        return
    
    def AddSPaths(self):
        
        # straight line paths joining the end points
        # headings = [] 
        # first two of the headings list are headings of line1_pt1 to line2_pt1 and line_pt2
        # third and fourth are  headings from line1_pt2 to line2_pt1 and line_pt2
        intPt = utils.IntersectionLineSegments2(self.line1, self.line2)
        if intPt[0] and  utils.CheckIntrvsIntersect(self.int1, self.int2):
            commonHead = utils.CommonAngle(self.int1, self.int2)
            self.VerAddTwoSegPath('S', intPt, commonHead , intPt, commonHead, 0, (0,0))            
            return
        
        for startPt in (self.line1.point1, self.line1.point2):
            for finalPt in (self.line2.point1, self.line2.point2):
                iniHead = np.arctan2(finalPt[1]-startPt[1], finalPt[0]-startPt[0])
                # headings.append(iniHead)
                lenS = np.sqrt((finalPt[1]-startPt[1])**2+(finalPt[0]-startPt[0])**2)
                # self.PlotScenario('S', (startPt[0], startPt[1], iniHead), (lenS,0), self.line2) 
                self.VerAddTwoSegPath('S', startPt, iniHead , finalPt, iniHead, lenS, (lenS,0))
        
        # straight line points from end points of first line to interval boundaries
        
        intrBoundaries = [self.int1[0], self.int1[1], self.int2[0], self.int2[1]]
        pathSList = self.AddSPathP1toT2(self.line1, self.line2, intrBoundaries)
        for pathL_params in pathSList:
            startPt = pathL_params[0]
            hdng = pathL_params[1]
            finalPt = pathL_params[2]
            lenS = pathL_params[3]
            # self.PlotScenario('S', (startPt[0], startPt[1], hdng), (lenS,0), self.line2) 
            self.VerAddTwoSegPath('S', startPt, hdng , finalPt, hdng, lenS, (lenS,0))

        
        intrBoundaries = [self.int1[0]+np.pi, self.int1[1]+np.pi, self.int2[0]+np.pi, self.int2[1]+np.pi]
        pathSList = self.AddSPathP1toT2(self.line2, self.line1, intrBoundaries)
        for pathL_params in pathSList:
            startPt = pathL_params[0]
            hdng = pathL_params[1]
            finalPt = pathL_params[2]
            lenS = pathL_params[3]
            # self.PlotScenario('S', (finalPt[0], finalPt[1], hdng+np.pi), (lenS,0), self.line2) 
            self.VerAddTwoSegPath('S', finalPt, hdng+np.pi , startPt, hdng+np.pi, lenS, (lenS,0))
        return
        
    def AddSPathP1toT2(self, lineSeg1, lineSeg2, intrBoundaries):
        
        headings = [] 
        for startPt in (lineSeg1.point1, lineSeg1.point2):
            for finalPt in (lineSeg2.point1, lineSeg2.point2):
                hdng = np.arctan2(finalPt[1]-startPt[1], finalPt[0]-startPt[0])
                headings.append(hdng)
        feasIntrs = [(headings[0], headings[1]), (headings[2], headings[3])]
        returnList = [] 
        for ind, startPt in enumerate((lineSeg1.point1, lineSeg1.point2)):
            feasInt = feasIntrs[ind]            
            p1 = np.array(startPt)
            p3 = np.array(lineSeg2.point1)
            p4 = np.array(lineSeg2.point2)
            crossprodP1P3P4 = np.cross(p3-p1, p4-p1)
            localMinHead = np.arctan2(p4[1]-p3[1], p4[0]-p3[0])-np.pi/2            
            if crossprodP1P3P4<0:
                feasInt = (feasInt[1], feasInt[0])
                localMinHead = localMinHead+np.pi
            if utils.InInt(feasInt[0], feasInt[1], localMinHead):
                    lenS = du.DistPtToLineSeg2(startPt, lineSeg2)
                    # self.PlotScenario('S', (startPt[0], startPt[1], localMinHead), (lenS,0), lineSeg2) 
                    # self.VerAddTwoSegPath('S', startPt, localMinHead , (0,0), localMinHead, lenS, (lenS,0))
                    finalPt = (startPt[0]+lenS*np.cos(localMinHead), startPt[1]+lenS*np.sin(localMinHead))
                    returnList.append((startPt, localMinHead, finalPt, lenS))
            
            for heading in intrBoundaries:
                if utils.InInt(feasInt[0], feasInt[1], heading):
                    lenS = du.DistPtHdngToLineSeg(startPt, heading, lineSeg2)
                    # self.PlotScenario('S', (startPt[0], startPt[1], heading), (lenS,0), lineSeg2) 
                    # self.VerAddTwoSegPath('S', startPt, heading , (0,0), heading, lenS, (lenS,0))
                    finalPt = (startPt[0]+lenS*np.cos(heading), startPt[1]+lenS*np.sin(heading))
                    returnList.append((startPt, heading, finalPt, lenS))
      
        return returnList
            
            
                