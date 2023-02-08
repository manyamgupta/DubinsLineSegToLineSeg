# Author: Satyanarayana Gupta Manyam
# Udpate Date: Feb 8th, 2023
# This example script shows how to use the Dubins LineSegment-Interval to LineSegment-Interval code
# DubinsL2l.py has the classes and most of the analytic functions
# dubutils and utils has the utility functions



from timeit import default_timer as timer
from types import SimpleNamespace
import utils
import dubutils as du
import DubinsL2L as dl2l
import matplotlib.pyplot as plt
import numpy as np    

rho = 8

# line1=dl2l.LineSegment(point1=(94, -54), point2=(18, -29))
# line2=dl2l.LineSegment(point1=(18, -33), point2=(67, -41))
# int1=(-5.348494555001529, -4.5630963916040805)
# int2=(-3.1499985717053094, -2.364600408307861)

pts = np.random.randint(0,80,8)
t_1_l = np.random.rand()*np.pi*2
t_2_l = np.random.rand()*np.pi*2
line1 = dl2l.LineSegment((pts[0], pts[1]), (pts[2], pts[3])) 
line2 = dl2l.LineSegment((pts[4], pts[5]), (pts[6], pts[7]))
int1 = (t_1_l, t_1_l+np.pi/4)
int2 = (t_2_l, t_2_l+np.pi/4)

plotFlag  = True
## computing shortest line to line dubins path using analytical results
tic = timer()
L2LDub = dl2l.Line2LineDubins(line1, int1, line2, int2, rho) 
minLength, minPath = L2LDub.MinDub_L2L()   #This returns the length of the minimum path, minimum path, and all the candidate paths
comp_time = timer()-tic
print(f"{minLength=}")
print("minPathType=", minPath.pathType)
print("minConfStart= ", minPath.iniPos, " ", minPath.iniHead)
print("minConfGoal= ", minPath.finalPos, " ", minPath.finalHead)
print("Segment lengths= ", minPath.segLengths)
print(f"{comp_time=}")

if plotFlag:
    plt.figure()
    linesegfmt = SimpleNamespace(color='m', linewidth=1.5, linestyle='-', marker='x',endPoints=True)        
    pathfmt = SimpleNamespace(color='b', linewidth=2, linestyle='-', marker='x', endPoints=False)    
    arrfmt = SimpleNamespace(color='g', linewidth=1, linestyle='--', marker='x', arrLen=10)       
    utils.PlotLineSeg(line1.point1, line1.point2, linesegfmt)
    utils.PlotLineSeg(line2.point1, line2.point2, linesegfmt)
    if minLength:   
        arrfmt = SimpleNamespace(color='g', linewidth=1, linestyle='--', marker='x', arrLen=10)                   
        utils.PlotInterval(minPath.iniPos, int1, arrfmt)
        utils.PlotInterval(minPath.finalPos, int2, arrfmt)
        du.PlotDubPathSegments((minPath.iniPos[0], minPath.iniPos[1], minPath.iniHead), minPath.pathType, minPath.segLengths, rho, pathfmt)     
        arrfmt = SimpleNamespace(color='c', linewidth=1, linestyle='--', marker='x', arrLen=10)           
        utils.PlotArrow(minPath.finalPos, minPath.finalHead, 10, arrfmt) 
        utils.PlotArrow(minPath.iniPos, minPath.iniHead, 10, arrfmt)            

        
    plt.axis('equal')
    plt.show()
