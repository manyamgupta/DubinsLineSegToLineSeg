from timeit import default_timer as timer
from types import SimpleNamespace
import utils
import dubutils as du
import DubinsL2L as dl2l
import matplotlib.pyplot as plt
    

line1 = dl2l.LineSegment((8.,30.), (35.,5.)) # arguments are two tuples, each tuple is an end point of the line segment
# line2 = dl2l.LineSegment((11.,-1.),(9.,6.))
# line2 = dl2l.LineSegment((7.,5.),(3.,9.))
# line2 = dl2l.LineSegment((0.,5.), (6.,5.5))
line2 = dl2l.LineSegment((35., 30.), (50., 8))
# line2 = dl2l.LineSegment((0,2.5), (5.5,5))
# int1 = (-.5, .5)
int1 = (-.4, .6)
# int2 = (1.2, 2.)
int2 = (-1.2, .5)

rho = 8

# line1=dl2l.LineSegment(point1=(94, 54), point2=(18, 29))
# line2=dl2l.LineSegment(point1=(18, 33), point2=(67, 41))
# int1=(4.5630963916040805, 5.348494555001529)
# int2=(2.364600408307861, 3.1499985717053094)

line1=dl2l.LineSegment(point1=(94, -54), point2=(18, -29))
line2=dl2l.LineSegment(point1=(18, -33), point2=(67, -41))
int1=(-5.348494555001529, -4.5630963916040805)
int2=(-3.1499985717053094, -2.364600408307861)

## computing shortest arc to arc dubins path using analytical results
tic = timer()
L2LDub = dl2l.Line2LineDubins(line1, int1, line2, int2, rho) 
minLength, minPath = L2LDub.MinDub_L2L()   #This returns the length of the minimum path, minimum path, and all the candidate paths
comp_time = timer()-tic
print(f"{minLength=}")
print("minPathType=", minPath.pathType)
print("minConfStart= ", minPath.iniPos, " ", minPath.iniHead)
print("minConfGoal= ", minPath.finalPos, " ", minPath.finalHead)
print(f"{comp_time=}")

plt.figure()
linesegfmt = SimpleNamespace(color='m', linewidth=1.5, linestyle='-', marker='x',endPoints=True)        
pathfmt = SimpleNamespace(color='b', linewidth=2, linestyle='-', marker='x', endPoints=False)    
arrfmt = SimpleNamespace(color='g', linewidth=1, linestyle='--', marker='x', arrLen=10)       
utils.PlotLineSeg(line1.point1, line1.point2, linesegfmt)
utils.PlotLineSeg(line2.point1, line2.point2, linesegfmt)
utils.PlotInterval(line2.point1, int2, arrfmt)
utils.PlotInterval(line1.point1, int1, arrfmt)
if minLength:   
    arrfmt = SimpleNamespace(color='c', linewidth=1, linestyle='--', marker='x', arrLen=10)           
    du.PlotDubPathSegments((minPath.iniPos[0], minPath.iniPos[1], minPath.iniHead), minPath.pathType, minPath.segLengths, rho, pathfmt)     
    utils.PlotArrow(minPath.finalPos, minPath.finalHead, 10, arrfmt)    
    utils.PlotInterval(minPath.iniPos, int1, arrfmt)
    
plt.axis('equal')
plt.show()
