from timeit import default_timer as timer
from types import SimpleNamespace
import utils
import dubutils as du
import DubinsL2L as dl2l
import matplotlib.pyplot as plt
import numpy as np
import DubinsLineToLine_V2 as dl2lv2
import utils_d2p as ut2    
import dubutils_d2p as du2
from collections import namedtuple
import DIP


# line1 = dl2l.LineSegment((1.,4.), (1.,.5)) # arguments are two tuples, each tuple is an end point of the line segment
# line2 = dl2l.LineSegment((2.5,0.), (2.5,-4.))
# int1 = (1.1, 2.5)
# int2 = (-np.pi/4+.1, np.pi/32)

for i in range(100):
    
    pts = np.random.randint(0,30,8)
    t_1_l = np.random.rand()*np.pi*2
    t_2_l = np.random.rand()*np.pi*2
    line1 = dl2l.LineSegment((pts[0], pts[1]), (pts[2], pts[3])) 
    line2 = dl2l.LineSegment((pts[4], pts[5]), (pts[6], pts[7]))
    int1 = (t_1_l, t_1_l+np.pi/4)
    int2 = (t_2_l, t_2_l+np.pi/4)
    rho = 8 #minimum turn radius
    
    # line1=dl2l.LineSegment(point1=(19, 39), point2=(22, 87))
    # line2=dl2l.LineSegment(point1=(49, 67), point2=(9, 85))
    # int1=(1.310365760358409, 2.0957639237558574)
    # int2=(2.071128811456242, 2.8565269748536903)
    
    print("\n************* Instance parameters ************\n")  
    print("instance: ", i)  
    print(f"{line1=}")
    print(f"{line2=}")
    print(f"{int1=}")
    print(f"{int2=}")

    ## computing shortest line to line dubins path using analytical method
    tic = timer()
    L2LDub = dl2l.Line2LineDubins(line1, int1, line2, int2, rho) 
    minLength_dl2l, minPath = L2LDub.MinDub_L2L()   #This returns the length of the minimum path, minimum path, and all the candidate paths
    comp_time = timer()-tic
    print("\n************* Results from Dubins line to line ************\n")
    print(f"{minLength_dl2l=}")
    print("minPathType=", minPath.pathType)
    print("minConfStart= ", minPath.iniPos, " ", minPath.iniHead)
    print("minConfGoal= ", minPath.finalPos, " ", minPath.finalHead)
    print(f"{comp_time=}")

    ## computing shortest line to line dubins path using numerical method
    tic = timer()
    numDisc = 400
    minLength_dl2l_num, minPath_dl2l_num = DIP.DubL2LNumSoln(line1, int1, line2, int2, rho, numDisc)
    comp_time = timer()-tic
    print("\n************* Results from Numerical Dubins line to line ************\n")
    print(f"{minLength_dl2l_num=}")
    print("minPathType=", minPath_dl2l_num.pathType)
    print("minConfStart= ", minPath_dl2l_num.iniPos, " ", minPath_dl2l_num.iniHead)
    print("minConfGoal= ", minPath_dl2l_num.finalPos, " ", minPath_dl2l_num.finalHead)
    print(f"{comp_time=}")
    
    
    #### Computing Dubins L2L with paralellogram code
    # pt1 = line1.point1
    # pt2 = line1.point2
    # pt3 = line2.point1
    # pt4 = line2.point2

    # line1_prlgrm = [pt1, pt2]
    # line2_prlgrm = [pt3, pt4]
    # sector1 = [int1[0], int1[1]]
    # sector2 = [int2[0], int2[1]]
    # tic = timer()
    # minLength_dp2prlgrm, minConfStart, minConfGoal, minPathType2, minPathSegLengths = dl2lv2.DubinsLineToLineV2(line1_prlgrm, sector1, line2_prlgrm, sector2, rho)
    # comp_time = timer()-tic

    # print("\n************* Results from Dubins to Paralellogram ************\n")
    # print(f"{minLength_dp2prlgrm=}")
    # print(f"{minConfStart=}")
    # print(f"{minConfGoal=}")
    # print(f"{minPathType2=}")
    # print(f"{comp_time=}")

    if minLength_dl2l-minLength_dl2l_num >.01 :
    # if True:
        ########### Plot from Dubins line to line #############
        plt.figure() 
        arrowLength = 10
        linesegfmt = SimpleNamespace(color='m', linewidth=1.5, linestyle='-', marker='x',endPoints=True)        
        pathfmt = SimpleNamespace(color='b', linewidth=2, linestyle='-', marker='x', endPoints=False)    
        arrfmt = SimpleNamespace(color='g', linewidth=1, linestyle='--', marker='x', arrLen=arrowLength)       
        utils.PlotLineSeg(line1.point1, line1.point2, linesegfmt)
        utils.PlotLineSeg(line2.point1, line2.point2, linesegfmt)
        utils.PlotInterval(line2.point1, int2, arrfmt)
        utils.PlotInterval(line1.point1, int1, arrfmt)
        if minLength_dl2l:   
            arrfmt = SimpleNamespace(color='c', linewidth=1, linestyle='--', marker='x', arrLen=arrowLength)           
            du.PlotDubPathSegments((minPath.iniPos[0], minPath.iniPos[1], minPath.iniHead), minPath.pathType, minPath.segLengths, rho, pathfmt)     
            utils.PlotArrow(minPath.finalPos, minPath.finalHead, arrowLength, arrfmt)    
            utils.PlotArrow(minPath.iniPos, minPath.iniHead, arrowLength, arrfmt)                
            arrfmt = SimpleNamespace(color='g', linewidth=1, linestyle='--', marker='x', arrLen=arrowLength)                   
            utils.PlotInterval(minPath.iniPos, int1, arrfmt)
            
        plt.axis('equal')

        ########### Plot from Dubins line to line numerical #############
        plt.figure() 
        arrowLength = 10
        linesegfmt = SimpleNamespace(color='m', linewidth=1.5, linestyle='-', marker='x',endPoints=True)        
        pathfmt = SimpleNamespace(color='b', linewidth=2, linestyle='-', marker='x', endPoints=False)    
        arrfmt = SimpleNamespace(color='g', linewidth=1, linestyle='--', marker='x', arrLen=arrowLength)       
        utils.PlotLineSeg(line1.point1, line1.point2, linesegfmt)
        utils.PlotLineSeg(line2.point1, line2.point2, linesegfmt)
        utils.PlotInterval(line2.point1, int2, arrfmt)
        utils.PlotInterval(line1.point1, int1, arrfmt)
        if minLength_dl2l_num:   
            arrfmt = SimpleNamespace(color='c', linewidth=1, linestyle='--', marker='x', arrLen=arrowLength)           
            du.PlotDubPathSegments((minPath_dl2l_num.iniPos[0], minPath_dl2l_num.iniPos[1], minPath_dl2l_num.iniHead), minPath_dl2l_num.pathType, minPath_dl2l_num.segLengths, rho, pathfmt)     
            utils.PlotArrow(minPath_dl2l_num.finalPos, minPath_dl2l_num.finalHead, arrowLength, arrfmt)    
            utils.PlotArrow(minPath_dl2l_num.iniPos, minPath_dl2l_num.iniHead, arrowLength, arrfmt)            
            arrfmt = SimpleNamespace(color='g', linewidth=1, linestyle='--', marker='x', arrLen=arrowLength)                   
            utils.PlotInterval(minPath_dl2l_num.iniPos, int1, arrfmt)
            
        plt.axis('equal')
        
        ########### Plot from Dubins to paralellogram  #############
        # plotformat = namedtuple("plotformat","color linewidth linestyle marker")
        # plt.figure()
        # ut2.PlotLineSeg(line1_prlgrm[0], line1_prlgrm[1], plotformat('g',2,'-',''))
        # ut2.PlotLineSeg(line2_prlgrm[0], line2_prlgrm[1], plotformat('g',2,'-',''))

        # if minPathType2 != 'None':            
        #     du2.PlotDubPathSegments(minConfStart, minPathType2, minPathSegLengths,rho, plotformat('b',2,'-',''))
        #     ut2.PlotArrow(minConfStart[0:2], sector1[0], arrowLength, plotformat('c',2, '-','x'))
        #     ut2.PlotArrow(minConfStart[0:2], sector1[1], arrowLength, plotformat('c',2,'-','x'))
        #     ut2.PlotArrow(minConfGoal[0:2], sector2[0], arrowLength, plotformat('c',2,'-','x'))
        #     ut2.PlotArrow(minConfGoal[0:2], sector2[1], arrowLength, plotformat('c',2,'-','x'))
        #     ut2.PlotArrow(minConfStart[0:2], minConfStart[2], arrowLength, plotformat('m',2,'dotted','x'))
        #     ut2.PlotArrow(minConfGoal[0:2], minConfGoal[2], arrowLength, plotformat('m',2,'dotted','x'))

        plt.axis('equal')

        plt.show()
