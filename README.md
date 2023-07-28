# Dubins line-segment To line-segment
Shortest Dubins Line-segment-Interval to Line-segment-Interval
Given the initial and final line segments and intervals, and this code finds the shortest Dubins path where the end positions (initial/final) lie on the given line-segments, and the initial/final headings lie in the given initial/final intervals.

The shortest path is computed by finding the candidate optimal paths for each mdoe of Dubins paths: LSL, LSR, RSL, RSR, LRL, RLR, LS, RS, SL, SR, LR, RL, L, R, and S.

Check the "Example_DubL2L.py" for usage.

Required packages:
numpy
matplotlib (pyplot for plotting)
dubins
dataclasses
timit (for checking computation time)

This code uses the shortest Dubins path from line-segment-interval to line-segment-interval from the paper below
https://www.roboticsproceedings.org/rss19/p059.html


Citation:
@INPROCEEDINGS{Manyam-RSS-23, 
    AUTHOR    = {Satyanarayana Gupta Manyam AND Abhishek Nayak AND Sivakumar  Rathinam}, 
    TITLE     = {{G*: A New Approach to Bounding Curvature Constrained Shortest Paths through Dubins Gates}}, 
    BOOKTITLE = {Proceedings of Robotics: Science and Systems}, 
    YEAR      = {2023}, 
    ADDRESS   = {Daegu, Republic of Korea}, 
    MONTH     = {July}, 
    DOI       = {10.15607/RSS.2023.XIX.059} 
}
