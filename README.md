# DubinsLineSegToLineSeg
Shortest Dubins Line-segment-Interval to Line-segment-Interval
Given the initial and final line segments and intervals, and this code finds the shortest Dubins path where the end positions (initial/final) lie on the given line-segments, and the initial/final headings lie in the given initial/final intervals.

The shortest path is computed by finding the candidate optimal paths for each mdoe of Dubins paths: LSL, LSR, RSL, RSR, LRL, RLR, LS, RS, SL, SR, LR, RL, L, R, and S.

Check the example script for usage.
