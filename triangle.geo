// Gmsh project created on Wed Jan 30 18:25:19 2019
//+
Point(1) = {0, 0, 0, 0.1};
//+
Point(2) = {1, 0.2, 0, 0.02};
//+
Point(3) = {0.4, 0.6, 0, 0.05};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 1};
//+
Curve Loop(1) = {3, 1, 2};
//+
Plane Surface(1) = {1};
