//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 0, 1, 1.0};
//+
Point(3) = {0, 1, 0, 1.0};
//+
Point(4) = {0, 1, 1, 1.0};
//+
Point(5) = {1, 0, 0, 1.0};
//+
Point(6) = {1, 0, 1, 1.0};
//+
Point(7) = {1, 1, 0, 1.0};
//+
Point(8) = {1, 1, 1, 1.0};
//+
Line(1) = {4, 3};
//+
Line(2) = {3, 1};
//+
Line(3) = {1, 2};
//+
Line(4) = {2, 4};
//+
Line(5) = {3, 7};
//+
Line(6) = {1, 5};
//+
Line(7) = {5, 6};
//+
Line(8) = {5, 7};
//+
Line(9) = {6, 8};
//+
Line(10) = {8, 7};
//+
Line(11) = {8, 4};
//+
Line(12) = {2, 6};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {6, 8, -5, 2};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {10, -5, -1, -11};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {3, 12, -7, -6};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {8, -10, -9, -7};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {12, 9, 11, -4};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {1, 3, 5, 2, 4, 6};
//+
Volume(1) = {1};
//+
Physical Surface("Sinusoidal1") = {4};
//+
Physical Volume("Domain") = {1};//+
Transfinite Surface {4};
//+
Transfinite Surface {6};
//+
Transfinite Surface {5};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
//Transfinite Curve {2, 1, 4, 3, 5, 6, 12, 11, 9, 7, 8, 10} = 2 Using Progression 1;
