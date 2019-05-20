Point(1) = {0,0,0};
Point(2) = {1,0,0};
Point(3) = {1,1,0};
Point(4) = {0,1,0};

Point(5) = {0.49,0,0};
Point(6) = {0.51,0,0};
Point(7) = {0.49,0.1,0};
Point(8) = {0.51,0.1,0};

Point(9) = {0.49,1,0};
Point(10) = {0.51,1,0};
Point(11) = {0.49,0.9,0};
Point(12) = {0.51,0.9,0};

//+
Line(1) = {4, 9};
//+
Line(2) = {9, 11};
//+
Line(3) = {11, 12};
//+
Line(4) = {12, 10};
//+
Line(5) = {10, 3};
//+
Line(6) = {3, 2};
//+
Line(7) = {2, 6};
//+
Line(8) = {6, 8};
//+
Line(9) = {8, 7};
//+
Line(10) = {7, 5};
//+
Line(11) = {5, 1};
//+
Line(12) = {1, 4};
//+
Curve Loop(1) = {12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
//+
Plane Surface(1) = {1};
//+
Physical Curve("WALL1") = {12};
//+
Physical Curve("WALL2") = {6};
//+
Physical Curve("WALL3") = {1};
//+
Physical Curve("WALL4") = {2};
//+
Physical Curve("WALL5") = {3};
//+
Physical Curve("WALL6") = {4};
//+
Physical Curve("WALL7") = {5};
//+
Physical Curve("WALL8") = {11};
//+
Physical Curve("WALL9") = {10};
//+
Physical Curve("WALL10") = {9};
//+
Physical Curve("WALL11") = {8};
//+
Physical Curve("WALL12") = {7};
//+
Physical Surface("MATERIAL1") = {1};
