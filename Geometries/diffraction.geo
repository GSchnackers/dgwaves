Point(1) = {0,0,0};
Point(2) = {5,0,0};
Point(3) = {5,1,0};
Point(4) = {0,1,0};

Point(5) = {3,0,0};
Point(6) = {3,1,0};
Point(7) = {3,0.4,0};
Point(8) = {3,0.6,0};

Point(9) = {3.1,0.4,0};
Point(10) = {3.1,0.6,0};

Point(11) = {3.1,0,0};
Point(12) = {3.1,1,0};

Point(13) = {2,0,0};
Point(14) = {2,1,0};


//+
Line(1) = {4, 14};
//+
Line(2) = {14, 13};
//+
Line(3) = {13, 1};
//+
Line(4) = {1, 4};
//+
Line(5) = {14, 6};
//+
Line(6) = {6, 8};
//+
Line(7) = {8, 10};
//+
Line(8) = {10, 12};
//+
Line(9) = {12, 3};
//+
Line(10) = {3, 2};
//+
Line(11) = {2, 11};
//+
Line(12) = {11, 9};
//+
Line(13) = {9, 7};
//+
Line(14) = {7, 5};
//+
Line(15) = {5, 13};
//+
Line(16) = {8, 7};
//+
Line(17) = {10, 9};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 6, 16, 14, 15, -2};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {16, -13, -17, -7};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {8, 9, 10, 11, 12, -17};
//+
Plane Surface(4) = {4};

Mesh.Algorithm = 6;

//+
Physical Curve("WALL1") = {4};
//+
Physical Curve("WALL2") = {1};
//+
Physical Curve("WALL3") = {5};
//+
Physical Curve("WALL4") = {6};
//+
Physical Curve("WALL5") = {7};
//+
Physical Curve("WALL6") = {8};
//+
Physical Curve("WALL7") = {9};
//+
Physical Curve("WALL8") = {3};
//+
Physical Curve("WALL9") = {15};
//+
Physical Curve("WALL10") = {14};
//+
Physical Curve("WALL11") = {13};
//+
Physical Curve("WALL12") = {12};
//+
Physical Curve("WALL13") = {11};
//+
Physical Curve("WALL14") = {10};
//+
Physical Surface("MATERIAL1") = {1};
//+
Physical Surface("MATERIAL2") = {2};
//+
Physical Surface("MATERIAL3") = {3};
//+
Physical Surface("MATERIAL4") = {4};
