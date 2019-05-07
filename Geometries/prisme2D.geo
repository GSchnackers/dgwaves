d = 0.05;
L = 1;
H = 1;
C = 0.2;
b = 1;
h = 2;
nx = 20;
ny = 20;
nxy = 20;

Point(1) = {0, 0, 0, d};
Point(2) = {C+L, 0, 0, d};
Point(3) = {C+L, H, 0, d};
Point(4) = {0, H, 0, d};
Point(5) = {C, 0, 0, d};
Line(11) = {1, 5};
Line(12) = {5, 3};
Line(13) = {3, 4};
Line(14) = {4, 1};
Line(15) = {5, 2};
Line(16) = {2, 3};
Line(17) = {3, 5};

Curve Loop(21) = {11, 12, 13, 14};
Curve Loop(22) = {15, 16, 17};
Plane Surface(31) = {21};
Plane Surface(32) = {22};

// transfinite mesh
//Transfinite Curve {13, 15} = ny+1 Using Progression 1;
//Transfinite Curve {12, 14} = nx+1 Using Progression 1;
//Transfinite Curve {11, 16} = nxy+1 Using Progression 1;

//Transfinite Surface {31} = {4, 1, 5, 3};
//Transfinite Surface {32} = {2, 3, 5};

// Physical curves taking into account the boundaries.
Physical Curve("WALL1") = {14};
Physical Curve("WALL2") = {11};
Physical Curve("WALL3") = {15};
Physical Curve("WALL4") = {16};
Physical Curve("WALL5") = {13};
Physical Surface("MATERIAL1") = {31};
Physical Surface("MATERIAL2") = {32};

Mesh.Algorithm = 6;