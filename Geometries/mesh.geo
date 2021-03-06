
d = 1; // useless if transfinite

nx = 1;
ny = 1;

Point(31) = {0, 0, 0, d};
Point(32) = {1, 0, 0, d};
Point(33) = {1, 1, 0, d};
Point(34) = {0, 1, 0, d};
Line(41) = {31, 32};
Line(42) = {32, 33};
Line(43) = {33, 34};
Line(44) = {34, 31};

Curve Loop(51) = {43, 44, 41, 42};
Plane Surface(61) = {51};

// transfinite mesh
//Transfinite Curve {44, 42} = ny+1 Using Progression 1;
//Transfinite Curve {43, 41} = nx+1 Using Progression 1;

//Transfinite Surface {61};

// Physical curves taking into account the boundaries.
Physical Curve("WALL1") = {44};
Physical Curve("WALL2") = {43};
Physical Curve("WALL3") = {42};
Physical Curve("WALL4") = {41};

//Physical Curve("Constant1", 1) = {44};
Physical Surface("MATERIAL1") = {61};

Mesh.Algorithm = 6;

Mesh.SaveAll = 1;