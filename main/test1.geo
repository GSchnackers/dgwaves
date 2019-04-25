
d = 1; // useless if transfinite

nx = 2;
ny = 2;

Point(31) = {0, 0, 0, d};
Point(32) = {0.5, 0, 0, d};
Point(33) = {1, 0, 0, d};
Point(34) = {1, 1, 0, d};
Point(35) = {0.5, 1, 0, d};
Point(36) = {0, 1, 0, d};
Line(41) = {31, 32};
Line(42) = {32, 35};
Line(43) = {35, 36};
Line(44) = {36, 31};
Line(45) = {32, 33};
Line(46) = {33, 34};
Line(47) = {34, 35};
Line(48) = {35, 32};

Curve Loop(51) = {41, 42, 43, 44};
Curve Loop(52) = {45, 46, 47, 48};
Plane Surface(61) = {51};
Plane Surface(62) = {52};

// transfinite mesh
Transfinite Curve {42, 44, 46, 48} = ny+1 Using Progression 1;
Transfinite Curve {41, 43, 45, 47} = nx+1 Using Progression 1;

Transfinite Surface {61};
Transfinite Surface {62};

// Physical curves taking into account the boundaries.
Physical Curve("Sinusoidal1", 1) = {44};
//Physical Curve("Constant1", 1) = {44};
Physical Surface("Domain1") = {61};
Physical Surface("Domain2") = {62};

Mesh.SaveAll = 1;