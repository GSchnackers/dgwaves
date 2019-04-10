
d = 1; // useless if transfinite

nx = 2;
ny = 2;

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
Transfinite Curve {44, 42} = ny+1 Using Progression 1;
Transfinite Curve {43, 41} = nx+1 Using Progression 1;

Transfinite Surface {61};

// Physical curves taking into account the boundaries.
Physical Curve("Sinusoidal1", 1) = {44};
Physical Curve("Wall1", 2) = {43};
Physical Curve("Wall2", 3) = {41};
Physical Curve("Output1", 4) = {42};

Physical Surface("Domain") = {61};

Mesh.SaveAll = 1;