
d = 1; // useless if transfinite
dy = 1;
dx = 5;
dz = 0.25;

nx = 2*dx;
ny = 2*dy;
nz = 1;

Point(301) = {0, 0, 0};
Point(302) = {dx, 0, 0};
Point(303) = {dx, dy, 0};
Point(304) = {0, dy, 0};
Point(305) = {0, 0, dz};
Point(306) = {dx, 0, dz};
Point(307) = {dx, dy, dz};
Point(308) = {0, dy, dz};
//face dessous

Line(401) = {301, 302};
Line(402) = {302, 303};
Line(403) = {303, 304};
Line(404) = {304, 301};
//face dessus
Line(405) = {305, 306};
Line(406) = {306, 307};
Line(407) = {307, 308};
Line(408) = {308, 305};
//lignes selon z
Line(409) = {301, 305};
Line(410) = {302, 306};
Line(411) = {303, 307};
Line(412) = {304, 308};

//face dessous
Curve Loop(51) = {401, 402, 403, 404};
//face dessus
Curve Loop(52) = {405, 406, 407, 408};
//faces des cotés
Curve Loop(53) = {401, 410, -405, -409};
Curve Loop(54) = {402, 411, -406, -410};
Curve Loop(55) = {403, 412, -407, -411};
Curve Loop(56) = {404, 409, -408, -412};

Plane Surface(61) = {51};
Plane Surface(62) = {52};
Plane Surface(63) = {53};
Plane Surface(64) = {54};
Plane Surface(65) = {55};
Plane Surface(66) = {56};

// Volume definition
Surface Loop(1) = {61, 62, 63, 64, 65, 66};
Volume(1) = {1};


// Physical groups definitions.
Physical Surface("WALL1") = {61}; //dessus  z=dz
Physical Surface("WALL2") = {62}; //dessous z=0
Physical Surface("WALL3") = {63}; //surface y=0
Physical Surface("WALL4") = {64}; //surface x=dx
Physical Surface("WALL5") = {65}; //surface y=dy
Physical Surface("WALL6") = {66}; //surface x=0

Physical Volume("MATERIAL1") = {1};

Mesh.Algorithm = 6;


Mesh.SaveAll = 1;
Mesh.Algorithm3D = 1;
Mesh.CharacteristicLengthMax = 0.1;