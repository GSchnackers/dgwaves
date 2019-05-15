
d = 0.2; // useless if transfinite
h = 3;
l = 2;

m1 = 0.3;
m2 = 0.4;

nx = 4;
ny = 4*h/d;

Point(301) = {0, 0, 0, d}; 
Point(302) = {m1*l, 0, 0, d};
Point(303) = {m2*l, 0, 0, d};
Point(304) = {(1-m2)*l, 0, 0, d};
Point(305) = {(1-m1)*l, 0, 0, d};
Point(306) = {l, 0, 0, d};
Point(307) = {l, h, 0, d};
Point(308) = {(1-m1)*l, h, 0, d};
Point(309) = {(1-m2)*l, h, 0, d};
Point(310) = {m2*l, h, 0, d};
Point(311) = {m1*l, h, 0, d};
Point(312) = {0, h, 0, d};
Line(401) = {301, 302}; //horz bas 1
Line(402) = {302, 311}; //vert 2
Line(403) = {311, 312}; //horz haut 1
Line(404) = {312, 301}; //vert 1
Line(405) = {302, 303}; //horz bas 2
Line(406) = {303, 310}; //vert 3
Line(407) = {310, 311}; //horz haut 2
// Line(408) = {311, 302}; // - Line(402)
Line(409) = {303, 304}; //horz bas 3
Line(410) = {304, 309}; //vert 4
Line(411) = {309, 310}; //horz haut 3
// Line(412) = {310, 303}; // - Line(406)
Line(413) = {304, 305}; //horz bas 4
Line(414) = {305, 308}; //vert 5
Line(415) = {308, 309}; //horz haut 4
// Line(416) = {309, 304}; // - Line(410)
Line(417) = {305, 306}; //horz bas 5
Line(418) = {306, 307}; //vert 6
Line(419) = {307, 308}; //horz haut 5
// Line(420) = {308, 305}; // - Line(414)

Curve Loop(51) = {401, 402, 403, 404};
Curve Loop(52) = {405, 406, 407, -402};
Curve Loop(53) = {409, 410, 411, -406};
Curve Loop(54) = {413, 414, 415, -410};
Curve Loop(55) = {417, 418, 419, -414};
Plane Surface(61) = {51};
Plane Surface(62) = {52};
Plane Surface(63) = {53};
Plane Surface(64) = {54};
Plane Surface(65) = {55};

// transfinite mesh
//Transfinite Curve {402, 404, 406, 410, 414, 418} = ny+1 Using Progression 1;
//Transfinite Curve {401, 403, 405, 407, 409, 411, 413, 415, 417, 419} = nx+1 Using Progression 1;

//Transfinite Surface {61};
//Transfinite Surface {62};
//Transfinite Surface {63};
//Transfinite Surface {64};
//Transfinite Surface {65};

// Physical curves taking into account the boundaries.
Physical Curve("WALL1") = {409}; //horz bas 3
Physical Curve("WALL2") = {401}; //horz bas 1
Physical Curve("WALL3") = {405}; //horz bas 2
Physical Curve("WALL4") = {413}; //horz bas 4
Physical Curve("WALL5") = {417}; //horz bas 5
Physical Curve("WALL6") = {418}; //vert 6 (à droite)
Physical Curve("WALL7") = {419}; //horz haut 5
Physical Curve("WALL8") = {415}; //horz haut 4
Physical Curve("WALL9") = {411}; //horz haut 3
Physical Curve("WALL10") = {407}; //horz haut 2
Physical Curve("WALL11") = {403}; //horz haut 1
Physical Curve("WALL12") = {404}; //vert 1 (à gauche)
//Physical Curve("Constant1") = {409};
Physical Surface("MATERIAL3") = {61, 65};
Physical Surface("MATERIAL2") = {62, 64};
Physical Surface("MATERIAL1") = {63};

Mesh.Algorithm = 6;