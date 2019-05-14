d = 1;

nx = 20;
ny = 20;
nz = 1;
L = 5;
W = 1;
H = 1;

Point(101) = {0, -W/2, -H/2, d};
Point(102) = {0, W/2, -H/2, d};
Point(103) = {0, W/2, H/2, d};
Point(104) = {0, -W/2, H/2, d};
Point(105) = {L/2, -W/2, -H/2, d};
Point(106) = {L/2, W/2, -H/2, d};
Point(107) = {L/2, W/2, H/2, d};
Point(108) = {L/2, -W/2, H/2, d};
Point(109) = {L, -W/2, -H/2, d};
Point(110) = {L, W/2, -H/2, d};
Point(111) = {L, W/2, H/2, d};
Point(112) = {L, -W/2, H/2, d};

Line(201) = {101, 102};
Line(202) = {102, 103};
Line(203) = {103, 104};
Line(204) = {104, 101};
Line(205) = {105, 106};
Line(206) = {106, 107};
Line(207) = {107, 108};
Line(208) = {108, 105};
Line(209) = {109, 110};
Line(210) = {110, 111};
Line(211) = {111, 112};
Line(212) = {112, 109};
Line(213) = {101, 105};
Line(214) = {105, 109};
Line(215) = {102, 106};
Line(216) = {106, 110};
Line(217) = {103, 107};
Line(218) = {107, 111};
Line(219) = {104, 108};
Line(220) = {108, 112};

Curve Loop(301) = {201, 202, 203, 204};
Plane Surface(401) = {301}; //Face

Curve Loop(302) = {209, 210, 211, 212};
Plane Surface(402) = {302}; //Arrière

Curve Loop(303) = {201, 215, -205, -213};
Plane Surface(403) = {303}; //Bas 1

Curve Loop(304) = {205, 216, -209, -214};
Plane Surface(404) = {304}; //Bas 2

Curve Loop(305) = {202, 217, -206, -215};
Plane Surface(405) = {305}; //Droite 1

Curve Loop(306) = {206, 218, -210, -216};
Plane Surface(406) = {306}; //Droite 2

Curve Loop(307) = {203, 219, -207, -217};
Plane Surface(407) = {307}; //Haut 1

Curve Loop(308) = {207, 220, -211, -218};
Plane Surface(408) = {308}; //Haut 2

Curve Loop(309) = {204, 213, -208, -219};
Plane Surface(409) = {309}; //Gauche 1

Curve Loop(310) = {208, 214, -212, -220};
Plane Surface(410) = {310}; //Gauche 2

Surface Loop(1) = {407, 401, 403, 405, 406, 408, 410, 409, 404, 402};
Volume(1) = {1};

//Transfinite Curve {213, 214, 215, 216, 217, 218, 219, 220} = nx + 1 Using Progression 1;
//Transfinite Curve {201, 203, 205, 207, 209, 211} = ny + 1 Using Progression 1;
//Transfinite Curve {202, 204, 206, 208, 210, 212} = nz + 1 Using Progression 1;

//Transfinite Surface {501} = {101, 102, 103, 104};
//Transfinite Surface {502} = {105, 106, 107, 108};
//Transfinite Surface {503} = {101, 102, 106, 105};
//Transfinite Surface {504} = {105, 106, 109, 108};
//Transfinite Surface {505} = {102, 103, 107, 106};
//Transfinite Surface {506} = {106, 107, 110, 109};
//Transfinite Surface {507} = {103, 104, 108, 107};
//Transfinite Surface {508} = {107, 108, 111, 110};
//Transfinite Surface {509} = {104, 101, 108, 105};
//Transfinite Surface {510} = {105, 108, 112, 109};

// Physical groups definitions.
Physical Surface("WALL1") = {401}; //Face
Physical Surface("WALL2") = {402}; //Arrière
Physical Surface("WALL3") = {403}; //Bas 1
Physical Surface("WALL4") = {404}; //Bas 2
Physical Surface("WALL5") = {405}; //Droite 1
Physical Surface("WALL6") = {406}; //Droite 2
Physical Surface("WALL7") = {407}; //Haut 1
Physical Surface("WALL8") = {408}; //Haut 2
Physical Surface("WALL9") = {409}; //Gauche 1
Physical Surface("WALL10") = {410}; //Gauche 2

Physical Volume("MATERIAL1") = {1};

Mesh.Algorithm = 6;