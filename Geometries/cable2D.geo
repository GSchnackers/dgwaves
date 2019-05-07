
d = 1; // useless if transfinite
h = 3;
l = 2;

nx = 2;
ny = 2;

Point(301) = {0, 0, 0, d};
Point(302) = {0.3*l, 0, 0, d};
Point(303) = {0.4*l, 0, 0, d};
Point(304) = {0.6*l, 0, 0, d};
Point(305) = {0.7*l, 0, 0, d};
Point(306) = {l, 0, 0, d};
Point(307) = {l, h, 0, d};
Point(308) = {0.7*l, h, 0, d};
Point(309) = {0.6*l, h, 0, d};
Point(310) = {0.4*l, h, 0, d};
Point(311) = {0.3*l, h, 0, d};
Point(312) = {0, h, 0, d};
Line(401) = {301, 302};
Line(402) = {302, 311};
Line(403) = {311, 312};
Line(404) = {312, 301};
Line(405) = {302, 303};
Line(406) = {303, 310};
Line(407) = {310, 311};
Line(408) = {311, 302};
Line(409) = {303, 304};
Line(410) = {304, 309};
Line(411) = {309, 310};
Line(412) = {310, 303};
Line(413) = {304, 305};
Line(414) = {305, 308};
Line(415) = {308, 309};
Line(416) = {309, 304};
Line(417) = {305, 306};
Line(418) = {306, 307};
Line(419) = {307, 308};
Line(420) = {308, 305};

Curve Loop(51) = {401, 402, 403, 404};
Curve Loop(52) = {405, 406, 407, 408};
Curve Loop(53) = {409, 410, 411, 412};
Curve Loop(54) = {413, 414, 415, 416};
Curve Loop(55) = {417, 418, 419, 420};
Plane Surface(61) = {51};
Plane Surface(62) = {52};
Plane Surface(63) = {53};
Plane Surface(64) = {54};
Plane Surface(65) = {55};

// transfinite mesh
Transfinite Curve {402, 404, 406, 408, 410, 412, 414, 416, 418, 420} = ny+1 Using Progression 1;
Transfinite Curve {401, 403, 405, 407, 409, 411, 413, 415, 417, 419} = nx+1 Using Progression 1;

Transfinite Surface {61};
Transfinite Surface {62};
Transfinite Surface {63};
Transfinite Surface {64};
Transfinite Surface {65};

// Physical curves taking into account the boundaries.
Physical Curve("Sinusoidal1", 1) = {409};
//Physical Curve("Constant1", 1) = {409};
Physical Surface("Domain1") = {61, 65};
Physical Surface("Domain2") = {62, 64};
Physical Surface("Domain3") = {63};

Mesh.SaveAll = 1;