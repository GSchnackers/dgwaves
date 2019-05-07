
d = 1; // useless if transfinite
h = 3; // height
l = 2; // length
k = 0.2; // shift in ]0;0.3[
m = h/2; // height of points for spline
n = 2*k*m/h; // shift of points for spline
c = 0.8; // coefficients for spline in ]0;1[

nx = 2;
ny = 2;

Point(301) = {0, 0, 0, d};
Point(302) = {(0.3-k)*l, 0, 0, d};
Point(303) = {(0.4-k)*l, 0, 0, d};
Point(304) = {(0.6-k)*l, 0, 0, d};
Point(305) = {(0.7-k)*l, 0, 0, d};
Point(306) = {l, 0, 0, d};
Point(307) = {l, h, 0, d};
Point(308) = {(0.7+k)*l, h, 0, d};
Point(309) = {(0.6+k)*l, h, 0, d};
Point(310) = {(0.4+k)*l, h, 0, d};
Point(311) = {(0.3+k)*l, h, 0, d};
Point(312) = {0, h, 0, d};
Point(313) = {(0.3-k+n*c/2)*l, m*c, 0, d};
Point(314) = {(0.4-k+n*c/2)*l, m*c, 0, d};
Point(315) = {(0.6-k+n*c/2)*l, m*c, 0, d};
Point(316) = {(0.7-k+n*c/2)*l, m*c, 0, d};
Point(317) = {(0.3+k-n*c/2)*l, h-m*c, 0, d};
Point(318) = {(0.4+k-n*c/2)*l, h-m*c, 0, d};
Point(319) = {(0.6+k-n*c/2)*l, h-m*c, 0, d};
Point(320) = {(0.7+k-n*c/2)*l, h-m*c, 0, d};
Line(401) = {301, 302};
Spline(402) = {302, 313, 317, 311};
Line(403) = {311, 312};
Line(404) = {312, 301};
Line(405) = {302, 303};
Spline(406) = {303, 314, 318, 310};
Line(407) = {310, 311};
Spline(408) = {311, 317, 313, 302};
Line(409) = {303, 304};
Spline(410) = {304, 315, 319, 309};
Line(411) = {309, 310};
Spline(412) = {310, 318, 314, 303};
Line(413) = {304, 305};
Spline(414) = {305, 316, 320, 308};
Line(415) = {308, 309};
Spline(416) = {309, 319, 315, 304};
Line(417) = {305, 306};
Line(418) = {306, 307};
Line(419) = {307, 308};
Spline(420) = {308, 320, 316, 305};

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
