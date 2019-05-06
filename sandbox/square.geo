d = 1;
L = 1;
n = 1;

// List of points
Point(1) = {0, 0, 0, d};
Point(2) = {L, 0, 0, d};
Point(3) = {L, L, 0, d};
Point(4) = {0, L, 0, d};

// List of lines
Line(11) = {1, 2};
Line(12) = {2, 3};
Line(13) = {3, 4};
Line(14) = {4, 1};

// Surface
Curve Loop(21) = {11, 12, 13, 14};
Plane Surface(31) = {21};

// Transfinite
Transfinite Curve {11, 12, 13, 14} = n+1 Using Progression 1;

Transfinite Surface {31};