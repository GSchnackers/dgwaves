// The following geometry is of the form.
// -----------------------------------------------
// |                   e3                         |
// -----------------------------------------------
// |                   e2                         | e
// -----------------------------------------------
// |                   e1                         | h    dy
// -----------------------------------------------
// |                   e2                         | e
// -----------------------------------------------
// |                   e3                         | 
// -----------------------------------------------
//
//                     dx


dx = 5; // Length of the cable
dy = 2; // width of the cable
h = 0.5; // width of the central cable
e = 0.25; // width of the dielectric layers.

// Points Delimitating the domain.

Point(1) = {0, -dy/2, 0};
Point(2) = {dx, -dy/2, 0};
Point(3) = {dx, dy/2, 0};
Point(4) = {0, dy/2, 0};

// Points delimitating the media

Point(5) = {0, h/2, 0};
Point(6) = {dx, h/2, 0};
Point(7) = {dx, h/2 + e, 0};
Point(8) = {0, h/2 + e, 0};

Point(9) = {0, -h/2, 0};
Point(10) = {dx, -h/2, 0};
Point(11) = {dx, -h/2 - e, 0};
Point(12) = {0, -h/2 - e, 0};

// Lines delimitating the central cable.

Line(1) = {5, 9}; // Entry of the domain.
Line(2) = {9, 10}; 
Line(3) = {10, 6};
Line(4) = {6, 5};

// Lines defining the upper e zone

Line(5) = {5, 8}; // Opening.
Line(6) = {8, 7};
Line(7) = {7, 6}; // Opening

// Lines defining the highest zone.

Line(8) = {8, 4}; // Opening.
Line(9) = {4, 3}; // Opening 
Line(10) = {3, 7}; // Opening

// Lines defining the lower e zone

Line(11) = {9, 12}; // Opening.
Line(12) = {12, 11};
Line(13) = {11, 10}; // Opening

// Lines defining the lowest zone.

Line(14) = {12, 1}; // Opening.
Line(15) = {1, 2}; // Opening 
Line(16) = {2, 11}; // Opening

// Surface of the center cable 
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Surface of the upper e zone.
Curve Loop(2) = {5, 6, 7, 4};
Plane Surface(2) = {2};

// Surface of the highest zone.
Curve Loop(3) = {8, 9, 10, -6};
Plane Surface(3) = {3};

// Surface of the lower e zone.
Curve Loop(4) = {2, -13, -12, -11};
Plane Surface(4) = {4};

// Surface of the lowest zone.
Curve Loop(5) = {12, -16, -15, -14};
Plane Surface(5) = {5};

// Physical groups.
Physical Curve("WALL1") = {1};
Physical Curve("WALL2") = {3, 5, 7, 8, 9, 10, 11, 13, 14, 15, 16};

Physical Surface("MATERIAL1") = {1};
Physical Surface("MATERIAL2") = {2, 4};
Physical Surface("MATERIAL3") = {3, 5};

Mesh.Algorithm = 6;

