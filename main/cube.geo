d = 1;

// Points definition
For i In {0:1}
    For j In {0:1}
        For k In {0:1}

            Point(4 * i + 2 * j + k + 1) = {i, j, k, d};

        EndFor
    EndFor
EndFor

// Edges definitions
count = 1;
For i In {0:1}
    For j In {0:1}

        Line(count) = {2 * i + j + 1, 4 + 2 * i + j + 1}; // Have the same y and z but differ from x.
        Line(count + 1) = {4 * i + j + 1, 4 * i + 2 + j + 1}; // Have the same x and z but differ from y
        Line(count + 2) = {4 * i + 2 * j + 1, 4 * i + 2 * j + 2}; // Have the same x and y but differ from z. 
        count = count + 3;

    EndFor
EndFor

// Definition of all surrounding surfaces.
Curve Loop(1) = {6, -5, -3, 2};
Plane Surface(1) = {1};

Curve Loop(2) = {1, 9, -4, -3};
Plane Surface(2) = {2};

Curve Loop(3) = {1, 8, -7, -2};
Plane Surface(3) = {3};

Curve Loop(4) = {5, 10, -11, -4};
Plane Surface(4) = {4};

Curve Loop(5) = {8, 12, -11, -9};
Plane Surface(5) = {5};

Curve Loop(6) = {10, -12, -7, 6};
Plane Surface(6) = {6};

// Volume definition
Surface Loop(1) = {1, 6, 4, 5, 3, 2};
Volume(1) = {1};

// Mesh method definition of the edges.
Transfinite Curve {6, 2, 3, 5, 4, 1, 9, 11, 8, 10, 12, 7} = 2 Using Progression 1;

// Mesh method definition on the surfaces.
For i In {1:6}
    Transfinite Surface {i};
EndFor

// Mesh definition on the volume.
Transfinite Volume{1} = {3, 4, 8, 7, 1, 2, 6, 5};

// Physical groups definitions.
Physical Surface("Sinusoidal1") = {1};
Physical Volume("Domain") = {1};
Physical Volume("Domain") += {1};
