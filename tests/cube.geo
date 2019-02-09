Lx=10;  // longueur 
Ly=10;   // largeur
Lz=10;   // hauteur
n=12;    // nbre d'elements 

Point(1) = {0, 0, 0, 1.0};
Point(2) = {0, Ly, 0, 1.0};
Point(3) = {Lx, Ly, 0, 1.0};
Point(4) = {Lx, 0, 0, 1.0};
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {2, 3};
Line(4) = {1, 2};
Line Loop(5) = {4, 3, -2, -1};
Plane Surface(6) = {5};
Transfinite Line {4, 2, 3, 1} = n+1 Using Progression 1;
Transfinite Surface {6} = {1, 4, 3, 2};
Recombine Surface {6};

Extrude {0, 0, Lz} {
  Surface{6}; Layers{n}; Recombine;
}
Physical Surface("Surface droite") = {23};
Physical Volume("Volume") = {1};
