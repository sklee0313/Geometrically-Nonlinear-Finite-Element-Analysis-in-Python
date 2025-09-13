// Gmsh project created on Mon Sep 01 21:10:03 2025
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 1, 100, 0};
//+
Transfinite Curve {3, 1} = 2 Using Progression 1;
//+
Transfinite Curve {4, 2} = 11 Using Progression 1;
//+
Transfinite Surface {1};
//+
Recombine Surface {1};
//+
Extrude {0, 0, 0.5} {
  Surface{1}; Layers {1}; Recombine;
}
//+
Physical Volume(13) = {1};
//+
Physical Surface(14) = {4};
