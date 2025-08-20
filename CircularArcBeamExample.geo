// Gmsh project created on Mon Aug 18 23:13:36 2025
SetFactory("OpenCASCADE");
//+

//+
Circle(2) = {0, 0, 0, 100, Pi/4, Pi/2};
//+
Circle(3) = {0, 0, 0, 99, Pi/4, Pi/2};
//+
Line(4) = {4, 2};
//+
Line(5) = {1, 3};
//+
Curve Loop(1) = {2, -4, -3, -5};
//+
Plane Surface(1) = {1};
//+
Transfinite Surface {1} = {4, 3, 1, 2};
//+
Transfinite Surface {1};
//+
Transfinite Curve {3, 2} = 9 Using Progression 1;
//+
Transfinite Curve {4, 5} = 2 Using Progression 1;
//+
Recombine Surface {1};
//+
Extrude {0, 0, 1} {
  Surface{1}; Layers {1}; Recombine;
}
//+
Physical Volume(14) = {1};
