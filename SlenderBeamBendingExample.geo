//+
SetFactory("OpenCASCADE");
Rectangle(3) = {0, 0, 0, 1, 100, 0};
//+
Physical Surface(6) = {3};
//+
Transfinite Surface {3};
//+
Transfinite Curve {3, 1} = 2 Using Progression 1;
//+
Transfinite Curve {4, 2} = 11 Using Progression 1;
//+
Recombine Surface {3};