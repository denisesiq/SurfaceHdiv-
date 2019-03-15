// Gmsh project created on Mon Mar 11 21:50:20 2019
SetFactory("OpenCASCADE");
//+
Sphere(1) = {0, 0, 0, 1, -Pi/2, Pi/2, 2*Pi};
//+
Physical Surface("domain") = {1};
//+
Physical Curve("boundary") = {2};
