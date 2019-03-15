// Gmsh project created on Mon Mar 11 22:41:02 2019
SetFactory("OpenCASCADE");
//+
Torus(1) = {0, 0, -0, 0.5, 0.2, 2*Pi};
//+
Physical Volume("domain") = {1};
//+
Physical Surface("dirichlet") = {1};
