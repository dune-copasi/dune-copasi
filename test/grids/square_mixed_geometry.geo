// Refinement value
lc = 0.5;

// Square points
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 0.5, 0, lc};
Point(4) = {1, 1, 0, lc};
Point(5) = {0, 1, 0, lc};
Point(6) = {0, 0.5, 0, lc};

// Square line
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,6};
Line(4) = {6,1};


Line(5) = {6,3};
Line(6) = {3,4};
Line(7) = {4,5};
Line(8) = {5,6};

Transfinite Curve {1, 3, 5} = 3 Using Progression 1;
Transfinite Curve {4, 2} = 3 Using Progression 1;

// Square curve
Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Physical Surface("Cubes") = {1};



Curve Loop(2) = {5,6,7,8};
Plane Surface(2) = {2};
Physical Surface("Triangles") = {2};

Transfinite Surface {1};
Recombine Surface{1};
// //+
// Transfinite Curve {4, 2} = 10 Using Progression 1;
// //+
// Transfinite Curve {3, 1} = 10 Using Progqression 1;
