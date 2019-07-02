lc = 1e-1;

Point(1) = {0, 0, 0, lc};
Point(2) = {2, 0, 0, lc};
Point(3) = {2, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {1,2,3,4};

Point(5) = {1.3, 0.1, 0, lc};
Point(6) = {1.3, 0.3, 0, lc};
Point(7) = {1.3, 0.5, 0, lc};
Point(8) = {1.3, 0.7, 0, lc};
Point(9) = {1.3, 0.9, 0, lc};

Circle(5) = {5,7,9};
Circle(6) = {9,7,5};

Curve Loop(2) = {5,6};

Circle(7) = {6,7,8};
Circle(8) = {8,7,6};

Curve Loop(3) = {7,8};

Plane Surface(1) = {1,2,3};
Plane Surface(2) = {2,3};
Plane Surface(3) = {3};

Physical Surface("Enviroment") = {1};
Physical Surface("Cell") = {2};
Physical Surface("Nucleoid") = {3};