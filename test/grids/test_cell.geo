// Refinement value
lc = 100;

// Cell and nucleioid points
Point(1) = {-2, 0, 0, lc};
Point(2) = {-1, 0, 0, lc};
Point(3) = {0, 0, 0, lc};
Point(4) = {1, 0, 0, lc};
Point(5) = {2, 0, 0, lc};

// Cell lines
Circle(1) = {1,3,5};
Circle(2) = {5,3,1};

// Cell curve
Curve Loop(1) = {1,2};

// Nucleoid lines
Circle(3) = {2,3,4};
Circle(4) = {4,3,2};

// Nuceloid curve
Curve Loop(2) = {3,4};


Plane Surface(1) = {1,2};
Physical Surface("Cell") = {1};

Plane Surface(2) = {2};
Physical Surface("Nucleoid") = {2};
