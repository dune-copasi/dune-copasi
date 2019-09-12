// Refinement value
r = 2e-1;
mr = 4e-2;


// Cytoplasm radious
cr = 2;
// Nucleus radious
nr = 0.7;
// Membrane width
mw = 0.05;

// Cytoplasm, Nucleus, and Membrane points
Point(1) = {-cr, 0, 0, r};
Point(2) = {-nr-mw, 0, 0, mr};
Point(3) = {-nr, 0, 0, mr};
Point(4) = {0, 0, 0, r};
Point(5) = {nr, 0, 0, mr};
Point(6) = {nr+mw, 0, 0, mr};
Point(7) = {cr, 0, 0, r};

// outer cytoplasm lines
Circle(1) = {1,4,7};
Circle(2) = {7,4,1};

// outer cytoplasm curve
Curve Loop(1) = {1,2};

// inner cytoplasm lines
Circle(3) = {2,4,6};
Circle(4) = {6,4,2};

// inner cytoplasm curve
Curve Loop(2) = {3,4};

// outer nucleus lines
Circle(5) = {3,4,5};
Circle(6) = {5,4,3};

// outer cytoplasm curve
Curve Loop(3) = {5,6};


Plane Surface(1) = {1,2};
Physical Surface("Cell") = {1};

Plane Surface(2) = {2,3};
Physical Surface("Membrane") = {2};

Plane Surface(3) = {3};
Physical Surface("Nucleus") = {3};

Point{4} In Surface{3};