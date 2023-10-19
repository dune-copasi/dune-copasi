// Refinement value
refinement = DefineNumber[ 2, Name "Parameters/refinement" ];

// Cytoplasm radious
cr = 2;
// Nucleus radious
nr = 0.7;
// Membrane width
mw = 0.01;

// Cytoplasm, Nucleus, and Membrane points
Point(1) = {-cr, 0, 0, refinement};
Point(2) = {-nr-mw, 0, 0, refinement};
Point(3) = {-nr, 0, 0, refinement};
Point(4) = {0, 0, 0, refinement};
Point(5) = {nr, 0, 0, refinement};
Point(6) = {nr+mw, 0, 0, refinement};
Point(7) = {cr, 0, 0, refinement};

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

membrane_size = DefineNumber[ 20, Name "Parameters/membrane_size" ];

// Transfinite Curve {5, 6, 4, 3} = membrane_size Using Progression 1;

// Recombine Surface {2};
