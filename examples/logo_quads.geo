Point (1) = {1.000000000000,1.000000000000,0.000000000000};
Point (2) = {1.000000000000,0.379757979263,-0.000000000000};
Point (3) = {0.801534880751,0.454963545532,-0.000000000000};
Point (4) = {0.656107331955,1.000000000000,-0.000000000000};
Point (5) = {0.600672553037,0.775245709538,-0.000000000000};
Point (6) = {0.000000000000,1.000000000000,0.000000000000};
Point (7) = {0.392825178821,0.672136259831,-0.000000000000};
Point (8) = {1.000000000000,0.000000000000,0.000000000000};
Point (9) = {0.547800422194,-0.000000000000,-0.000000000000};
Point (10) = {0.488710023938,0.224380304618,-0.000000000000};
Point (11) = {0.000000000000,0.000000000000,0.000000000000};
Point (12) = {-0.000000000000,0.324566579562,-0.000000000000};
Point (13) = {0.172066723668,0.367888021869,-0.000000000000};
Line (1) = {2,1};
Line (2) = {1,4};
Line (3) = {2,3};
Line (4) = {3,5};
Line (5) = {5,4};
Line (6) = {4,6};
Line (7) = {7,6};
Line (8) = {5,7};
Line (9) = {8,2};
Line (10) = {9,8};
Line (11) = {9,10};
Line (12) = {10,3};
Line (13) = {12,11};
Line (14) = {11,9};
Line (15) = {12,13};
Line (16) = {13,10};
Line (17) = {6,12};
Line (18) = {7,13};
Line Loop (1) = {-1,3,4,5,-2};
Plane Surface (1) = {1}; Physical Surface (1) = {1};
Line Loop (2) = {7,-6,-5,8};
Plane Surface (2) = {2}; Physical Surface (2) = {2};
Line Loop (3) = {-10,11,12,-3,-9};
Plane Surface (3) = {3}; Physical Surface (3) = {3};
Line Loop (4) = {-11,-14,-13,15,16};
Plane Surface (4) = {4}; Physical Surface (4) = {4};
Line Loop (5) = {-7,18,-15,-17};
Plane Surface (5) = {5}; Physical Surface (5) = {5};
Line Loop (6) = {-16,-18,-8,-4,-12};
Plane Surface (6) = {6}; Physical Surface (6) = {6};
Recombine Surface{5};

ReverseMesh Surface{1,2,3,4,5,6};

Physical Curve("left") = {17,13};
Physical Curve("bottom") = {14,10};
Physical Curve("right") = {9,1};
Physical Curve("top") = {2,6};

Physical Surface("green") = {6};
