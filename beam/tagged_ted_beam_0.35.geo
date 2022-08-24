// Gmsh project created on Wed May 18 17:44:33 2022
//+
Point(1) = {0, 0, 0, 0.35};
//+
Point(2) = {0, 1, 0, 0.35};
//+
Point(3) = {0, 1, 1, 0.35};
//+
Point(4) = {0, 0, 1, 0.35};
//+
Point(5) = {4, 0, 0, 0.35};
//+
Point(6) = {4, 1, 0, 0.35};
//+
Point(7) = {4, 1, 1, 0.35};
//+
Point(8) = {4, 0, 1, 0.35};
//+
Point(9) = {8, 0, 0, 0.35};
//+
Point(10) = {8, 1, 0, 0.35};
//+
Point(11) = {8, 1, 1, 0.35};
//+
Point(12) = {8, 0, 1, 0.35};
//+
Line(1) = {2, 3};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 1};
//+
Line(4) = {1, 2};
//+
Line(5) = {6, 7};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 5};
//+
Line(8) = {5, 6};
//+
Line(9) = {2, 6};
//+
Line(10) = {5, 1};
//+
Line(11) = {4, 8};
//+
Line(12) = {7, 3};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {6, 7, 8, 5};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {12, -1, 9, 5};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {10, 4, 9, -8};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {10, -3, 11, 7};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {12, 2, 11, -6};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {4, 5, 1, 3, 6, 2};
//+
Volume(1) = {1};
//+
Line(13) = {6, 10};
//+
Line(14) = {5, 9};
//+
Line(15) = {7, 11};
//+
Line(16) = {8, 12};
//+
Line(17) = {11, 10};
//+
Line(18) = {10, 9};
//+
Line(19) = {11, 12};
//+
Line(20) = {12, 9};
//+
Curve Loop(7) = {8, 13, 18, -14};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {6, 16, -19, -15};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {7, 14, -20, -16};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {5, 15, 17, -13};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {17, 18, -20, -19};
//+
Plane Surface(11) = {11};
//+
Surface Loop(2) = {10, 8, 9, 7, 11, 2};
//+
Volume(2) = {2};
//+
Physical Surface("Fixed", 21) = {1};
//+
Physical Surface("Force", 22) = {11};
//+
Physical Volume("MaterialOne", 23) = {1};
//+
Physical Volume("MaterialTwo", 24) = {2};
//+
Physical Surface("Other", 25) = {5, 4, 3, 6, 9, 10, 7, 8};
