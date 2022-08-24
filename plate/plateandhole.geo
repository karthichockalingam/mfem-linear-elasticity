//+
Point(1) = {1, 0, 0, 0.05};
//+
Point(2) = {4, 0, 0, 0.05};
//+
Point(3) = {4, 4, 0, 0.05};
//+
Point(4) = {0, 4, 0, 0.05};
//+
Point(5) = {0, 1, 0, 0.05};
//+
Point(6) = {0, 0, 0, 0.05};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Circle(5) = {5, 6, 1};
//+
Curve Loop(1) = {4, 5, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Curve("bottom", 6) = {1};
//+
Physical Curve("right", 7) = {2};
//+
Physical Curve("top", 8) = {3};
//+
Physical Curve("left", 9) = {4};
//+
Physical Curve("arch", 10) = {5};
//+
Physical Surface("material", 11) = {1};
