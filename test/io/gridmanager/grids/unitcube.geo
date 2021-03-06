size = 0.2;
Point(1) = {0.0, 0.0, 0.0, size};
Point(2) = {1.0, 0.0, 0.0, size};
Point(3) = {1.0, 1.0, 0.0, size};
Point(4) = {0.0, 1.0, 0.0, size};
Point(5) = {0.0, 0.0, 1.0, size};
Point(6) = {1.0, 0.0, 1.0, size};
Point(7) = {1.0, 1.0, 1.0, size};
Point(8) = {0.0, 1.0, 1.0, size};

Line(1) = {4, 8};
Line(2) = {8, 7};
Line(3) = {7, 3};
Line(4) = {3, 4};
Line(5) = {4, 1};
Line(6) = {1, 5};
Line(7) = {5, 6};
Line(8) = {6, 2};
Line(9) = {2, 1};
Line(10) = {2, 3};
Line(11) = {8, 5};
Line(12) = {7, 6};
Line Loop(13) = {4, 1, 2, 3};
Plane Surface(14) = {13};
Line Loop(15) = {11, 7, -12, -2};
Plane Surface(16) = {15};
Line Loop(17) = {5, -9, 10, 4};
Plane Surface(18) = {17};
Line Loop(19) = {1, 11, -6, -5};
Plane Surface(20) = {19};
Line Loop(21) = {10, -3, 12, 8};
Plane Surface(22) = {21};
Line Loop(23) = {8, 9, 6, 7};
Plane Surface(24) = {23};
Surface Loop(25) = {14, 18, 20, 16, 24, 22};
Volume(26) = {25};

Physical Surface(3) = {14, 20, 22, 24};
Physical Surface(2) = {18};
Physical Surface(1) = {16};
Physical Volume(0) = {26};
