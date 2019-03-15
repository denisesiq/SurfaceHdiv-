

p1 = newp; Point(p1) = {0,0,0};
p2 = newp; Point(p2) = {1,0,0};
p3 = newp; Point(p3) = {1,1,0};
p4 = newp; Point(p4) = {0,1,0};


l1=newl; Line(l1) = {p1,p2};
l2=newl; Line(l2) = {p2,p3};
l3=newl; Line(l3) = {p3,p4};
l4=newl; Line(l4) = {p4,p1};

ll1 = newll; Line Loop(ll1) = {l1,l2,l3,l4};
s1 =news; Surface(s1) = {ll1};

Physical Surface("domain") = {s1};
Physical Line("boundary") = {l1,l2,l3,l4};

