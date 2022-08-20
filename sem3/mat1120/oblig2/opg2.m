%DEL 1
v1 = [1;2;1];
v2 = [2;1;1];
v3 = [3;3;1];
lambda1 = 3;
lambda2 = 6;
lambda3 = 9;
P = [v1 v2 v3];
D = [lambda1, 0,0; 0,lambda2,0 ; 0,0,lambda3];

A= P*D*inv(P)



%DEL2
C = [1,2,3,1;2,1,3,1;1,1,1,1];
rref(C)
%Gir oss [1,0,0.6667 ; 0,1,0.6667 ; 0,0,1,0.3333
%Ergo har vi C1 = 2/3, C2 = 2/3, C3 = -1/3
%og løsningen finnes på papir