%Del A
v1 = [1;2;0;1;5;1];
v2 = [1;1;0;1;3;1];
v3 = [1;2;1;2;1;2];
v4 = [2;6;3;5;2;5];

Aa = [v1 v2 v3 v4];
rref(Aa)

%Ser at de 3 første søylene i A er lin.uavhengig.
%Ergo har vi at dim(A) = 3
%Og vi kan bruke v1,v2,v3 til å lage en basis for H.

%Del B
%Bruker de 3 lin.uavhengige vektorene fra del A.
%For å få en invertibel (6x6)-matrise som avhbilder H på L6,d
%kan vi bare legge på enhetsvektorene opp til e6.
%Altså får vi: 

e4 = [0;0;0;1;0;0];
e5 = [0;0;0;0;1;0];
e6 = [0;0;0;0;0;1];
A = [v1 v2 v3 e4 e5 e6];

%Tester at det virker
rref(A);
Ainv = inv(A);
Ainv*v1
%Ser at det gjør det