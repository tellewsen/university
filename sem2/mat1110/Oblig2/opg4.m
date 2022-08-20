a1=[1;5;7;6;3];
a2=[2;1;7;5;1];
a3=[3;4;1;2;3];
A= [a1 a2 a3];
rref(A)


%For at vektorene a1,a2,a3 skal være lineært uavhengige må hver søyle være
%pivotsøyle. Når vi kjører programmet ser vi at dette stemmer. Altså er
%vektorene a1,a2,a3 lineært uavhengige.

%Kjøreeksempel:
%>>opg4
%
%ans =
%
%     1     0     0
%     0     1     0
%     0     0     1
%     0     0     0
%     0     0     0