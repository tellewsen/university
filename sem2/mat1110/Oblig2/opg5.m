a1=[1;2;3;4];
a2=[5;1;2;3];
a3=[1;2;9;1];
a4=[5;1;-4;6;];
A= [a1 a2 a3 a4];
rref(A)

%For at vektorene a1,a2,a3,a4 skal utspenne R4 må hver rad ha pivotelement
%når vi radreduserer matrisen A = [a1 a2 a3 a4]. Vi ser når vi kjører
%programmet at dette ikke er tilfellet. Altså utspenner ikke vektorene R4.

%>>opg5
%
%ans =
%
%     1     0     0     1
%     0     1     0     1
%     0     0     1    -1
%     0     0     0     0