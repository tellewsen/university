%x1=[1,0,2,1]
%x2=[0,1,1,3]

b1 = [1; 0; 0; 0];
b2 = [0; 0; 0; 1];
b3 = [8; 4; 2; 1];
b4 = [1; 3; 9; 27];
A = [b1, b2, b3, b4];

rref(A)

%Velger verdier for x1 og x2 slik at vi f�r 4 vektorer ut.
%Setter disse som s�yler i en matrise A.
%Radreduserer denne og ser at vi f�r ut I4.
%Ergo er vektorene line�rt uavhengige og ogs� en basis for P2,3.
