a1=[1;2;3;4];
a2=[5;1;2;3];
a3=[1;2;9;1];
a4=[5;1;-4;6;];
A_start= [a1 a2 a3 a4];
rref(A_start);

%Kj�rer programmet og ser at de 3 f�rste vektorene er line�rt uavhengige.
%Sl�yfer derfor a4 og legger til I4 for � kunne utspenne R4 senere.

A_new = [a1 a2 a3 eye(4)];
%Hvis jeg s� radreduserer denne f�r jeg B bakerst.
A_new_rref= rref(A_new);

%Trekker ut B av A_new_rref
B = [A_new_rref(:,4:7)]

%Tester at B*A =rref(A)
Test1 = B*A_start;
Test2 = rref(A_start);
%Ser at det fungerer.

%Skal s� bruke B og legge til en s�yle i A for � utspenne R4
%Velger siste s�yle i B
b4= [B(:,4)];

A_R4 = [ a1 a2 a3 b4]
%Skal n� ha vektorer som utspenner R4. 
%Dobbeltsjekker ved � radredusere og ser om jeg f�r I4 som er m�let.
A_R4_rref = rref(A_R4);

%Kj�reeksempel:
%>>opg6
%B =
%
%         0    1.6667   -0.3333   -0.3333
%         0   -2.2000    0.4000    0.8000
%         0   -0.0667    0.1333   -0.0667
%    1.0000    9.4000   -1.8000   -3.6000
%
%
%A_R4 =
%
%    1.0000    5.0000    1.0000   -0.3333
%    2.0000    1.0000    2.0000    0.8000
%    3.0000    2.0000    9.0000   -0.0667
%   4.0000    3.0000    1.0000   -3.6000