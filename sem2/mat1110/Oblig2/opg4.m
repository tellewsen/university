a1=[1;5;7;6;3];
a2=[2;1;7;5;1];
a3=[3;4;1;2;3];
A= [a1 a2 a3];
rref(A)


%For at vektorene a1,a2,a3 skal v�re line�rt uavhengige m� hver s�yle v�re
%pivots�yle. N�r vi kj�rer programmet ser vi at dette stemmer. Alts� er
%vektorene a1,a2,a3 line�rt uavhengige.

%Kj�reeksempel:
%>>opg4
%
%ans =
%
%     1     0     0
%     0     1     0
%     0     0     1
%     0     0     0
%     0     0     0