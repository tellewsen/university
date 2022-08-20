%Setter inn grensene på omtrent samme metode som i opg 7.

triplequad(@(x,y,z)1.*(x<=(5-2.*y-z)/3).*(x>=0).*(y<=(5-3.*x-z)/2).*(y>=0).*(z<=5-3.*x-2.*y).*(z>=0), 0,5/3, 0,5/2 ,0,5)

%Får samme svar som i opg 7. Antar dette er riktig da.

%Kjøreeksempel:
%>>opg10
%
%ans =
%
%    3.4722