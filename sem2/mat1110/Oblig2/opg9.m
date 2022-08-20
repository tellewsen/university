%Bruker kalkulator for � finne grensene
xmin  = 0;
xmid  = 0.191717;
xmax  = 0.567143;
ymin  = 0;
ymax1 = 2*exp(5*xmid);
ymax2 = exp(xmax);


%Finner arealet ved � regne ut omr�det under grafen frem til de �verste krysser
%hverandre og s� etter de har krysset. Trekker s� fra omr�det
%under exp(x)
Area1 = dblquad(@(x,y)(x.^2.*y.*sin(x.*y))  .*(y<=2.*exp(5.*x)).*(y>=0),xmin,xmid,ymin,ymax1);
Area2 = dblquad(@(x,y)(x.^2.*y.*sin(x.*y))  .*(y<=1./x)        .*(y>=0),xmid,xmax,ymin,ymax1);
Area3 = dblquad(@(x,y)(x.^2.*y.*sin(x.*y))  .*(y<=exp(x))      .*(y>=0),xmin,xmax,ymin,ymax2);

Area= Area1 + Area2 - Area3

%Her setter jeg inn grensene slik at hele omr�det blir regnet ut ved hjelp
%av 1 linje med kode.
testArea = dblquad(@(x,y)(x.^2.*y.*sin(x.*y)).*(y<=2*exp(5.*x)).*(y>=exp(x)).*(y<=1./x),0,1,0,1)

%Denne b�r fungere bedre enn l�sningen over, men usikker p� hva som er
%korrekt


%Kj�reeksempel:
%>> opg9
%
%Area =
%
%    0.0897
%
%
%testArea =
%
%    0.0873