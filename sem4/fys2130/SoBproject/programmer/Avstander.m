function AB = Avstander(a,s,h)

%Genererer en matrise AB,
%med avstanden fra punkt A til kilde 1-11 i rad 1,
% og avstanden fra punkt B til kilde 1-11 i rad 2

AB = zeros(11,2);

%A distances
for i=1:5
    AB(6-i,1) = sqrt((i/10*s)*(i/10*s)+(h)*(h));
end
AB(6,1)  = h;
AB(7,1)  = AB(5,1);
AB(8,1)  = AB(4,1);
AB(9,1)  = AB(3,1);
AB(10,1) = AB(2,1);
AB(11,1) = AB(1,1);

%B distances
for i=1:11
   b = s/2 - (i - 1)*s/10;
   c = a - b;
   AB(12-i,2) = sqrt((c)*(c)+(h)*(h));
end
end