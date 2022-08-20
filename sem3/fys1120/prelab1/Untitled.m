%prelab 2

current = [5,7,10,12,15]; %M�lt str�m
voltage = [1,1,1,1,1]; %M�lt spenning

n=1; %Hvilken grad av polynom vi �nsker
p=polyfit(current,voltage,n); %n+1 koeffisienter i poly som passer best.
U_poly = polyval(p,current); % Spennings verdier for m�lt str�m
pol = p(1).*current + p(2); %Funksjonen vi lager oss med polyfit/polyval
plot(current,voltage,'r*') %Plotter m�lte verdier
hold on
plot(current,pol) %Plotter tiln�rmingen med polyfit/polyval
%plot(x,U_poly)

%2.2
R = [500,700,1000,1200,1500]; %Liste med R p� prelaben
R_poly = U_poly./current;
r_int = R_poly - R;%Indre resistans

Avg_r_int = mean(r_int) %Gjennomsnitlig indre resistans
SD_r_int = std(r_int)   %Standardavvik

