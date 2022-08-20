%prelab 2

current = [5,7,10,12,15]; %Målt strøm
voltage = [1,1,1,1,1]; %Målt spenning

n=1; %Hvilken grad av polynom vi ønsker
p=polyfit(current,voltage,n); %n+1 koeffisienter i poly som passer best.
U_poly = polyval(p,current); % Spennings verdier for målt strøm
pol = p(1).*current + p(2); %Funksjonen vi lager oss med polyfit/polyval
plot(current,voltage,'r*') %Plotter målte verdier
hold on
plot(current,pol) %Plotter tilnærmingen med polyfit/polyval
%plot(x,U_poly)

%2.2
R = [500,700,1000,1200,1500]; %Liste med R på prelaben
R_poly = U_poly./current;
r_int = R_poly - R;%Indre resistans

Avg_r_int = mean(r_int) %Gjennomsnitlig indre resistans
SD_r_int = std(r_int)   %Standardavvik

