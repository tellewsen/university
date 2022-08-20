[x,y] = meshgrid(-2:.1:2);
z1 = 3.*y.^2 + 2.*x.*y;
z2 = 2.*x + 8.*y - 5;
surf(x,y,z1)
hold on
surf(x,y,z2)