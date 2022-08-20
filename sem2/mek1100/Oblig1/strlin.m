[x,y,psi] = streamfun(5);
figure(1)
contour(x,y,psi)
axis equal
legend('n=5')
xlabel('x[m]')
ylabel('y[m]')

[x1,y1,psi1] = streamfun(30);
figure(2)
contour(x1,y1,psi1)
axis equal
legend('n=30')
xlabel('x[m]')
ylabel('y[m]')