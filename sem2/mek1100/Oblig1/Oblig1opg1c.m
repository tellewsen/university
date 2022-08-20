[tscaled,xscaled,yscaled]= kast(0,0,10,pi./5);
plot(xscaled,yscaled,'r')
axis([0 1 0 1])
hold on

[tscaled,xscaled,yscaled]= kast(0,0,10,pi./4);
plot(xscaled,yscaled,'g')
axis([0 1 0 1])

[tscaled,xscaled,yscaled]= kast(0,0,10,pi./3);
plot(xscaled,yscaled,'b')
axis([0 1 0 1])
legend('kast1','kast2','kast3')
xlabel('x[m]')
ylabel('y[m]')
