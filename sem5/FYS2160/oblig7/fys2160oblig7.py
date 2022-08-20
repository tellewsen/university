# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 21:34:27 2014

@author: Andreas Ellewsen
"""
#Start of part c
Vhat = np.linspace(0.4,20,1001)
That = array([1.15,1.0,0.85])
phat = zeros([len(That),len(Vhat)])

figure(0)
hold('on')

for j in range(len(That)):
    for i in range(len(Vhat)):
        phat[j,i] = 8*That[j]/(3*Vhat[i] -1) -3/(Vhat[i]**2)
    plot(Vhat,phat[j])
legend(['That = 1.15','That = 1.0','That = 0.85'])
xlabel('Vhat')
ylabel('phat')
savefig('c.png')
#End of part c

#Start of part e
rhohat = np.linspace(0,2,1001)
figure(1)
for j in range(len(That)):
    for i in range(len(rhohat)):
        phat[j,i] = 3*rhohat[i]*That[j]/(3-rhohat[i]) - 3*rhohat[i]**2
    plot(rhohat,phat[j])
legend(['That = 1.15','That = 1.0','That = 0.85'])
xlabel('Vhat')
ylabel('rhohat')
savefig('e.png')
#End of part e

#Start of part k
That = .9
rhohat = np.linspace(0.2,2,1001)
phat = 3*rhohat*That/(3-rhohat) - 3*rhohat**2

figure(2)
subplot(3,1,1)
plot(rhohat,phat)
xlabel('rhohat')
ylabel('phat')
legend(['phat(rhohat)'])


Vhat = 1./rhohat

subplot(3,1,2)
plot(rhohat,Vhat)
xlabel('rhohat')
ylabel('Vhat')
legend(['Vhat(rhohat)'])

ghat = -3*rhohat-8./3*That*log(3/rhohat-1) + phat/rhohat

subplot(3,1,3)
plot(rhohat,ghat)
xlabel('rhohat')
ylabel('ghat')
legend(['ghat(rhohat)'])
savefig('k.png')
#To do: Pick some points and realize stuff"
#End of part k

#Start of part m

#End of part m