import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R = np.loadtxt('profsaveR.txt')
T = np.loadtxt('profsaveT.txt')
D = np.loadtxt('profsaveD.txt')
OD = np.loadtxt('profsaveOD.txt')
DM = np.loadtxt('profsaveDM.txt')
DMOD = np.loadtxt('profsaveDMOD.txt')

#DM ONLY sim data
DMOD = np.loadtxt('profsaveDMOD.txt')
DMsims = np.loadtxt('profsaveDMsims.txt')
DMODsims = np.loadtxt('profsaveDMODsims.txt')
DM = np.concatenate((DM,DMsims),axis=0)
DMOD = np.concatenate((DMOD,DMODsims),axis=0)

simtype = ['vanilla','cool star','cool star SN','DM']
labels    = [r'$\Lambda$CDM'    ,r'SymA'     , r'SymB'     , r'SymC'     ,r'SymD',]



#NO DM ONLY SIMS
# testlab = [r'$\Lambda$CDM vanilla','Symmetron A vanilla','Bv','Cv','Dv',r'$\Lambda$CDM cool star','Symmetron A cool star','Bc','Cc','Dc',r'$\Lambda$CDM cool star sn','Symmetron A cool star sn','Bcs','Ccs','Dcs']
# lcdmsymaonly = [0,5,10,1,6,11]
# colors = ['b-','g-','r-','c-','m-','b--','g--','r--','c--','m--','b-.','g-.','r-.','c-.','m-.']

#WITH DM ONLY SIMS
testlab = [r'$\Lambda$CDM vanilla','Symmetron A vanilla','Bv','Cv','Dv',
           r'$\Lambda$CDM cool star','Symmetron A cool star','Bc','Cc','Dc',
           r'$\Lambda$CDM cool star sn','Symmetron A cool star sn','Bcs','Ccs','Dcs',
           r'$\Lambda$CDM DM','Symmetron A DM','Sym B DM','SymC DM','SymD DM']
lcdmsymaonly = [0,5,10,15,1,6,11,16]
colors = ['b-','g-','r-','c-','m-','b--','g--','r--','c--','m--','b-.','g-.','r-.','c-.','m-.','b.','g.','r-','c-','m.']


#Fig resolution
width = 10.24
height = 7.68

#Limits
rlim = [1e-1,1e1]
ylimDM = np.log10([np.min(DM)-.1*np.min(DM),np.max(DM)+.1*np.max(DM)])
ylimDMOD = [np.min(DMOD)-.1*np.min(DMOD),np.max(DMOD)+.1*np.max(DMOD)]

#DM density
fig = plt.figure(9)#,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
ax = fig.add_subplot(111) #subplot over GRAVITY
for j,i in enumerate(lcdmsymaonly): #GRAVITY
    plt.plot(R,DM[i],colors[i],label=testlab[i])

#ax.set_ylim(ylimDM)
ax.set_yscale('log')
ax.set_xlim(rlim)
ax.set_xscale('log')
plt.legend(loc='best')
figax.set_title('Average DM halo density profile',y=1.05)
figax.set_ylabel(r'$\rho$ [M$_\odot$ kpc$^{-3}$]',labelpad=10)
figax.set_xlabel(r'$R/R_{vir}$')
plt.savefig('AvgHaloProfDM.png')

# DM OverDensity
fig = plt.figure(10)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
ax = fig.add_subplot(111) #subplot over GRAVITY

for j,i in enumerate(lcdmsymaonly): #GRAVITY
        plt.plot(R,DMOD[i],colors[i],label=testlab[i])#)str(labels[i])+' '+str(simtype[j]))

ax.set_xscale('log')
ax.set_xlim(rlim)
#ax.set_ylim(ylimOD)
ax.set_yscale('log')

plt.legend(loc='best')
figax.set_title('Average DM halo overdensity profile',y=1.05)
figax.set_ylabel(r'Overdensity $\frac{\rho}{\bar{\rho}}-1$',labelpad=10)
figax.set_xlabel(r'$R/R_{vir}$')
plt.savefig('AvgHaloProfDMOD.png')
