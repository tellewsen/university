import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif'}) # This is for Latex writing
# import seaborn as sns
# sns.set()

# folder = 'Figures/profs135up/'
folder = 'Figures/profs012/'

R = np.loadtxt('%sprofsaveR.txt'%folder)
T = np.loadtxt('%sprofsaveT.txt'%folder)
D = np.loadtxt('%sprofsaveD.txt'%folder)
OD = np.loadtxt('%sprofsaveOD.txt'%folder)
DM = np.loadtxt('%sprofsaveDM.txt'%folder)
DMOD = np.loadtxt('%sprofsaveDMOD.txt'%folder)

#simtype = ['vanilla','cool star','cool star SN','DM']
#labels    = [r'$\Lambda$CDM'    ,r'SymA'     , r'SymB'     , r'SymC'     ,r'SymD',]

labels = [r'$\Lambda$CDM vanilla','Sym A vanilla','Sym B vanilla','Sym C vanilla','Sym D vanilla',
r'$\Lambda$CDM cool star','Sym A cool star','Sym B cool star','Sym C cool star ','Sym D cool star',
r'$\Lambda$CDM cool star sn','Sym A cool star sn','Sym B cool star sn','Sym C cool star sn ','Sym D cool star sn']
nolcdmprofs = [1,2,3,4,6,7,8,9,11,12,13,14]
#labellist = ['vanilla','cool star','cool star sn']


colors = ['b-','g-','r-','c-','m-','b--','g--','r--','c--','m--','b-.','g-.','r-.','c-.','m-.']


#Fig resolution
width = 10.24
height = 7.68

#Limits
rlim = [1e-1,1e1]
ylimD  = np.log10([np.min(D)-.1*np.min(D),np.max(D)+.1*np.max(D)])
ylimT  = np.log10([np.min(T)-.1*np.min(T),np.max(T)+.1*np.max(T)])
ylimDM = np.log10([np.min(DM)-.1*np.min(DM),np.max(DM)+.1*np.max(DM)])
ylimOD = [np.min(OD)-.1*np.min(OD),np.max(OD)+.1*np.max(OD)]
ylimDMOD = [np.min(DMOD)-.1*np.min(DMOD),np.max(DMOD)+.1*np.max(DMOD)]

"""
#Gas Temperature
fig = plt.figure(6)#,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
ax = fig.add_subplot(111) #subplot over GRAVITY

for j,i in enumerate(nolcdmprofs): #GRAVITY
    plt.plot(R,T[i],colors[i],label=labels[i])#str(labels[i])+' '+str(simtype[j]))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

#Formatted such that vanilla has whole line, c has --, and cs has -.
plt.legend(loc='best')
#ax.yscale('log')
#ax.set_ylim(ylimT)
ax.set_xlim(rlim)
ax.set_xscale('log')

figax.set_title('Average halo temperature profile',y=1.05)
figax.set_ylabel('Temperature [K]',labelpad=10)
figax.set_xlabel(r'$R/R_{vir}$')
plt.savefig('AvgHaloProfT.png')

#Gas density
fig = plt.figure(8)#,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
ax = fig.add_subplot(111) #subplot over GRAVITY
for j,i in enumerate(nolcdmprofs): #GRAVITY
    plt.plot(R,D[i],colors[i],label=labels[i])

#ax.set_ylim(ylimD)
ax.set_yscale('log')
ax.set_xlim(rlim)
ax.set_xscale('log')

plt.legend(loc='best')
figax.set_title('Average halo density profile',y=1.05)
figax.set_ylabel(r'$\rho$ [M$_\odot$ kpc$^{-3}$]',labelpad=10)
figax.set_xlabel(r'$R/R_{vir}$')
plt.savefig('AvgHaloProfD.png')

# gas OverDensity
fig = plt.figure(7)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
ax = fig.add_subplot(111) #subplot over GRAVITY

for j,i in enumerate(nolcdmprofs): #GRAVITY
        plt.plot(R,OD[i],colors[i],label=labels[i])#)str(labels[i])+' '+str(simtype[j]))

ax.set_xscale('log')
ax.set_xlim(rlim)
#ax.set_ylim(ylimOD)
ax.set_yscale('log')

plt.legend(loc='best')
figax.set_title('Average halo overdensity profile',y=1.05)
figax.set_ylabel(r'Overdensity $\frac{\rho}{\bar{\rho}}-1$',labelpad=10)
figax.set_xlabel(r'$R/R_{vir}$')
plt.savefig('AvgHaloProfOD.png')

#DM density
fig = plt.figure(9)#,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
ax = fig.add_subplot(111) #subplot over GRAVITY


for j,i in enumerate(nolcdmprofs): #GRAVITY
    plt.plot(R,DM[i],colors[i],label=labels[i])

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

for j,i in enumerate(nolcdmprofs): #GRAVITY
        plt.plot(R,DMOD[i],colors[i],label=labels[i])#)str(labels[i])+' '+str(simtype[j]))

ax.set_xscale('log')
ax.set_xlim(rlim)
#ax.set_ylim(ylimOD)
ax.set_yscale('log')

plt.legend(loc='best')
figax.set_title('Average DM halo overdensity profile',y=1.05)
figax.set_ylabel(r'Overdensity $\frac{\rho}{\bar{\rho}}-1$',labelpad=10)
figax.set_xlabel(r'$R/R_{vir}$')
plt.savefig('AvgHaloProfDMOD.png')
"""


####################RATIOS BETWEEN SAME HYDRO#################
nolcdmprofs = [2,3,4,7,8,9,12,13,14]
reference   = [1,1,1,6,6,6,11,11,11] # reference sims for ratios

#####Gas density####
fig = plt.figure(17)#,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
ax = fig.add_subplot(111) #subplot over GRAVITY



for j,i in enumerate(nolcdmprofs): #GRAVITY
    plt.plot(R,D[i]/D[reference[j]]-1,colors[i],label=labels[i])

#ax.set_ylim(ylimD)
#ax.set_yscale('log')
ax.set_xlim(rlim)
ax.set_xscale('log')

plt.legend(loc='best')
#figax.set_title('Average halo density profile',y=1.05)
figax.set_ylabel(r'$\rho/\rho_{Sym A}-1$',labelpad=10,size='xx-large')
figax.set_xlabel(r'$R/R_{200c}$',size='xx-large')
plt.savefig('ProfDratioHydro.png')

#####Temperature####
fig = plt.figure(18)#,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
ax = fig.add_subplot(111) #subplot over GRAVITY

for j,i in enumerate(nolcdmprofs): #GRAVITY
    plt.plot(R,T[i]/T[reference[j]]-1,colors[i],label=labels[i])

#ax.set_ylim(ylimD)
#ax.set_yscale('log')
ax.set_xlim(rlim)
ax.set_xscale('log')

plt.legend(loc='best')
#figax.set_title('Average halo temperature profile',y=1.05)
figax.set_ylabel(r'$T/T_{Sym A}-1$',labelpad=10,size='xx-large')
figax.set_xlabel(r'$R/R_{200c}$',size='xx-large')
plt.savefig('ProfTratioHydro.png')

#####DM density####
fig = plt.figure(19)#,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
ax = fig.add_subplot(111) #subplot over GRAVITY

for j,i in enumerate(nolcdmprofs): #GRAVITY
    plt.plot(R,DM[i]/DM[reference[j]]-1,colors[i],label=labels[i])

#ax.set_ylim(ylimD)
#ax.set_yscale('log')
ax.set_xlim(rlim)
ax.set_xscale('log')

plt.legend(loc='best')
#figax.set_title('Average halo DM density profile',y=1.05)
figax.set_ylabel(r'$\rho/\rho_{Sym A}-1$',labelpad=10,size='xx-large')
figax.set_xlabel(r'$R/R_{200c}$',size='xx-large')
plt.savefig('ProfDMratioHydro.png')
