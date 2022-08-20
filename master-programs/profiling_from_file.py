import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

folder = 'Figures/profs135up/'

R       = np.loadtxt('%sprofsaveR.txt'%folder)
T       = np.loadtxt('%sprofsaveT.txt'%folder)
D       = np.loadtxt('%sprofsaveD.txt'%folder)
OD      = np.loadtxt('%sprofsaveOD.txt'%folder)
DM      = np.loadtxt('%sprofsaveDM.txt'%folder)
DMOD    = np.loadtxt('%sprofsaveDMOD.txt'%folder)

# Tstd = np.loadtxt('%sprofsaveTstd.txt'%folder)
Dstd = np.loadtxt('%sprofsaveDstd.txt'%folder)
# ODstd = np.loadtxt('%sprofsaveODstd.txt'%folder)
# DMstd = np.loadtxt('%sprofsaveDMstd.txt'%folder)
# DMODstd = np.loadtxt('%sprofsaveDMODstd.txt'%folder)

#DMsims = np.loadtxt('profsaveDMsims.txt')
#DMODsims = np.loadtxt('profsaveDMODsims.txt')

#simtype = ['vanilla','cool star','cool star SN','DM']
#labels    = [r'$\Lambda$CDM'    ,r'SymA'     , r'SymB'     , r'SymC'     ,r'SymD',]

#No DM ONLY sims
labels = ['Lv','Av','Bv','Cv','Dv','Lc','Ac','Bc','Cc','Dc','Lcs','Acs','Bcs','Ccs','Dcs']
colors = ['b-','g-','r-','c-','m-','b--','g--','r--','c--','m--','b-.','g-.','r-.','c-.','m-.']

# #WITH DM ONLY SIMS
# labels = ['Lv','Av','Bv','Cv','Dv','Lc','Ac','Bc','Cc','Dc','Lcs','Acs','Bcs','Ccs','Dcs']
# colors = ['b-','g-','r-','c-','m-','b--','g--','r--','c--','m--','b-.','g-.','r-.','c-.','m-.']

#resolutions
width = 10.24
height = 7.68

#limits
rlim = [1e-1,1e1]
ylimD  = np.log10([np.min(D)-.1*np.min(D),np.max(D)+.1*np.max(D)])
ylimT  = np.log10([np.min(T)-.1*np.min(T),np.max(T)+.1*np.max(T)])
ylimDM = np.log10([np.min(DM)-.1*np.min(DM),np.max(DM)+.1*np.max(DM)])

ylimOD = [np.min(OD)-.1*np.min(OD),np.max(OD)+.1*np.max(OD)]
ylimDMOD = [np.min(DMOD)-.1*np.min(DMOD),np.max(DMOD)+.1*np.max(DMOD)]


"""
#gas Temperature
fig = plt.figure(6)#,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
ax = fig.add_subplot(111) #subplot over GRAVITY

for i in range(len(labels)): #GRAVITY
    plt.plot(R,T[i],colors[i],label=labels[i])#str(labels[i])+' '+str(simtype[j]))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

#plt.legend()
#ax.yscale('log')
#ax.set_ylim(ylimT)
ax.set_xlim(rlim)
ax.set_xscale('log')

figax.set_title('Average halo temperature profile',y=1.05)
figax.set_ylabel('Temperature [K]',labelpad=10)
figax.set_xlabel(r'$R/R_{vir}$')
plt.savefig('AvgHaloProfT.png')

"""

#Gas density
fig = plt.figure(8)#,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
ax = fig.add_subplot(111) #subplot over GRAVITY
for i in range(len(labels)): #GRAVITY
    plt.plot(R,D[i],colors[i],label=labels[i])
    plt.fill_between(R,D[i]+Dstd[i],D[i]-Dstd[i],facecolor=colors[i][0],alpha=0.5)
#ax.set_ylim(ylimD)
#ax.set_ylim(np.min(D-Dstd),np.max(D+Dstd))
#ax.set_yscale('log')
ax.set_xlim(rlim)
ax.set_xscale('log')
#plt.legend()
figax.set_title('Average halo density profile',y=1.05)
figax.set_ylabel(r'$\rho$ [M$_\odot$ kpc$^{-3}$]',labelpad=10)
figax.set_xlabel(r'$R/R_{vir}$')
plt.savefig('AvgHaloProfD.png')

sys.exit()


#gas OverDensity
fig = plt.figure(7)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
ax = fig.add_subplot(111) #subplot over GRAVITY

for i in range(len(labels)): #GRAVITY
        plt.plot(R,OD[i],colors[i],label=labels[i])#)str(labels[i])+' '+str(simtype[j]))

ax.set_xscale('log')
ax.set_xlim(rlim)
#ax.set_ylim(ylimOD)
ax.set_yscale('log')

#plt.legend()
figax.set_title('Average halo overdensity profile',y=1.05)
figax.set_ylabel(r'Overdensity $\frac{\rho}{\bar{\rho}}-1$',labelpad=10)
figax.set_xlabel(r'$R/R_{vir}$')
plt.savefig('AvgHaloProfOD.png')

#DM density
fig = plt.figure(9)#,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
ax = fig.add_subplot(111) #subplot over GRAVITY
for i in range(len(labels)): #GRAVITY
    plt.plot(R,DM[i],colors[i],label=labels[i])

#ax.set_ylim(ylimDM)
ax.set_yscale('log')
ax.set_xlim(rlim)
ax.set_xscale('log')
#plt.legend()
figax.set_title('Average DM halo density profile',y=1.05)
figax.set_ylabel(r'$\rho$ [M$_\odot$ kpc$^{-3}$]',labelpad=10)
figax.set_xlabel(r'$R/R_{vir}$')
plt.savefig('AvgHaloProfDM.png')

# DM OverDensity
fig = plt.figure(10)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
ax = fig.add_subplot(111) #subplot over GRAVITY

for i in range(len(labels)): #GRAVITY
        plt.plot(R,DMOD[i],colors[i],label=labels[i])#)str(labels[i])+' '+str(simtype[j]))

ax.set_xscale('log')
ax.set_xlim(rlim)
#ax.set_ylim(ylimOD)
ax.set_yscale('log')

#plt.legend()
figax.set_title('Average DM halo overdensity profile',y=1.05)
figax.set_ylabel(r'Overdensity $\frac{\rho}{\bar{\rho}}-1$',labelpad=10)
figax.set_xlabel(r'$R/R_{vir}$')
plt.savefig('AvgHaloProfDMOD.png')

####RATIOS####





















####OLD STUFF###
"""
#Temperature hydro vs gravity
fig = plt.figure(0,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
# Turn off axis lines and ticks of the big subplot
figax.spines['top'].set_color('none')
figax.spines['bottom'].set_color('none')
figax.spines['left'].set_color('none')
figax.spines['right'].set_color('none')
figax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
for i in range(3): #hydro
    ax = fig.add_subplot(221+i) #subplot over GRAVITY
    ax.set_title('%s'%simtype[i])
    for j in range(5): #gravity
        #plt.semilogy(R,T[i+3*j],label=labels[j])
        plt.plot(R,np.log10(T[i+3*j]),label=labels[j])
    ax.set_ylim(ylimT)
    #ax.set_yscale('symlog')
    ax.set_xlim(rlim)
    #ax.set_xscale('log')
    plt.legend()

figax.set_title('Average halo temperature profile',y=1.05)
figax.set_ylabel('Temperature (log) [K]',labelpad=10)
figax.set_xlabel(r'$R/R_{vir}$')
plt.savefig('AvgHaloProfT_HvG.png')
#plt.show()

#overdenstiy hydro vs gravity
fig = plt.figure(1,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
# Turn off axis lines and ticks of the big subplot
figax.spines['top'].set_color('none')
figax.spines['bottom'].set_color('none')
figax.spines['left'].set_color('none')
figax.spines['right'].set_color('none')
figax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
for i in range(3): #hyrd
    ax = fig.add_subplot(221+i) #subplot over GRAVITY
    ax.set_title('%s'%simtype[i])
    for j in range(5): #gravity
        plt.plot(R,OD[i+3*j],label=labels[j])
    ax.set_ylim(ylimOD)
    #ax.set_yscale('symlog')
    ax.set_xlim(rlim)
    #ax.set_xscale('log')
    plt.legend()
figax.set_title('Average halo overdensity profile',y=1.05)
figax.set_ylabel(r'Overdensity $\frac{\rho}{\bar{\rho}} -1$',labelpad=10)
figax.set_xlabel(r'$R/R_{vir}$')
plt.savefig('AvgHaloProfOD_HvG.png')
#plt.show()

#density hydro vs gravity
fig = plt.figure(2,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
# Turn off axis lines and ticks of the big subplot
figax.spines['top'].set_color('none')
figax.spines['bottom'].set_color('none')
figax.spines['left'].set_color('none')
figax.spines['right'].set_color('none')
figax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
for i in range(3): #hyrdo
    ax = fig.add_subplot(221+i) #subplot over GRAVITY
    ax.set_title('%s'%simtype[i])
    for j in range(5): #gravity
        plt.plot(R,np.log10(D[i+3*j]),label=labels[j])
    ax.set_ylim(ylimD)
    #ax.set_yscale('symlog')
    ax.set_xlim(rlim)
    #ax.set_xscale('log')
    plt.legend()
figax.set_title('Average halo density profile',y=1.05)
figax.set_ylabel(r'$\rho$ [M$_\odot$ kpc$^{-3}$] (log)',labelpad=10)
figax.set_xlabel(r'$R/R_{vir}$')
plt.savefig('AvgHaloProfD_HvG.png')

#Temperature gravity v hydro
fig = plt.figure(3,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
# Turn off axis lines and ticks of the big subplot
figax.spines['top'].set_color('none')
figax.spines['bottom'].set_color('none')
figax.spines['left'].set_color('none')
figax.spines['right'].set_color('none')
figax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
for i in range(5): #GRAVITY
    ax = fig.add_subplot(321+i) #subplot over GRAVITY
    ax.set_title('%s'%labels[i])
    for j in range(3): #HYDRO
        #plt.semilogy(R,T[i+3*j],label=labels[i+3*j])
        plt.plot(R,np.log10(T[j+3*i]),label=simtype[j])
    ax.set_ylim(ylimT)
    #ax.set_yscale('symlog')
    ax.set_xlim(rlim)
    #ax.set_xscale('log')
    plt.legend()
figax.set_title('Average halo temperature profile',y=1.05)
figax.set_ylabel('Temperature (log) [K]',labelpad=10)
figax.set_xlabel(r'$R/R_{vir}$')
plt.savefig('AvgHaloProfT_GvH.png')
#plt.show()

#overdnsity gravity v hydro
fig = plt.figure(4,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
# Turn off axis lines and ticks of the big subplot
figax.spines['top'].set_color('none')
figax.spines['bottom'].set_color('none')
figax.spines['left'].set_color('none')
figax.spines['right'].set_color('none')
figax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
for i in range(5): #GRAVITY
    ax = fig.add_subplot(321+i) #subplot over GRAVITY
    ax.set_title('%s'%labels[i])
    for j in range(3): #HYDRO
        #plt.semilogy(R,T[i+3*j],label=labels[i+3*j])
        plt.plot(R,OD[j+3*i],label=simtype[j])
    ax.set_ylim(ylimOD)
    #ax.set_yscale('symlog')
    ax.set_xlim(rlim)
    #ax.set_xscale('log')
    plt.legend()
figax.set_title('Average halo overdensity profile',y=1.05)
figax.set_ylabel(r'Overdensity $\frac{\rho}{\bar{\rho}}-1$',labelpad=10)
figax.set_xlabel(r'$R/R_{vir}$')
plt.savefig('AvgHaloProfOD_GvH.png')
#plt.show()

#dnsity gravity v hydro
fig = plt.figure(5,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
# Turn off axis lines and ticks of the big subplot
figax.spines['top'].set_color('none')
figax.spines['bottom'].set_color('none')
figax.spines['left'].set_color('none')
figax.spines['right'].set_color('none')
figax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
for i in range(5): #GRAVITY
    ax = fig.add_subplot(321+i) #subplot over GRAVITY
    ax.set_title('%s'%labels[i])
    for j in range(3): #HYDRO
        #plt.semilogy(R,T[i+3*j],label=labels[i+3*j])
        plt.plot(R,np.log10(D[j+3*i]),label=simtype[j])
    ax.set_ylim(ylimD)
    #ax.set_yscale('symlog')
    ax.set_xlim(rlim)
    #ax.set_xscale('log')
    plt.legend()
figax.set_title('Average halo density profile',y=1.05)
figax.set_ylabel(r'$\rho$ [M$_\odot$ kpc$^{-3}$] (log)',labelpad=10)
figax.set_xlabel(r'$R/R_{vir}$')fig = plt.figure(5,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
# Turn off axis lines and ticks of the big subplot
figax.spines['top'].set_color('none')
figax.spines['bottom'].set_color('none')
figax.spines['left'].set_color('none')
figax.spines['right'].set_color('none')
figax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
for i in range(5): #GRAVITY
    ax = fig.add_subplot(321+i) #subplot over GRAVITY
    ax.set_title('%s'%labels[i])
    for j in range(3): #HYDRO
        #plt.semilogy(R,T[i+3*j],label=labels[i+3*j])
        plt.plot(R,np.log10(D[j+3*i]),label=simtype[j])
    ax.set_ylim(ylimD)
    #ax.set_yscale('symlog')
    ax.set_xlim(rlim)
    #ax.set_xscale('log')
    plt.legend()
figax.set_title('Average halo density profile',y=1.05)
figax.set_ylabel(r'$\rho$ [M$_\odot$ kpc$^{-3}$] (log)',labelpad=10)
figax.set_xlabel(r'$R/R_{vir}$')
plt.savefig('AvgHaloProfD_GvH.png')
plt.savefig('AvgHaloProfD_GvH.png')

"""
