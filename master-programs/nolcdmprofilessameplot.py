import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif'}) # This is for Latex writing
# import seaborn as sns
# sns.set()

labels = [r'$\Lambda$CDM vanilla','A vanilla','B vanilla','C vanilla','D vanilla',
r'$\Lambda$CDM cool star','A cool star','B cool star','C cool star ','D cool star',
r'$\Lambda$CDM cool star sn','A cool star sn','B cool star sn','C cool star sn ','D cool star sn']

colors = ['b-','g-','r-','c-','m-','b--','g--','r--','c--','m--','b-.','g-.','r-.','c-.','m-.']

#Fig resolution
width = 10.24
height = 7.68

#Limits
rlim = [1e-1,1e1]

####################RATIOS BETWEEN SAME HYDRO#################
nolcdmprofs = [2,3,4,7,8,9,12,13,14]
reference   = [1,1,1,6,6,6,11,11,11] # reference sims for ratios


#folder1 = 'Figures/profs012/'
#folder2 = 'Figures/profs135up/'
#folder1 = 'Figures/profs012-morethan100/'
#folder2 = 'Figures/profs135up-morethan100/'

folder1 = 'Figures/profs012withstars/'
folder2 = 'Figures/profs135upwithstars/'

folders = (folder1,folder2)

#Prepare figure
fig = plt.figure(17)#,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111)
# Turn off axis lines and ticks of the big subplot
figax.spines['top'].set_color('none')
figax.spines['bottom'].set_color('none')
figax.spines['left'].set_color('none')
figax.spines['right'].set_color('none')
figax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

#loop over folders
for k in range(len(folders)):
    #load files
    R = np.loadtxt('%sprofsaveR.txt'%folders[k])
    T = np.loadtxt('%sprofsaveT.txt'%folders[k])
    D = np.loadtxt('%sprofsaveD.txt'%folders[k])
    #OD = np.loadtxt('%sprofsaveOD.txt'%folders[k])
    DM = np.loadtxt('%sprofsaveDM.txt'%folders[k])
    #DMOD = np.loadtxt('%sprofsaveDMOD.txt'%folders[k])


    #####Gas density####
    ax1 = fig.add_subplot(321+k) #subplot over GRAVITY
    if k==0:
        ax1.set_title(r'$M\leq 10^{12}M_\odot/h$')
    if k==1:
        ax1.set_title(r'$M\geq 10^{13.5}M_\odot/h$')

    for j,i in enumerate(nolcdmprofs): #GRAVITY
        plt.plot(R,D[i]/D[reference[j]]-1,colors[i],label=labels[i])
    #plt.legend(loc='best')

    #####Temperature####
    ax2 = fig.add_subplot(323+k,sharex=ax1) #subplot over GRAVITY
    for j,i in enumerate(nolcdmprofs): #GRAVITY
        plt.plot(R,T[i]/T[reference[j]]-1,colors[i],label=labels[i])
    #plt.legend(loc='best')

    #####DM density####
    ax3 = fig.add_subplot(325+k,sharex=ax1) #subplot over GRAVITY
    for j,i in enumerate(nolcdmprofs): #GRAVITY
        plt.plot(R,DM[i]/DM[reference[j]]-1,colors[i],label=labels[i])

    if k==0:
        plt.legend(loc='right')

    ax1.set_xlim(rlim)
    ax1.set_xscale('log')

    ax1.grid(color='gray',linestyle='-',linewidth=1,alpha = 0.25)
    ax2.grid(color='gray',linestyle='-',linewidth=1,alpha = 0.25)
    ax3.grid(color='gray',linestyle='-',linewidth=1,alpha = 0.25)
        # ax1.set_ylabel(r'$\frac{\rho_{Sym A}}{\rho_{\Lambda CDM}}-1$',labelpad=10,size='xx-large')
        # ax2.set_ylabel(r'$\frac{T_{Sym A}}{T_{\Lambda CDM}}-1$',labelpad=10,size='xx-large')
        # ax3.set_ylabel(r'$\frac{\rho_{DM,Sym A}}{\rho_{DM,\Lambda CDM}}-1$',labelpad=10,size='xx-large')
    if k==0:
        ax1.set_ylabel(r'$\frac{\rho }{\rho_{Sym A}}-1$',labelpad=10,size='xx-large')
        ax2.set_ylabel(r'$\frac{T}{T_{Sym A}}-1$',labelpad=10,size='xx-large')
        ax3.set_ylabel(r'$\frac{\rho_{DM}}{\rho_{DM,Sym A}}-1$',labelpad=10,size='xx-large')

fig.subplots_adjust(hspace=0)
for i in [1,2,4,5]:
    plt.setp(fig.axes[i].get_xticklabels(), visible=False)
figax.set_xlabel(r'$R/R_{200c}$',size='xx-large')
plt.savefig('ProfSymcomp.png')
