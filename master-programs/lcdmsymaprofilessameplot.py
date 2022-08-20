import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif'}) # This is for Latex writing
# import seaborn as sns
# sns.set()

#Fig resolution
width = 10.24
height = 7.68

#Limits
rlim = [1e-1,1e1]

SymAlist = [1,6,11]
labellist = ['vanilla','cool star','cool star sn']

#simtype = ['vanilla','cool star','cool star SN','DM']
#labels    = [r'$\Lambda$CDM'    ,r'SymA'     , r'SymB'     , r'SymC'     ,r'SymD',]
#testlab = ['Lv','Lc','Lcs','Av','Ac','Acs','Bv','Bc','Bcs','Cv','Cc','Ccs','Dv','Dc','Dcs']
# testlab =  [r'$\Lambda$CDM vanilla','Symmetron A vanilla','Bv','Cv','Dv',
#             r'$\Lambda$CDM cool star','Symmetron A cool star','Bc','Cc','Dc',
#             r'$\Lambda$CDM cool star sn','Symmetron A cool star sn','Bcs','Ccs','Dcs']
# lcdmsymaonly = [0,5,10,1,6,11]
# colors = ['b-','g-','r-','c-','m-','b--','g--','r--','c--','m--','b-.','g-.','r-.','c-.','m-.']


#folder1 = 'Figures/profs012/'
#folder2 = 'Figures/profs135up/'
#folder1 = 'Figures/profs012-morethan100/'
#folder1 = 'Figures/profs135up-morethan100/'
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

for k in range(len(folders)):
    #Loading files
    R = np.loadtxt('%sprofsaveR.txt'%folders[k])
    T = np.loadtxt('%sprofsaveT.txt'%folders[k])
    D = np.loadtxt('%sprofsaveD.txt'%folders[k])
    #OD = np.loadtxt('%sprofsaveOD.txt'%folders[k])
    DM = np.loadtxt('%sprofsaveDM.txt'%folders[k])
    #DMOD = np.loadtxt('%sprofsaveDMOD.txt'%folders[k])
    Dstd  = np.loadtxt('%sprofsaveDstd.txt'%folders[k])


    #####Gas density####
    ax1 = fig.add_subplot(321+k) #subplot over GRAVITY
    if k==0:
        ax1.set_title(r'$M\leq 10^{12}M_\odot/h$')
    if k==1:
        ax1.set_title(r'$M\geq 10^{13.5}M_\odot/h$')

    for j,i in enumerate(SymAlist): #GRAVITY
        plt.plot(R,D[i]/D[i-1]-1,label=labellist[j])
    #plt.legend(loc='best')

    ax1.set_xlim(rlim)
    ax1.set_xscale('log')

    #####Temperature####
    ax2 = fig.add_subplot(323+k,sharex=ax1) #subplot over GRAVITY
    for j,i in enumerate(SymAlist): #GRAVITY
        plt.plot(R,T[i]/T[i-1]-1,label=labellist[j])

    #plt.legend(loc='best')


    #####DM density####
    ax3 = fig.add_subplot(325+k,sharex=ax1) #subplot over GRAVITY
    for j,i in enumerate(SymAlist): #GRAVITY
        plt.plot(R,DM[i]/DM[i-1]-1,label=labellist[j])

    if k==0:
        plt.legend(loc='right')




    ax1.grid(color='gray',linestyle='-',linewidth=1,alpha = 0.25)
    ax2.grid(color='gray',linestyle='-',linewidth=1,alpha = 0.25)
    ax3.grid(color='gray',linestyle='-',linewidth=1,alpha = 0.25)

    if k==0:
        ax1.set_ylabel(r'$\frac{\rho_{Sym A}}{\rho_{\Lambda CDM}}-1$',labelpad=10,size='xx-large')
        ax2.set_ylabel(r'$\frac{T_{Sym A}}{T_{\Lambda CDM}}-1$',labelpad=10,size='xx-large')
        ax3.set_ylabel(r'$\frac{\rho_{DM,Sym A}}{\rho_{DM,\Lambda CDM}}-1$',labelpad=10,size='xx-large')


fig.subplots_adjust(hspace=0)
for i in [1,2,4,5]:
    plt.setp(fig.axes[i].get_xticklabels(), visible=False)

figax.set_xlabel(r'$R/R_{200c}$',size='xx-large')
plt.savefig('ProfALL.png')
