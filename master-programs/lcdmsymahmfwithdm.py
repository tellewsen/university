import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif'}) # This is for Latex writing
# import seaborn as sns
# sns.set()


from matplotlib.font_manager import FontProperties

#Halo mass histograms Same hydroeffect in same plot.

simtypes   = ['LCDM_vanilla'     ,'SymA_vanilla'     ,'SymB_vanilla'     ,'SymC_vanilla'     ,'SymD_vanilla',
              'LCDM_cool_star'   ,'SymA_cool_star'   ,'SymB_cool_star'   ,'SymC_cool_star'   ,'SymD_cool_star',
              'LCDM_cool_star_SN','SymA_cool_star_SN','SymB_cool_star_SN','SymC_cool_star_SN','SymD_cool_star_SN',
              'LCDM_DM'          ,'SymA_DM'          ,'SymB_DM'          ,'SymC_DM'          ,'SymD_DM']

simtype = ['vanilla','cool star','cool star SN','DM']
labels    = [r'$\Lambda$CDM'    ,r'SymA'     , r'SymB'     , r'SymC'     ,r'SymD']

testlab = [r'$\Lambda$CDM vanilla','Symmetron A vanilla','Bv','Cv','Dv'
          ,r'$\Lambda$CDM cool star','Symmetron A cool star','Bc','Cc','Dc'
          ,r'$\Lambda$CDM cool star sn','Symmetron A cool star sn','Bcs','Ccs','Dcs'
          ,r'$\Lambda$CDM DM','Symmetron A DM','BDM','CDM','DDM']

lcdmsymaonly = [0,5,10,15,1,6,11,16]
colors = ['b-','b--','b-','b-','b-','g-','g--','g--','g--','g--','r-','r--','r-.','r-.','r-','c-','c--','c:','c:','c:']


Volume = 64.*64.*64. #( h(1+z)Mpc ) ^3
binrange = (11,14) #bin range in log space

try:
    num_bins = int(sys.argv[1])#200 #Number of bins
except:
    print 'Correct usage: python thisfile.py num_bins'
    sys.exit(1)


#window size
width = 10.24 #1024 pixels
height = 7.68 #768 pixels

##############################RATIO PART##############################

###Make the standard to take ratios from###
LCDMmvir = np.loadtxt('Halos/LCDM_vanilla.txt', usecols=8,unpack=True)
LCDMmvir = np.log10(LCDMmvir)
LCDMcounts,LCDMbins = np.histogram(LCDMmvir,num_bins,range=binrange)
LCDMbins = LCDMbins[:-1]
LCDMbins = LCDMbins + (LCDMbins[1]-LCDMbins[0])/2
LCDMcounts = LCDMcounts/Volume
for i in range(len(LCDMcounts)):
    LCDMcounts[i] = np.sum(LCDMcounts[i:])
#############################################


#Plot ratios of each to LCDM vanilla over gravity
fig = plt.figure(4,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111)
# Turn off axis lines and ticks of the big subplot
figax.spines['top'].set_color('none')
figax.spines['bottom'].set_color('none')
figax.spines['left'].set_color('none')
figax.spines['right'].set_color('none')
figax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

#####All HMFs in one plot, then ratio in another plot

#HMFs
ax = fig.add_subplot(211)
for i,j in enumerate(lcdmsymaonly):
    mvir = np.loadtxt('Halos/%s.txt'%simtypes[j], usecols=8,unpack=True)
    mvir = np.log10(mvir)
    counts,bins = np.histogram(mvir,num_bins,range=binrange)
    bins = bins[:-1] #Remove last edge

    #shift bins half a bin to the right
    bins += (bins[1]-bins[0])/2.

    # turn counts into counts per comoving Mpc
    counts = counts/Volume
    #Turn into sum of larger than
    for k in range(len(counts)):
        counts[k] = np.sum(counts[k:])

    plt.plot(bins,counts,colors[j],label=testlab[j])#,'o')
ax.set_xlim([11,14])
ax.set_ylabel(r'n(>M) $[h^3 Mpc^{-3}]$',labelpad=10,size='x-large')
plt.yscale('log')
plt.legend(loc=3)

#RATIOS
ax = fig.add_subplot(212)
for i,j in enumerate(lcdmsymaonly):
    mvir = np.loadtxt('Halos/%s.txt'%simtypes[j], usecols=8,unpack=True)
    mvir = np.log10(mvir)
    counts,bins = np.histogram(mvir,num_bins,range=binrange)
    bins = bins[:-1] #Remove last edge

    #shift bins half a bin to the right
    bins += (bins[1]-bins[0])/2.

    # turn counts into counts per comoving Mpc
    counts = counts/Volume
    #Turn into sum of larger than
    for k in range(len(counts)):
        counts[k] = np.sum(counts[k:])

    plt.plot(bins,counts/LCDMcounts-1.,colors[j],label=testlab[j])#,'o')

#ax.set_xlim([10,15])
ax.set_ylabel(r'$\frac{n(>M)}{n_{\Lambda CDM}(>M)} -1 $',labelpad=10,size='xx-large')
ax.grid(color='gray',linestyle='-',linewidth=.2,alpha=1)
#figax.set_title(r'Ratio HMF to $\Lambda$CDM vanilla HMF, $N_{bins}$=%i'%num_bins,y=1.05)
figax.set_xlabel(r'M($M_\odot h^{-1}$)(log)',size='x-large')
ax.set_xlim([11,14])
plt.savefig('HMF_symaonly_%i.png'%num_bins)
plt.show()
