import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif'}) # This is for Latex writing
# import seaborn as sns
# sns.set()
#Halo mass histograms Same hydroeffect in same plot.

simtypes   = ['LCDM_vanilla'     ,'SymA_vanilla'     ,'SymB_vanilla'     ,'SymC_vanilla'     ,'SymD_vanilla',
              'LCDM_cool_star'   ,'SymA_cool_star'   ,'SymB_cool_star'   ,'SymC_cool_star'   ,'SymD_cool_star',
              'LCDM_cool_star_SN','SymA_cool_star_SN','SymB_cool_star_SN','SymC_cool_star_SN','SymD_cool_star_SN',
              'LCDM_DM'          ,'SymA_DM'          ,'SymB_DM'          ,'SymC_DM'          ,'SymD_DM']

#simtype = ['vanilla','cool star','cool star SN','DM']
#labels    = [r'$\Lambda$CDM'    ,r'SymA'     , r'SymB'     , r'SymC'     ,r'SymD']

# testlab = [     r'$\Lambda$CDM vanilla','Symmetron A vanilla','Symmetron B vanilla','Symmetron C vanilla','Symmetron D vanilla',
#                 r'$\Lambda$CDM cool star','Symmetron A cool star','Symmetron B cool star','Symmetron C cool star','Symmetron D cool star',
#                 r'$\Lambda$CDM cool star sn','Symmetron A cool star sn','Symmetron B cool star sn','Symmetron C cool star sn','Symmetron D cool star sn']
testlab = [     r'$\Lambda$CDM vanilla','A vanilla','B vanilla','C vanilla','D vanilla',
                r'$\Lambda$CDM cool star','A cool star','B cool star','C cool star','D cool star',
                r'$\Lambda$CDM cool star sn','A cool star sn','B cool star sn','C cool star sn','D cool star sn',
                r'$\Lambda$CDM DM','A DM','B DM','C DM','D DM']


nolcdm = [1,6,11,16,2,7,12,17,3,8,13,18,4,9,14,19]
colors = ['b-','g-','r-','c-','m-','b--','g--','r--','c--','m--','b-.','g-.','r-.','c-.','m-.','b:','g:','r:','c:','m:']




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
SymAmvir = np.loadtxt('Halos/SymA_vanilla.txt', usecols=8,unpack=True)
SymAmvir = np.log10(SymAmvir)
SymAcounts,SymAbins = np.histogram(SymAmvir,num_bins,range=binrange)
SymAbins = SymAbins[:-1]
SymAbins = SymAbins + (SymAbins[1]-SymAbins[0])/2
SymAcounts = SymAcounts/Volume
for i in range(len(SymAcounts)):
    SymAcounts[i] = np.sum(SymAcounts[i:])
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
for i,j in enumerate(nolcdm):
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

ax.set_ylabel(r'n(>M) $[h^3 Mpc^{-3}]$',labelpad=10,size='x-large')
ax.set_xlim([11,14])
plt.legend(loc=3,ncol=2)
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.2)
plt.yscale('log')

#RATIOS
ax = fig.add_subplot(212)
for i,j in enumerate(nolcdm[1:]):
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

    plt.plot(bins,counts/SymAcounts -1.,colors[j],label=testlab[j])#,'o')

#ax.set_xlim([10,15])
ax.set_ylabel(r'$\frac{n(>M)}{n_{SymA}(>M)} -1 $',labelpad=10,size='xx-large')
ax.grid(color='gray',linestyle='-',linewidth=.2,alpha=1)
ax.set_xlim([11,14])
#figax.set_title(r'Ratio HMF to $\Lambda$CDM vanilla HMF, $N_{bins}$=%i'%num_bins,y=1.05)

figax.set_xlabel(r'M($M_\odot h^{-1}$)(log)',size='x-large')
plt.savefig('HMF_nolcdm_%i.png'%num_bins)
#plt.show()
