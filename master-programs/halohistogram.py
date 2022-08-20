import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif'}) # This is for Latex writing

#Halo mass histograms Same hydroeffect in same plot.

simtypes   = ['LCDM_vanilla'     ,'SymA_vanilla'     ,'SymB_vanilla'     ,'SymC_vanilla'     ,'SymD_vanilla',
              'LCDM_cool_star'   ,'SymA_cool_star'   ,'SymB_cool_star'   ,'SymC_cool_star'   ,'SymD_cool_star',
              'LCDM_cool_star_SN','SymA_cool_star_SN','SymB_cool_star_SN','SymC_cool_star_SN','SymD_cool_star_SN',
              'LCDM_DM'          ,'SymA_DM'          ,'SymB_DM'          ,'SymC_DM'          ,'SymD_DM']

simtype = ['vanilla','cool star','cool star SN','DM']
labels    = [r'$\Lambda$CDM'    ,r'SymA'     , r'SymB'     , r'SymC'     ,r'SymD']

Volume = 64.*64.*64. #( h(1+z)Mpc ) ^3
binrange = (10,15) #bin range in log space

try:
    num_bins = int(sys.argv[1])#200 #Number of bins
except:
    print 'Correct usage: python thisfile.py num_bins'
    sys.exit(1)


#window size
width = 10.24 #1024 pixels
height = 7.68 #768 pixels

#RUNS OVER GRAVITY SPIS OUT 5 PLOTS with 4 hydromodels in each
fig = plt.figure(0,frameon=False)

fig.set_size_inches(width,height)

figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
# Turn off axis lines and ticks of the big subplot
figax.spines['top'].set_color('none')
figax.spines['bottom'].set_color('none')
figax.spines['left'].set_color('none')
figax.spines['right'].set_color('none')
figax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

for k in range(len(labels)): #RUN OVER GRAVITY
    ax = fig.add_subplot(321+k) #subplot over GRAVITY
    for l in range(len(simtype)): #RUN OVER HYDRO

        mvir = np.loadtxt('Halos/%s.txt'%simtypes[5*l+k], usecols=8,unpack=True)
        mvir = np.log10(mvir)
        counts,bins = np.histogram(mvir,num_bins,range=binrange)
        bins = bins[:-1] #Remove last edge

        #shift bins half a bin to the right
    	bins += (bins[3]-bins[2])/2.

        # turn counts into counts per comoving Mpc
        counts = counts/Volume
        for i in range(len(bins)):
	        counts[i] = np.sum(counts[i:])

        plt.plot(bins,counts,label=simtype[l])#,'o')
        #plt.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
        ax.set_title('%s'%labels[k])
        ax.set_yscale('log')
        ax.set_ylim([1e-6,1e-1])
        ax.set_xlim([10,15])

        #ax.set_xscale('log')
        #ax.set_xlim([1e10,1e15])

    plt.legend()
figax.set_title(r'Halo mass function, $N_{bins}$=%i'%num_bins,y=1.05)
figax.set_ylabel(r'n(>M) $[h^3 Mpc^{-3}]$',labelpad=10)
figax.set_xlabel(r'Mass($M_\odot h^{-1}$)(log)')
plt.savefig('HMFGvH_%i.png'%num_bins)
#plt.show()


#Runs over hydromodel, spits out 4 plots with 5 gravitymodels in each.
fig = plt.figure(1,frameon=False)
width = 10.24
height = 7.68
fig.set_size_inches(width,height)
figax = fig.add_subplot(111)
# Turn off axis lines and ticks of the big subplot
figax.spines['top'].set_color('none')
figax.spines['bottom'].set_color('none')
figax.spines['left'].set_color('none')
figax.spines['right'].set_color('none')
figax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

for k in range(len(simtype)): #RUN OVER hydrotype
    ax = fig.add_subplot(221+k) #subplot over hydro
    for l in range(len(labels)): #RUN OVER gravity

        mvir = np.loadtxt('Halos/%s.txt'%simtypes[5*k+l], usecols=8,unpack=True)
        mvir = np.log10(mvir)
        counts,bins = np.histogram(mvir,num_bins,range=binrange)
        bins = bins[:-1] #Remove last edge

        #shift bins half a bin to the right
        bins += (bins[1]-bins[0])/2

        # turn counts into counts per comoving Mpc
        counts = counts/Volume # n -> n/(Mpc/h) = n [h Mpc^-1]^3 = n [h^3 Mpc^-3]
        for i in range(len(bins)):
	        counts[i] = np.sum(counts[i:])

        plt.plot(bins,counts,label=labels[l])#,'o')
        #plt.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
        ax.set_title('%s'%simtype[k])
        ax.set_yscale('log')
        ax.set_ylim([1e-6,1e-1])

        #ax.set_xscale('log')
        #ax.set_xlim([1e10,1e15])
        ax.set_xlim([10,15])

    plt.legend()

figax.set_title(r'Halo mass function, $N_{bins}$=%i'%num_bins,y=1.05)
figax.set_ylabel(r'n(>M) $ [h^3 Mpc^{-3}]$',labelpad=10)
figax.set_xlabel(r'Mass($M_\odot h^{-1}$)(log)')
#plt.show()
plt.savefig('HMFHvG_%i.png'%num_bins)


##############################RATIO PART##############################

#Make the standard to take ratios from
LCDMmvir = np.loadtxt('Halos/LCDM_vanilla.txt', usecols=8,unpack=True)
LCDMmvir = np.log10(LCDMmvir)
LCDMcounts,LCDMbins = np.histogram(LCDMmvir,num_bins,range=binrange)
LCDMbins = LCDMbins[:-1]
LCDMbins = LCDMbins + (LCDMbins[1]-LCDMbins[0])/2
LCDMcounts = LCDMcounts/Volume
for i in range(len(LCDMcounts)):
    LCDMcounts[i] = np.sum(LCDMcounts[i:])

fig = plt.figure(3,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111) #Big plot used for centerd x and y axis labels
# Turn off axis lines and ticks of the big subplot
figax.spines['top'].set_color('none')
figax.spines['bottom'].set_color('none')
figax.spines['left'].set_color('none')
figax.spines['right'].set_color('none')
figax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

#Plot ratios of each to LCDM
for k in range(len(labels)): #RUN OVER GRAVITY
    ax = fig.add_subplot(321+k) #subplot over GRAVITY
    for l in range(len(simtype)): #RUN OVER HYDRO

        mvir = np.loadtxt('Halos/%s.txt'%simtypes[5*l+k], usecols=8,unpack=True)
        mvir = np.log10(mvir)
        counts,bins = np.histogram(mvir,num_bins,range=binrange)
        bins = bins[:-1] #Remove last edge

        #shift bins half a bin to the right
        bins += (bins[1]-bins[0])/2. #define length of half a bin

        # turn counts into counts per comoving Mpc
        counts = counts/Volume
        for i in range(len(counts)):
	        counts[i] = np.sum(counts[i:])

        plt.plot(bins,LCDMcounts/counts-1.,label=simtype[l])#,'o')
        #plt.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
        ax.set_title('%s'%labels[k])
        #ax.set_yscale('symlog')
        #ax.set_ylim([-1,0.05])
        ax.set_xlim([10,15])

        #ax.set_xscale('log')
        #ax.set_xlim([1e10,1e15])

    plt.legend()
figax.set_title(r'Ratio HMF to $\Lambda$CDM vanilla HMF, $N_{bins}$=%i'%num_bins,y=1.05)
figax.set_ylabel(r'$\frac{n_{\Lambda CDM}(>M)}{n(>M)} -1$',labelpad=10)
figax.set_xlabel(r'Mass($M_\odot h^{-1}$)(log)')
plt.savefig('HMFGvH_ratio_%i.png'%num_bins)
#plt.show()



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


for k in range(len(simtype)): #RUN OVER HYDRO
    ax = fig.add_subplot(221+k) #subplot over HYDRO
    for l in range(len(labels)): #RUN OVER GRAVITY
        mvir = np.loadtxt('Halos/%s.txt'%simtypes[5*k+l], usecols=8,unpack=True)
        mvir = np.log10(mvir)
        counts,bins = np.histogram(mvir,num_bins,range=binrange)
        bins = bins[:-1] #Remove last edge

        #shift bins half a bin to the right
        bins += (bins[1]-bins[0])/2.

        # turn counts into counts per comoving Mpc
        counts = counts/Volume
        for i in range(len(counts)):
	        counts[i] = np.sum(counts[i:])

        plt.plot(bins,LCDMcounts/counts-1.,label=labels[l])#,'o')
        #plt.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
        ax.set_title('%s'%simtype[k])
        #ax.set_yscale('symlog')
        #ax.set_ylim([-1,0.05])
        ax.set_xlim([10,15])

        #ax.set_xscale('log')
        #ax.set_xlim([1e10,1e15])

    plt.legend(loc='lower left')
figax.set_title(r'Ratio HMF to $\Lambda$CDM vanilla HMF, $N_{bins}$=%i'%num_bins,y=1.05)
figax.set_ylabel(r'$\frac{n_{\Lambda CDM}(>M)}{n(>M)} -1 $',labelpad=10)
figax.set_xlabel(r'Mass($M_\odot h^{-1}$)(log)')
plt.savefig('HMFHvG_ratio_%i.png'%num_bins)
#plt.show()
