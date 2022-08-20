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
try:
    num_bins = int(sys.argv[1])#200 #Number of bins
except:
    print 'Correct usage: python thisfile.py num_bins'
    sys.exit(1)

binrange = (10,15) #bin range in log space

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

"""
alphaer = np.zeros(len(labels)*len(simtype))#GJENNOMSIKTIGHET FOR HISTOGRAMMER

#alphaer = [0,0,0,0,0]
count = 0
multip =0
for i in range(len(alphaer)):
    value = 1-0.2*multip
    alphaer[i] = value
    count += 1
    if count > 4:
        multip +=1
        count = 0

print alphaer
"""
lowhalo = 1e9
highhalo = 0


for k in range(len(labels)): #RUN OVER GRAVITY
    ax = fig.add_subplot(321+k) #subplot over GRAVITY
    for l in range(len(simtype)): #RUN OVER HYDRO
        mvir = np.loadtxt('Halos/%s.txt'%simtypes[5*l+k], usecols=8,unpack=True)
        mvir = np.log10(mvir)
        counts,bins = np.histogram(mvir,num_bins)#,range=binrange)
        #counts = np.array(counts,dtype=float)/len(mvir)

        #Generate a normal distribution with the same std and mean as mvir
        #[i for i,v in enumerate(a) if v > 4]
        halolist = []
        masslimitup = 14.
        masslimitdown = 10.
        for index,mass in enumerate(mvir):
            if masslimitdown <= mass <= masslimitup :
                halolist.append(index)

        #dist = mvir[halolist]
        #print dist
        #numchoice = 100000
        print simtypes[5*l+k],np.average(mvir)
        print simtypes[5*l+k],np.size(mvir)
        if np.size(mvir) < lowhalo:
            lowhalo = np.size(mvir)
        if np.size(mvir> highhalo):
            highhalo = np.size(mvir)
        #test1 = np.sum(mvir< np.average(mvir))
        #test2 = np.sum(mvir>= np.average(mvir))
        #print test1+test2
        print np.sum(mvir >= 14)
        #dist = np.random.choice(mvir,numchoice)
        #countsdist,bins = np.histogram(dist,num_bins,range=binrange)
        #countsdist = np.array(countsdist,dtype=float)/numchoice*len(halolist)

        #bins = bins[:-1] #Remove last edge
        #shift bins half a bin to the right
    	#bins += (bins[3]-bins[2])/2.

        #plt.plot(bins,countsdist,label='distrib')
        #plt.plot(bins,counts,label=simtype[l])#,'o')
        plt.hist(mvir,num_bins,label = simtype[l],histtype = 'step')#,alpha=alphaer[5*l+k])
        ax.set_title('%s'%labels[k])
        #ax.set_yscale('log')

        #ax.set_ylim([1e-6,1e-1])
        ax.set_ylim([0,3800])

        #ax.set_xscale('log')
        #ax.set_xlim([1e10,1e15])

    plt.legend()

print 'lowest halo number',lowhalo
print 'highest halo number',highhalo

figax.set_ylabel(r'$\#$ of halos')
figax.set_xlabel(r'Mass($M_\odot h^{-1}$)(log)')
figax.yaxis.set_label_coords(-0.075,0.5)
plt.savefig('HMFGvH_%i.png'%num_bins)
plt.show()
