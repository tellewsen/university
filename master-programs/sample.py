import yt
yt.enable_parallelism()
import numpy as np
from yt.analysis_modules.halo_analysis.api import HaloCatalog
import tempfile
import shutil
import os
from datetime import datetime

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from params import *

# Create temporary directory for storing files
#tmpdir = tempfile.mkdtemp()

#########DEFINE functions#########
#Define the overdensity field, not the meanOD parameter that is taken from the dataset
def _overdensity(field,data):
    global rhobar
    """Whatever you do, do not calculate rhobar inside this funciton.
       It will make a 2 second operation into a 30minute one."""
    return data[('ramses','Density')]/rhobar - 1.
####################################

#########ALLOCATE ARRAYS AND PARAMETERS NEEDED#########
#MAKE THE AVERAGE PROFILES FOR EACH OF THESE
nbins = 64
numsims = len(datapath)
numhalos = 100 #this should be a list generated from a fit to the halo mass distribution.

#Arrays for radii, odens, and temp for each halodir
#Will be overwritten each loop
rrvir = np.empty((numhalos,nbins))
odens = np.empty((numhalos,nbins))
temp  = np.empty((numhalos,nbins))

#Arrays for storing the average values of the halos in each sim for plotting
avgodens = np.empty((numsims,nbins))
avgtemp  = np.empty((numsims,nbins))

########################################################################


##########Run through all the sources and make profiels##################
for j in range(numsims):
    print datetime.now().strftime('%H:%M:%S'),'Loading data from %s_%s_11'%(simtype[j][2],simtype[j][0])
    #load sim
    ds = yt.load(datapath[j])

    #Define overdensity function (over density w/ respect to average density not critical density.
    ad = ds.all_data()
    #This is the first funciton to read the dataset, and mean is a costly operation

    #ad.set_field_parameter("meanOD", ad.mean(('ramses','Density')))
    #rhobar = ad.get_field_parameter("meanOD")
    rhobar = ad.mean(('ramses','Density'))
    #add overdensity field to use later, this is actually calculates
    ds.add_field(("gas","overdensity"), function=_overdensity,take_log=False,force_override=True)#, units="auto"), dimensions="auto")

    #Load halo coords, virial radius, and assign units to both
    halos    = ds.arr(np.loadtxt('Halos/Halos/%s.txt'%simtype[j][3],skiprows=1,unpack=False,usecols=(1,2,3)), 'Mpc/h')
    virial   = ds.arr(np.loadtxt('Halos/Halos/%s.txt'%simtype[j][3],skiprows=1,unpack=False,usecols=(7)),'Mpc/h')

    #Calculate the profiles for each halo
    for i in range(numhalos):
        print datetime.now().strftime('%H:%M:%S'),'Halo #',i
        cen = halos[i]
        sph = ds.sphere(cen, 5.*virial[i])
        #Calculate the overdensity and temperature profile for each halo
        profileD  = yt.create_profile(sph,'radius',[('gas','overdensity')],logs={'radius':False},weight_field='cell_volume',extrema={'radius':(0,5.*virial[i])})
        profileT  = yt.create_profile(sph,'radius',[('gas','temperature')],logs={'radius':False},weight_field='cell_mass',extrema={'radius':(0,5.*virial[i])})
        #Pick out the radius scaled to virial radius for each halo
        rrvir[i] = np.array(profileD.x/virial[i])
        #Pick out field values at these radii
        odens[i] = np.array(profileD.field_data['gas','overdensity'])
        temp[i]  = np.array(profileT.field_data['gas','temperature'])

    #Calculate the average values of temp and overdensity at each radii
    for i in range(np.size(rrvir,1)):
        avgtemp[j][i]  = np.average(temp[:,i])
        avgodens[j][i] = np.average(odens[:,i])
######################################################





##################PLOTTING THE PROFILES##################

#Format the plots
simtype = ['vanilla','cool star','cool star SN','DM']
labels    = [r'$\Lambda$CDM'    ,r'SymA'     , r'SymB'     , r'SymC'     ,r'SymD']
width = 10.24 #1024 pixels
height = 7.68 #768 pixels


#GvH Odens
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
    for l in range(len(simtype)-1): #RUN OVER HYDRO
        plt.plot(rrvir[0],avgodens[k+5*l],label=simtype[l])#,'o')
        #plt.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
        ax.set_title('%s'%labels[k])
        #ax.set_yscale('log')
        #ax.set_ylim([1e-6,1e-1])
        ax.set_xlim([0,5])
        #ax.set_xscale('log')
        #ax.set_xlim([1e10,1e15])
    plt.legend()
#figax.set_title(r'Halo mass function, $N_{bins}$=%i'%num_bins,y=1.05)
figax.set_ylabel(r'Overdensity $\frac{\rho-\bar{\rho}}{\bar{\rho}}$',labelpad=10)
figax.set_xlabel(r'$\rm{R/R_{vir}}$')
plt.savefig('AvgHaloProfOD_GvH.png')
#plt.show()

#GvH Temp
fig = plt.figure(1,frameon=False)
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
    for l in range(len(simtype)-1): #RUN OVER HYDRO
        plt.plot(rrvir[0],avgtemp[k+5*l],label=simtype[l])#,'o')
        #plt.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
        ax.set_title('%s'%labels[k])
        #ax.set_yscale('log')
        #ax.set_ylim([1e-6,1e-1])
        ax.set_xlim([0,5])
        #ax.set_xscale('log')
        #ax.set_xlim([1e10,1e15])
    plt.legend()
#figax.set_title(r'Halo mass function, $N_{bins}$=%i'%num_bins,y=1.05)
figax.set_ylabel(r'Temperature [K]',labelpad=10)
figax.set_xlabel(r'$\rm{R/R_{vir}}$')
plt.savefig('AvgHaloProfTemp_GvH.png')
#plt.show()

#HvG Overdensity
fig = plt.figure(2,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111)
# Turn off axis lines and ticks of the big subplot
figax.spines['top'].set_color('none')
figax.spines['bottom'].set_color('none')
figax.spines['left'].set_color('none')
figax.spines['right'].set_color('none')
figax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

for k in range(len(simtype)-1): #RUN OVER hydrotype
    ax = fig.add_subplot(221+k) #subplot over hydro
    for l in range(len(labels)): #RUN OVER gravity
        plt.plot(rrvir[0],avgodens[l+5*k],label=labels[l])#,'o')
        #plt.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
        ax.set_title('%s'%simtype[k])
        #ax.set_yscale('log')
        #ax.set_ylim([1e-6,1e-1])
        #ax.set_xscale('log')
        #ax.set_xlim([1e10,1e15])
        ax.set_xlim([0,5])
    plt.legend()
#figax.set_title(r'Halo mass function, $N_{bins}$=%i'%num_bins,y=1.05)
figax.set_ylabel(r'Overdensity $\frac{\rho-\bar{\rho}}{\bar{\rho}}$',labelpad=10)
figax.set_xlabel(r'$\rm{R/R_{vir}}$')
plt.savefig('AvgHaloProfOD_HvG.png')

#HvG temp
fig = plt.figure(3,frameon=False)
fig.set_size_inches(width,height)
figax = fig.add_subplot(111)
# Turn off axis lines and ticks of the big subplot
figax.spines['top'].set_color('none')
figax.spines['bottom'].set_color('none')
figax.spines['left'].set_color('none')
figax.spines['right'].set_color('none')
figax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

for k in range(len(simtype)-1): #RUN OVER hydrotype
    ax = fig.add_subplot(221+k) #subplot over hydro
    for l in range(len(labels)): #RUN OVER gravity
        plt.plot(rrvir[0],avgtemp[l+5*k],label=labels[l])#,'o')
        #plt.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
    ax.set_title('%s'%simtype[k])
    #ax.set_yscale('log')
    #ax.set_ylim([1e-6,1e-1])
    #ax.set_xscale('log')
    #ax.set_xlim([1e10,1e15])
    ax.set_xlim([0,5])
    plt.legend()
#figax.set_title(r'Halo mass function, $N_{bins}$=%i'%num_bins,y=1.05)
figax.set_ylabel(r'Temperature [K]',labelpad=10)
figax.set_xlabel(r'$\rm{R/R_{vir}}$')
plt.savefig('AvgHaloProfTemp_HvG.png')




"""
#ALL IN THE SAME PLOT /TERRIBLE/
for j in range(numsims):
    plt.figure(0)
    plt.plot(rrvir[0],avgodens[j],label=simtype[j][3])

    plt.figure(1)
    plt.semilogy(rrvir[0],avgtemp[j],label=simtype[j][3])

    plt.figure(2)
    plt.plot(rrvir[0],avgtemp[j],label=simtype[j][3])



plt.figure(0)
plt.xlabel(r'$\rm{R/R_{vir}}$')
plt.ylabel(r'Overdensity $\frac{\rho-\bar{\rho}}{\bar{\rho}}$')
plt.legend()
plt.savefig('AvgHaloProfOD.png')

plt.figure(1)
plt.xlabel(r'$\rm{R/R_{vir}}$')
plt.ylabel('Temperature [K]')
plt.legend()
plt.savefig('AvgHaloProfTemp.png')

plt.figure(2)
plt.xlabel(r'$\rm{R/R_{vir}}$')
plt.ylabel('Temperature [K]')
plt.legend()
plt.savefig('AvgHaloProfTempnonlog.png')
"""
######################################################
