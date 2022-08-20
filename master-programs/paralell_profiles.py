import yt
import numpy as np

yt.enable_parallelism()
#import glob
from yt.utilities.parallel_tools.parallel_analysis_interface \
    import communication_system
import params
from datetime import datetime
import sys

from yt.utilities.physical_constants import \
    boltzmann_constant_cgs, \
    mass_hydrogen_cgs#, mh


#Mass limits for profiling
masslimitdown = np.log10(10**0) #Msun/h
masslimitup   = np.log10(10**12) #Msun/h


#parameters for parallelizing
num_procs = 15
comm = communication_system.communicators[-1]

##################parameters for computation###########################
#Remove dark matter sims
datapath = params.datapath[:-5]
simtype = params.simtype[:-5]


nbins = 64
numsims = len(datapath)

######################################################

#Define fields
def _overdensity(field,data):
    global rhobar
    """Whatever you do, do not put ad and mean inside this funciton.
       It will make a 2 second operation into a 30minute one."""
    return data[('ramses','Density')]/rhobar - 1.

def _temperature(field,data):
    T  = data[('ramses','Pressure')]/data[('ramses','Density')]
    T *= .59*mass_hydrogen_cgs/boltzmann_constant_cgs# from diatomic gas
    return T

def _dmoverdensity(field,data):
    global dmrhobar
    return data[('deposit','all_density')]/dmrhobar -1.


#Defining particle filters for RAMSES
def Stars(pfilter,data):
    filter  = data[(pfilter.filtered_type,"particle_age")] != 0
    return filter

def DM(pfilter,data):
    filter = data[(pfilter.filtered_type,"particle_age")] == 0
    return filter

yt.add_particle_filter("stars",function=Stars,filtered_type='all',
                    requires=['particle_age'])

yt.add_particle_filter("dm",function=DM,filtered_type='all',
                    requires=['particle_age'])


#####################Data arrays############################
container= {}

#Averages
avgdens      = np.zeros((numsims,nbins)) #the global containers for density
avgtemp      = np.zeros((numsims,nbins)) #the global containers temperature
avgdmdens    = np.zeros((numsims,nbins)) #the global containers dm dens
avgstarsdens = np.zeros((numsims,nbins)) #the global containers dm dens
avgodens     = np.zeros((numsims,nbins)) #the global containers for overdesnity
avgdmodens   = np.zeros((numsims,nbins)) #the global containers for overdesnity

#Standard deviation
stddens      = np.zeros((numsims,nbins))
stdtemp      = np.zeros((numsims,nbins))
stddmdens    = np.zeros((numsims,nbins))
stdstarsdens = np.zeros((numsims,nbins))
stdodens     = np.zeros((numsims,nbins))
stddmodens   = np.zeros((numsims,nbins))
###############################################################


#Run over the simulations in parallel
for sto,path in yt.parallel_objects(datapath,num_procs,container):
    print datetime.now().strftime('%H:%M:%S'),'Loading Sim #',sto.result_id

    #Load data
    ds = yt.load(path)

    #pick the halos to be used in this simulation according to the mass limits
    simid    = simtype[sto.result_id][3] #sto.result_id works like enumerate
    halos    = ds.arr(np.loadtxt('Halos/%s.txt'%simid,skiprows=1,unpack=False,usecols=(1,2,3)),'Mpc/h')
    virial   = ds.arr(np.loadtxt('Halos/%s.txt'%simid,skiprows=1,unpack=False,usecols=(7))    ,'Mpc/h')
    mvir     = ds.arr(np.loadtxt('Halos/%s.txt'%simid,skiprows=1,unpack=False,usecols=(8))    ,'Msun/h')
    mvir     = np.log10(mvir)

    #Remove halos with less than x number of particles
    npart    = np.loadtxt('Halos/%s.txt'%simid,skiprows=1,unpack=False,usecols=0)
    minparts = 100
    for i in range(len(npart)):
        if npart[i] > minparts:
            minindex = i #find index of last halo with more than 100 particles
    halos = halos[:minindex] #remove the halos with too few
    virial = virial[:minindex]#remove the halos with too few
    mvir = mvir[:minindex]#remove the halos with too few

    #Fill halo list according to mass criterion
    halolist = []
    for index,mass in enumerate(mvir):
        if masslimitdown <= mass <= masslimitup :
            halolist.append(index)

    numhalos = len(halolist)
    if numhalos < 1:
        print 'no halos to examine, exiting'
        sys.exit()

    #arrays for haloprofiles
    rrvir     = np.zeros((numhalos,nbins))
    dens      = np.zeros((numhalos,nbins))
    temp      = np.zeros((numhalos,nbins))
    dmdens    = np.zeros((numhalos,nbins))
    starsdens = np.zeros((numhalos,nbins))
    odens     = np.zeros((numhalos,nbins))
    dmodens   = np.zeros((numhalos,nbins))

    #Add overdenstiy field to dataset
    #ad       = ds.all_data()
    #rhobar   = ad.mean(('ramses','Density'),weight='cell_volume')
    #dmrhobar = ad.mean(('deposit','all_density'),weight='cell_volume')
    #ds.add_field(("gas","overdensity"), function=_overdensity,take_log=False,force_override=True)#, units="auto"), dimensions="auto")
    #ds.add_field(("deposit","dmoverdensity"), function=_dmoverdensity,take_log=False,force_override=True)#, units="auto"), dimensions="auto")

    #Add temperature field with the right average mass
    ds.add_field(("gas","temperature"), function=_temperature,take_log=False,force_override=True, units="K")#, dimensions="auto")

    #Add particle fitlers for stars and dark matter
    if ( ('all','particle_age') in ds.field_list):
        ds.add_particle_filter("stars")
        ds.add_particle_filter("dm")

    #Calculate the profiles for each halo
    print datetime.now().strftime('%H:%M:%S'),'Starting halo loop'
    for j,i in enumerate(halolist): #run through the halolist
        cen = halos[i]
        sph = ds.sphere(cen, 10.*virial[i])

        #Profiles of all quantities and save field values
        profileD     = yt.create_profile(sph,'radius',[('gas','density')]            ,logs={'radius':True},weight_field='cell_volume',extrema={'radius':(1e-1*virial[i],10.*virial[i])})
        dens[j]      = np.array(profileD.field_data['gas','density'].in_units('Msun/kpc**3'))
        rrvir[j]     = np.array(profileD.x/virial[i]) #Pick out the radius scaled to virial radius for each halo

        profileT     = yt.create_profile(sph,'radius',[('gas','temperature')]        ,logs={'radius':True},weight_field='cell_mass'  ,extrema={'radius':(1e-1*virial[i],10.*virial[i])})
        temp[j]      = np.array(profileT.field_data['gas','temperature'].in_units('K'))

        #profileOD   = yt.create_profile(sph,'radius',[('gas','overdensity')]        ,logs={'radius':True},weight_field='cell_volume',extrema={'radius':(1e-1*virial[i],10.*virial[i])})
        #odens[j]     = np.array(profileOD.field_data['gas','overdensity'])

        if ( ('all','particle_age') in ds.field_list):
            profileDM    = yt.create_profile(sph,'radius',[('deposit','dm_density')] ,logs={'radius':True},weight_field='cell_volume',extrema={'radius':(1e-1*virial[i],10.*virial[i])})
            dmdens[j]    = np.array(profileDM.field_data['deposit','dm_density'].in_units('Msun/kpc**3'))

            profilestars = yt.create_profile(sph,'radius',[('deposit','stars_density')] ,logs={'radius':True},weight_field='cell_volume',extrema={'radius':(1e-1*virial[i],10.*virial[i])})
            starsdens[j] = np.array(profilestars.field_data['deposit','stars_density'].in_units('Msun/kpc**3'))
        else:
            profileDM    = yt.create_profile(sph,'radius',[('deposit','all_density')] ,logs={'radius':True},weight_field='cell_volume',extrema={'radius':(1e-1*virial[i],10.*virial[i])})
            dmdens[j]    = np.array(profileDM.field_data['deposit','all_density'].in_units('Msun/kpc**3'))

            #Note that stars are left as zeros

        #profileODDM = yt.create_profile(sph,'radius',[('deposit','dmoverdensity')]  ,logs={'radius':True},weight_field='cell_volume',extrema={'radius':(1e-1*virial[i],10.*virial[i])})
        #dmodens[j]   = np.array(profileODDM.field_data['deposit','dmoverdensity'])

    #Calculate the average values of temp and overdensity at each radii
    print datetime.now().strftime('%H:%M:%S'),'Starting averaging loop'
    for i in range(np.size(rrvir,1)):
        avgtemp[sto.result_id][i]      = np.mean(temp[:,i])
        avgdens[sto.result_id][i]      = np.mean(dens[:,i])
        avgdmdens[sto.result_id][i]    = np.mean(dmdens[:,i])
        avgstarsdens[sto.result_id][i] = np.mean(starsdens[:,i])
        #avgodens[sto.result_id][i]   = np.mean(odens[:,i])
        #avgdmodens[sto.result_id][i] = np.mean(dmodens[:,i])

        stdtemp[sto.result_id][i]      = np.std(temp[:,i])
        stddens[sto.result_id][i]      = np.std(dens[:,i])
        stddmdens[sto.result_id][i]    = np.std(dmdens[:,i])
        stdstarsdens[sto.result_id][i] = np.std(starsdens[:,i])
        #stdodens[sto.result_id][i]     = np.std(odens[:,i])
        #stddmodens[sto.result_id][i]   = np.std(dmodens[:,i])


print datetime.now().strftime('%H:%M:%S'),'Combine all local arrays to global'

#Write results to file using cpu 0
comm.mpi_allreduce(avgtemp)
comm.mpi_allreduce(avgdens)
comm.mpi_allreduce(avgdmdens)
comm.mpi_allreduce(avgstarsdens)
#comm.mpi_allreduce(avgdmodens)
#comm.mpi_allreduce(avgodens)

comm.mpi_allreduce(stdtemp)
comm.mpi_allreduce(stddens)
comm.mpi_allreduce(stddmdens)
comm.mpi_allreduce(stdstarsdens)
#comm.mpi_allreduce(stddmodens)
#comm.mpi_allreduce(stdodens)

if comm.rank == 0:
    print 'Saving to file using root process'
    np.savetxt('profsaveR.txt'     ,rrvir[0]     ,fmt='%.18e',delimiter=' ',header='Masslimitdown = %s, Masslimitup = %s, numhalos = %i'%(masslimitdown,masslimitup,numhalos))
    np.savetxt('profsaveT.txt'     ,avgtemp      ,fmt='%.18e',delimiter=' ',header='Masslimitdown = %s, Masslimitup = %s, numhalos = %i'%(masslimitdown,masslimitup,numhalos))
    np.savetxt('profsaveD.txt'     ,avgdens      ,fmt='%.18e',delimiter=' ',header='Masslimitdown = %s, Masslimitup = %s, numhalos = %i'%(masslimitdown,masslimitup,numhalos))
    np.savetxt('profsaveDM.txt'    ,avgdmdens    ,fmt='%.18e',delimiter=' ',header='Masslimitdown = %s, Masslimitup = %s, numhalos = %i'%(masslimitdown,masslimitup,numhalos))
    np.savetxt('profsavestars.txt' ,avgstarsdens ,fmt='%.18e',delimiter=' ',header='Masslimitdown = %s, Masslimitup = %s, numhalos = %i'%(masslimitdown,masslimitup,numhalos))
    #np.savetxt('profsaveOD.txt'  ,avgodens   ,fmt='%.18e',delimiter=' ',header='Masslimitdown = %s, Masslimitup = %s, numhalos = %i'%(masslimitdown,masslimitup,numhalos))
    #np.savetxt('profsaveDMOD.txt',avgdmodens ,fmt='%.18e',delimiter=' ',header='Masslimitdown = %s, Masslimitup = %s, numhalos = %i'%(masslimitdown,masslimitup,numhalos))

    np.savetxt('profsaveTstd.txt'     ,stdtemp      ,fmt='%.18e',delimiter=' ',header='Masslimitdown = %s, Masslimitup = %s, numhalos = %i'%(masslimitdown,masslimitup,numhalos))
    np.savetxt('profsaveDstd.txt'     ,stddens      ,fmt='%.18e',delimiter=' ',header='Masslimitdown = %s, Masslimitup = %s, numhalos = %i'%(masslimitdown,masslimitup,numhalos))
    np.savetxt('profsaveDMstd.txt'    ,stddmdens    ,fmt='%.18e',delimiter=' ',header='Masslimitdown = %s, Masslimitup = %s, numhalos = %i'%(masslimitdown,masslimitup,numhalos))
    np.savetxt('profsavestarsstd.txt' ,stdstarsdens ,fmt='%.18e',delimiter=' ',header='Masslimitdown = %s, Masslimitup = %s, numhalos = %i'%(masslimitdown,masslimitup,numhalos))
    #np.savetxt('profsaveODstd.txt'  ,stdodens   ,fmt='%.18e',delimiter=' ',header='Masslimitdown = %s, Masslimitup = %s, numhalos = %i'%(masslimitdown,masslimitup,numhalos))
    #np.savetxt('profsaveDMODstd.txt',stddmodens ,fmt='%.18e',delimiter=' ',header='Masslimitdown = %s, Masslimitup = %s, numhalos = %i'%(masslimitdown,masslimitup,numhalos))
