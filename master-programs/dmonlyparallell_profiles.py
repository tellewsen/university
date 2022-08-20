import yt
yt.enable_parallelism()
import glob
from yt.utilities.parallel_tools.parallel_analysis_interface \
    import communication_system
from params import *
from datetime import datetime
import sys

#Mass limits for profiling
try:
    masslimitdown = np.log10(1e12) #Msun/h
    masslimitup   = np.log10(1e120) #Msun/h
except:
    print 'Correct usage: "python thisfile.py lowerlimit upperlimit" \
           Both limits are in units log10 Msun/h e.g 1e14 Msun/h -> 14'


#parameters for parallelizing
num_procs = 5
comm = communication_system.communicators[-1]

##################parameters for computation###########################
#Keep only dark matter sims
datapath = datapath[15:]
simtype = simtype[15:]

nbins = 64
numsims = len(datapath)

######################################################
def _dmoverdensity(field,data):
    global dmrhobar
    return data[('deposit','all_density')]/dmrhobar -1.

container= {}
avgdmdens = np.zeros((numsims,nbins)) #the global containers dm dens
avgdmodens = np.zeros((numsims,nbins)) #the global containers for overdesnity

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
    halolist = []
    for index,mass in enumerate(mvir):
        if masslimitdown <= mass <= masslimitup :
            halolist.append(index)

    #halolist = np.random.randint(len(halos),size=int(len(halos)*.2)) # random selection of 20% of the totalt halos in each sim
    numhalos = len(halolist)
    if numhalos < 1:
        print 'no halos to examine, exiting'
        sys.exit()

    #arrays for haloprofiles
    rrvir   = np.zeros((numhalos,nbins))
    dmdens  = np.zeros((numhalos,nbins))
    dmodens = np.zeros((numhalos,nbins))

    #Add overdenstiy field to dataset
    ad       = ds.all_data()
    dmrhobar = ad.mean(('deposit','all_density'),weight='cell_volume')

    ds.add_field(("deposit","dmoverdensity"), function=_dmoverdensity,take_log=False,force_override=True)#, units="auto"), dimensions="auto")

    #Calculate the profiles for each halo
    print datetime.now().strftime('%H:%M:%S'),'Starting halo loop'
    for j,i in enumerate(halolist): #run through the halolist
        cen = halos[i]
        sph = ds.sphere(cen, 10.*virial[i])

        #log radius
        profileDM   = yt.create_profile(sph,'radius',[('deposit','all_density')]    ,logs={'radius':True},weight_field='cell_volume',extrema={'radius':(1e-1*virial[i],10.*virial[i])})
        profileODDM = yt.create_profile(sph,'radius',[('deposit','dmoverdensity')]  ,logs={'radius':True},weight_field='cell_volume',extrema={'radius':(1e-1*virial[i],10.*virial[i])})


        #Pick out the radius scaled to virial radius for each halo
        rrvir[j]   = np.array(profileDM.x/virial[i])
        #Pick out field values at these radii
        dmdens[j]  = np.array(profileDM.field_data['deposit','all_density'].in_units('Msun/kpc**3'))
        dmodens[j] = np.array(profileODDM.field_data['deposit','dmoverdensity'])
    #Calculate the average values of temp and overdensity at each radii
    print datetime.now().strftime('%H:%M:%S'),'Starting averaging loop'
    for i in range(np.size(rrvir,1)):
        avgdmdens[sto.result_id][i]  = np.average(dmdens[:,i])
        avgdmodens[sto.result_id][i] = np.average(dmodens[:,i])

print datetime.now().strftime('%H:%M:%S'),'Combine all local arrays to global'

#Write results to file using cpu 0
comm.mpi_allreduce(avgdmdens)
comm.mpi_allreduce(avgdmodens)

if comm.rank == 0:
    print 'Saving to file using root process'
    np.savetxt('profsaveR.txt'   ,rrvir[0]   ,fmt='%.18e',delimiter=' ',header='Masslimitdown = %s, Masslimitup = %s, numhalos = %i'%(masslimitdown,masslimitup,numhalos))
    np.savetxt('profsaveDM.txt'  ,avgdmdens  ,fmt='%.18e',delimiter=' ',header='Masslimitdown = %s, Masslimitup = %s, numhalos = %i'%(masslimitdown,masslimitup,numhalos))
    np.savetxt('profsaveDMOD.txt',avgdmodens ,fmt='%.18e',delimiter=' ',header='Masslimitdown = %s, Masslimitup = %s, numhalos = %i'%(masslimitdown,masslimitup,numhalos))
