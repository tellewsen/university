import yt
yt.enable_parallelism()
import glob
from yt.utilities.parallel_tools.parallel_analysis_interface \
    import communication_system
from params import *
from datetime import datetime
import sys

from yt.utilities.physical_constants import \
    boltzmann_constant_cgs, \
    mass_hydrogen_cgs, \
    mh

#parameters for parallelizing
num_procs = 15
comm = communication_system.communicators[-1]

##################parameters for computation###########################
#Remove dark matter sims
datapath = datapath[:-5]
simtype = simtype[:-5]
numsims = len(datapath)

######################################################

container= {}
avggdens  = np.zeros(numsims) #the global containers for density
avgdmdens = np.zeros(numsims) #the global containers dm dens

#Run over the simulations in parallel
    for sto,path in yt.parallel_objects(datapath,num_procs,container):
        print datetime.now().strftime('%H:%M:%S'),'Loading Sim #',sto.result_id

    #Load data
    ds = yt.load(path)

    ad = ds.all_data()
    avggdens[sto.result_id] = ad.mean(('ramses','Density'),weight='cell_volume').in_units('Msun/kpc**3')
    avgdmdens[sto.result_id] = ad.mean(('deposit','all_density'),weight='cell_volume').in_units('Msun/kpc**3')

print datetime.now().strftime('%H:%M:%S'),'Combine all local arrays to global'

#Write results to file using cpu 0
comm.mpi_allreduce(avggdens)
comm.mpi_allreduce(avgdmdens)

if comm.rank == 0:
    print 'Saving to file using root process'
    np.savetxt('avggasdens.txt',avggdens    ,fmt='%.18e',delimiter=' ',header='Units = Msun/kpc**3')
    np.savetxt('avgdmdens.txt',avgdmdens ,fmt='%.18e',delimiter=' ',header='Units = Msun/kpc**3')
