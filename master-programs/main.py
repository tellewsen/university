"""
This program reads data from a RAMSES cosmological simulation along with a halo catalog and produces various plots for it.
There is also implemented the possibility of reading a halomatching list matching halos from two simulations.

Volume rendering was not supported by yt for RAMSES simulations at the end of summer 2016 so that is not included.
Made by Thor Andreas Seiff Ellewsen summer 2016.

Updated 13/3/17
"""

#Import useful modules
from params import *
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import yt
from yt.units.yt_array import YTQuantity

#from yt.fields import cosmology_fields
from functions import *

#Define custom fields for yt
####################Define derived fields####################
def _overdensity(field,data):
    ad = data.ds.all_data()
    mean = ad.mean(('ramses','Density'))
    return data[('ramses','Density')]/mean -1
#Fields are not added here. They are instead added to the dataset when it is loaded. See main program after sim is loaded



"""
Here we load the data for each output in each simulation, and plot various things.
Note that this takes a VERY LONG time. I reccommend using a computer with a huge amount of memory.
In my experience, it takes about 50 minutes to make all the possible plots for on output folder.
Of that time, 20-30 minutes is the time it takes to load the output into memory before you can do anything with it.

For 2 simulations with 9 folders each, you end up with 15HOURS.
In other words, you usually do one or two folders from each sim when working on things, and when you're happy with it, you start the program when you leave work.

This example is for a 64cMpc**3 box.
"""

#Set number of sims and outputs for each
numsims  = len(simtype)
numsource = len(sourcelist)
#If no specific halos are chosen we use the 4 largest
if numsims > 1:
  if halolist == None:
    numhalos = 4 #
    halolist = np.empty((numsims,numhalos),dtype=int)
    for l,simname in enumerate(simtype):
        for i in range(numhalos):
            if match[l][i] != 'None':
                halolist[l][i] = match[l][i] #load halos for right sim
            else:
                halolist[l][i] = 9999 #Mark with this number if no matches were found
  else:
      numhalos = len(halolist)
else:
    numhalos = 4
    halolist = np.linspace(0,numhalos-1,numhalos,dtype=int)

#print 'halolist = ',halolist

#If no specific datasource is chosen we use the latest output for each sim.
if sourcelist == None:
    sourcelist = [len(ts)-1]

######################Prepare plot containers######################
#Plots of entire simulation
dens_whole_sim_plot            = np.empty((numsims,numsource,3),dtype=object)
slic_whole_plot                = np.empty((numsims,numsource,3),dtype=object)
part_mass_whole_plot           = np.empty((numsims,numsource,3),dtype=object)
temp_whole_sim_plot            = np.empty((numsims,numsource,3),dtype=object)
#Plots of halos
proj_plot_halo_plot            = np.empty((numsims,numsource,numhalos,3),dtype=object)
cube_proj_halo_plot            = np.empty((numsims,numsource,numhalos,3),dtype=object)
slic_halo_plot                 = np.empty((numsims,numsource,numhalos,3),dtype=object)
part_mass_halo_plot            = np.empty((numsims,numsource,numhalos,3),dtype=object)
temp_plot_halo_plot            = np.empty((numsims,numsource,numhalos,3),dtype=object)
pres_plot_halo_plot            = np.empty((numsims,numsource,numhalos,3),dtype=object)

halo_over                      = np.empty((numsims,numsource,numhalos,3),dtype=object)
halo_over_vel                  = np.empty((numsims,numsource,numhalos,3),dtype=object)
halo_over_conv                 = np.empty((numsims,numsource,numhalos,3),dtype=object)
halo_over_grid                 = np.empty((numsims,numsource,numhalos,3),dtype=object)
halo_over_stream               = np.empty((numsims,numsource,numhalos,3),dtype=object)

#Profiles of halos
prof_density_halo_plot_labels  = np.empty((numsims,numsource,numhalos),dtype=object)
prof_density_halo_plot_specs   = np.empty((numsims,numsource,numhalos),dtype=object)
prof_density_halo_plot         = np.empty((numsims,numsource,numhalos),dtype=object)

prof_temp_halo_plot_labels     = np.empty((numsims,numsource,numhalos),dtype=object)
prof_temp_halo_plot_specs      = np.empty((numsims,numsource,numhalos),dtype=object)
prof_temp_halo_plot            = np.empty((numsims,numsource,numhalos),dtype=object)

DM_dens_prof_plot_labels       = np.empty((numsims,numsource,numhalos),dtype=object)
DM_dens_prof_plot_specs        = np.empty((numsims,numsource,numhalos),dtype=object)
DM_dens_prof_plot              = np.empty((numsims,numsource,numhalos),dtype=object)

redshift_list                  = np.empty((numsims,numsource),dtype=object)
#Container for current extremas of output
gas_density_extrema            = np.zeros((numsims,numsource,2))
all_particle_mass_extrema      = np.zeros((numsims,numsource,2))
#Containers for the global extrema of all outputs
gas_density_extrema_glob       = np.empty((numsims,2),dtype=object)
all_particle_mass_extrema_glob = np.empty((numsims,2),dtype=object)
#Periodic boundary projection test container
PBCtest                        = np.empty((numsims,numsource),dtype=object)


########################LOAD DATA and run functions on it########################
for l in range(numsims):

    #Load data
    #print datetime.now().strftime('%H:%M:%S'),'Loading data from %s_%s'%(simtype[l][2],simtype[l][0])
    #ts = yt.load(simpath) #POSSIBLE REASON FOR THE MEMORY ISSUE

    #RUN VARIOUS FUNCTIONS TO SAVE PLOTS OF INTERESTING THINGS FOR ALL SIMULATIONTYPES,OUTPUTS,HALOS,DIMENSIONS.
    for i in range(numsource):
        #Define various things that are useful later
        print datetime.now().strftime('%H:%M:%S'),'Loading data from %s_%s_%s'%(simtype[l][2],simtype[l][0],sourcelist[i])
        ds = yt.load(datapath[l+i])

        #Add fields to dataset
        if ('ramses', 'Density') in ds.field_list: #Check for gas density in sim
            ds.add_field(("gas","overdensity"), function=_overdensity,take_log=True,force_override=True)#, units="auto"), dimensions="auto")

        #Test of DM particle mass
        #L = ds.length_unit.in_units('Mpc')

        #Gravconstant = YTQuantity(4.302*10**-9,'Mpc/Msun*(km/s)**2')
        #Gravconstant = YTQuantity(6.674*10**-11,'m**3/kg/s**2')
        #H = YTQuantity(100*ds.hubble_constant,'km/Mpc/s')
        #rho_c = 3*H**2/(8*np.pi*Gravconstant)
        #partmass = ds.omega_matter*rho_c*L**3/256**3


        #ad = ds.all_data()
        curr_src = str(ds)[5:11] #string used when saving plot

        #Add current redshift to list. This is the one for this sim and this output. Will use for profile plot of same halo in diff sim
        redshift_list[l,i] = 'z = %.2f'%ds.current_redshift

        """
        #############Save various extrema for simulations. Used to scale colorbars on plots.#############
        print datetime.now().strftime('%H:%M:%S'),'Calculating various extrema for current simulation'
        #Used for densities. Does not work on projections since that is not density anymore.
        if ('ramses', 'Density') in ds.field_list: #Check for gas density in sim
            gas_density_extrema[l][i]       = ad.quantities.extrema(('gas','density'))
            if(gas_density_extrema[l][i][0] < gas_density_extrema_glob[l][0] or gas_density_extrema_glob[l][0] == None):
                gas_density_extrema_glob[l][0] = gas_density_extrema[l][i][0]
            if(gas_density_extrema[l][i][1] > gas_density_extrema_glob[l][1] or gas_density_extrema_glob[l][1] == None):
                gas_density_extrema_glob[l][1] = gas_density_extrema[l][i][1]

        #Used for particle mass. Dark matter
        if ('io','particle_mass') in ds.field_list: #Check for DM particles in sim
            all_particle_mass_extrema[l][i] = ad.quantities.extrema(('all','particle_mass'))
            if(all_particle_mass_extrema[l][i][0] < all_particle_mass_extrema_glob[l][0] or all_particle_mass_extrema_glob[l][0] == None):
                all_particle_mass_extrema_glob[l][0] = all_particle_mass_extrema[l][i][0]
            if(all_particle_mass_extrema[l][i][1] > all_particle_mass_extrema_glob[l][1] or all_particle_mass_extrema_glob[l][1] == None):
                all_particle_mass_extrema_glob[l][1] = all_particle_mass_extrema[l][i][1]


        ##################PLOTS START HERE###################
        #Plots of whole simulation for each dimension
        if ('ramses', 'Density') in ds.field_list: #Check for gas density in sim
            for k,dim in enumerate(['x','y','z']):
                dens_whole_sim_plot[l,i,k]  = dens_whole_sim(ds,dim,'%s_%s'%(simtype[l][2],simtype[l][0]),curr_src)
                slic_whole_plot[l,i,k]      =     slic_whole(ds,dim,gas_density_extrema[l][i],'%s_%s'%(simtype[l][2],simtype[l][0]),curr_src)
        if ('io','particle_mass') in ds.field_list: #Check for DM particles in sim
            for k,dim in enumerate(['x','y','z']):
                part_mass_whole_plot[l,i,k] = part_mass_whole(ds,dim,'%s_%s'%(simtype[l][2],simtype[l][0]),curr_src)
        if ('gas','temperature') in ds.derived_field_list:
            for k,dim in enumerate(['x','y','z']):
                temp_whole_sim_plot[l,i,k] = temp_whole_sim(ds,dim,'%s_%s'%(simtype[l][2],simtype[l][0]),curr_src)
        """

        #################Make plots for each of the halos chosen#################
        for j in range(numhalos):
            curr_halo = halolist[l][j]
            if curr_halo == 9999: #Skip halos in sims where no match to original was found
                continue
            rad = haloradius[l][curr_halo]
            cen = halos[l][curr_halo] #Center of current halo
            wid = 2.*rad # diameter of current halo

            """
            #Slice plots for each simension
            if ("gas","baryon_overdensity") in ds.derived_field_list: #Check that overdensity is defined
                for k,dim in enumerate(['x','y','z']):
                  halo_over[l,i,j,k]       = halo_over_plot(ds,dim,cen,wid,'%s_%s'%(simtype[l][2],simtype[l][0]),curr_src,curr_halo)
                  halo_over_vel[l,i,j,k]   = halo_over_vel_plot(ds,dim,cen,wid,'%s_%s'%(simtype[l][2],simtype[l][0]),curr_src,curr_halo)
                  halo_over_conv[l,i,j,k]  = halo_over_conv_plot(ds,dim,cen,wid,'%s_%s'%(simtype[l][2],simtype[l][0]),curr_src,curr_halo)
                  halo_over_grid[l,i,j,k]  = halo_over_grid_plot(ds,dim,cen,wid,'%s_%s'%(simtype[l][2],simtype[l][0]),curr_src,curr_halo)
                  halo_over_stream[l,i,j,k]= halo_over_stream_plot(ds,dim,cen,wid,'%s_%s'%(simtype[l][2],simtype[l][0]),curr_src,curr_halo)

            #Plots of each dimension for each halos in each output
            if ('ramses', 'Density') in ds.field_list: #Check for gas density in sim
                left_edge  = cen - wid #A bit larger than the halo but we want that
                right_edge = cen + wid #A bit larger than the halo but we want that
                cube       = ds.region(cen, left_edge, right_edge, fields=None, ds=ds, field_parameters=None, data_source=None)
                for k,dim in enumerate(['x','y','z']):
                    proj_plot_halo_plot[l,i,j,k] = proj_plot_halo(ds,dim,cen,wid,'%s_%s'%(simtype[l][2],simtype[l][0]),curr_src,curr_halo)
                    cube_proj_halo_plot[l,i,j,k] = cube_proj_halo(ds,dim,cen,wid,'%s_%s'%(simtype[l][2],simtype[l][0]),curr_src,curr_halo,cube)
                    #slic_halo_plot[l,i,j,k]      =      slic_halo(ds,dim,cen,wid,'%s_%s'%(simtype[l][2],simtype[l][0]),curr_src,curr_halo)
                    temp_plot_halo_plot[l,i,j,k] = temp_plot_halo(ds,dim,cen,wid,'%s_%s'%(simtype[l][2],simtype[l][0]),curr_src,curr_halo)
                    pres_plot_halo_plot[l,i,j,k] = pres_plot_halo(ds,dim,cen,wid,'%s_%s'%(simtype[l][2],simtype[l][0]),curr_src,curr_halo)

            if ('io','particle_mass') in ds.field_list: #Check for DM particles in sim
                for k,dim in enumerate(['x','y','z']):
                    part_mass_halo_plot[l,i,j,k] = part_mass_halo(ds,dim,cen,wid,'%s_%s'%(simtype[l][2],simtype[l][0]),curr_src,curr_halo)
            """


            #Radial density profiles of halos
            sph = ds.sphere(cen, 2.*rad)
            """
            if ('ramses', 'Density') in ds.field_list: #Check for gas density in sim
                prof_density_halo_plot[l,i,j]        = prof_density_halo(ds,sph,'%s_%s'%(simtype[l][2],simtype[l][0]),curr_src,curr_halo)
                prof_density_halo_plot_specs[l,i,j]  = dict(linewidth=2, alpha=0.7)
                prof_density_halo_plot_labels[l,i,j] = 'Halo #%d'%curr_halo
            """
            """
            if ('gas','temperature') in ds.derived_field_list:
                prof_temp_halo_plot[l,i,j]        = prof_temp_halo(ds,sph,'%s_%s'%(simtype[l][2],simtype[l][0]),curr_src,curr_halo)
                prof_temp_halo_plot_specs[l,i,j]  = dict(linewidth=2, alpha=0.7)
                prof_temp_halo_plot_labels[l,i,j] = 'Halo #%d'%curr_halo
            """
            """
            if ('io','particle_mass') in ds.field_list: #Check for DM particles in sim
                DM_dens_prof_plot[l,i,j]        = DM_dens_prof(ds,sph,'%s_%s'%(simtype[l][2],simtype[l][0]),curr_src,curr_halo)
                DM_dens_prof_plot_specs[l,i,j]  = prof_density_halo_plot_specs[l,i,j]
                DM_dens_prof_plot_labels[l,i,j] = prof_density_halo_plot_labels[l,i,j]
            """
        """
        #Convert arrays with profiles to lists to make profile function happy.
        prof_density_halo_plotlist              = prof_density_halo_plot[l,i].tolist()
        prof_density_halo_plot_labelslist       = prof_density_halo_plot_labels[l,i].tolist()
        prof_density_halo_plot_specslist        = prof_density_halo_plot_specs[l,i].tolist()
        DM_dens_prof_plotlist                   = DM_dens_prof_plot[l,i].tolist()
        DM_dens_prof_plot_labelslist            = DM_dens_prof_plot_labels[l,i].tolist()
        DM_dens_prof_plot_specslist             = DM_dens_prof_plot_specs[l,i].tolist()

        #Remove all the None values to make profile funcion happy
        if None in prof_density_halo_plotlist: prof_density_halo_plotlist.remove(None)
        if None in prof_density_halo_plot_labelslist: prof_density_halo_plot_labelslist.remove(None)
        if None in prof_density_halo_plot_specslist: prof_density_halo_plot_specslist.remove(None)
        if None in DM_dens_prof_plotlist: DM_dens_prof_plotlist.remove(None)
        if None in DM_dens_prof_plot_labelslist: DM_dens_prof_plot_labelslist.remove(None)
        if None in DM_dens_prof_plot_specslist: DM_dens_prof_plot_specslist.remove(None)

        #Plot profileplot of the halos in this output in same plot. Not very useful.
        #Skips halos in simulations that have no match to the original ones
        gas_halos_plot = yt.ProfilePlot.from_profiles(profiles   = prof_density_halo_plotlist,
                                                      labels     = prof_density_halo_plot_labelslist,
                                                      plot_specs = prof_density_halo_plot_specslist)
        gas_halos_plot.set_unit('radius', 'kpc')
        gas_halos_plot.set_unit('density', 'Msun/kpc**3')
        gas_halos_plot.set_log('radius',False)
        gas_halos_plot.save('yt_%s_%s_allprof_gas.png'%('%s_%s'%(simtype[l][2],simtype[l][0]),curr_src))

        DM_halos_plot = yt.ProfilePlot.from_profiles(profiles   = DM_dens_prof_plotlist,
                                                     labels     = DM_dens_prof_plot_labelslist,
                                                     plot_specs = DM_dens_prof_plot_specslist)
        DM_halos_plot.set_unit('radius', 'kpc')
        DM_halos_plot.set_unit('io_density', 'Msun/kpc**3')
        DM_halos_plot.set_log('radius',False)
        DM_halos_plot.save('yt_%s_%s_allprof_DM.png'%('%s_%s'%(simtype[l][2],simtype[l][0]),curr_src))
        """



        #PBC testing
        #This plots the bottom and top middle part of the gas density projected through the whole y axis.
        #look at this together with the corresponding projection of the whole sim.
        #testcen = np.array([1.0,0.5,0.5])
        #testwid = 10./64
        #PBCtest[l][i]  = yt.ProjectionPlot(ds, "y", 'density', center=testcen, width=testwid)
        #PBCtest[l][i].set_axes_unit('Mpccm/h')
        #PBCtest[l][i].set_unit('density', 'Msun/kpc**2')
        #PBCtest[l][i].save('PBCtest_%s_%s.png'%('%s_%s'%(simtype[l][2],simtype[l][0]),curr_src))


        del ds






"""
#Combine various plots to make plots to compare halos at different redshifts for each model.
figs = np.zeros(numhalos,dtype=object)
haloevo = np.zeros(numsource,dtype=object)

for l in range(numsims):
    for j in range(numhalos):
        figs[j] = plt.figure()
        #Make multipanel evolution plots of halos for one sim
        grlena = int(np.sqrt(numsource))
        grlenb = int(np.ceil(numsource/float(grlena)))
        grid = AxesGrid(figs[j], (0.075,0.075,0.85,0.85),
                nrows_ncols = (grlenb, grlena),
                axes_pad = 0.05,
                add_all = True,
                label_mode = "1",
                share_all = True,
                aspect = True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="3%",
                cbar_pad="0%")

        for i,k in enumerate(sourcelist):
            if cube_proj_halo_plot[l,i,j,0] == None:
                continue
            else:
                #for k in range(3): Currently uses only x axis, can be generalized to z if needed.
                cube_proj_halo_plot[l,i,j,0].set_zlim('density',1e5,1e8) #Needs tweaking
                cube_proj_halo_plot[l,i,j,0].annotate_timestamp(corner='upper_left', time=False,redshift=True, draw_inset_box=False)
                haloevo[i]        = cube_proj_halo_plot[l,i,j,0].plots['density']
                haloevo[i].figure = figs[j]
                haloevo[i].axes   = grid[i].axes
                haloevo[i].cax    = grid.cbar_axes[i]
                cube_proj_halo_plot[l,i,j,0]._setup_plots()
        figs[j].savefig('multiplot_2x2_time_series_%s_%d.png'%('%s_%s'%(simtype[l][2],simtype[l][0]),halolist[l][j]))
"""

"""
#Plot projecitons of same halo in different sim at different timesteps.
figs = np.zeros(numhalos,dtype=object)
haloevo = np.zeros((numsims,numsource),dtype=object)
for j in range(numhalos):
    figs[j] = plt.figure()
    #Make multipanel evolution plots of halos for one sim
    grlena = numsource
    grlenb = numsims
    grid = AxesGrid(figs[j], (0.075,0.075,0.85,0.85),
                nrows_ncols = (grlena, grlenb),
                axes_pad = 0.05,
                add_all = True,
                label_mode = "1",
                share_all = True,
                aspect = True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="3%",
                cbar_pad="0%")
    for l in range(numsims):
        for i,k in enumerate(sourcelist):
            #for k in range(3): Currently uses only x axis, can be generalized to z if needed.
            cube_proj_halo_plot[l,i,j,0].set_zlim('density',1e5,1e9)
            cube_proj_halo_plot[l,i,j,0].annotate_timestamp(corner='upper_left', time=False,redshift=True, draw_inset_box=False)
            haloevo[l][i]        = cube_proj_halo_plot[l,i,j,0].plots['density']
            haloevo[l][i].figure = figs[j]
            haloevo[l][i].axes   = grid[l+2*i].axes
            haloevo[l][i].cax    = grid.cbar_axes[l+2*i]
            cube_proj_halo_plot[l,i,j,0]._setup_plots()
    figs[j].savefig('multiplot_2x2_time_series_%d_compare.png'%(halolist[l][j]))
"""

"""
#Plot profiles of same halos in different sims at one redshift for each figure
compareprofiles_g  = np.zeros((numhalos,numsims),dtype=object)
comparelabels_g    = np.zeros((numhalos,numsims),dtype=object)
comparespecs_g     = np.zeros((numhalos,numsims),dtype=object)
compareprofiles_dm = np.zeros((numhalos,numsims),dtype=object)
comparelabels_dm   = np.zeros((numhalos,numsims),dtype=object)
comparespecs_dm    = np.zeros((numhalos,numsims),dtype=object)
compareprofiles_temp = np.zeros((numhalos,numsims),dtype=object)
comparelabels_temp   = np.zeros((numhalos,numsims),dtype=object)
comparespecs_temp    = np.zeros((numhalos,numsims),dtype=object)


for i,k in enumerate(sourcelist): #output(redshift)
    currsource = sourcelist[i]
    for j in range(numhalos):#halo
        for l in range(numsims):#simulation type
            #pick out profiles for the current output(redshift)
            compareprofiles_g[j][l]  = prof_density_halo_plot[l,i,j]
            comparelabels_g[j][l]    = '%s %s'%(simtype[l][2],simtype[l][0])
            comparespecs_g[j][l]     = prof_density_halo_plot_specs[l,i,j]

            compareprofiles_dm[j][l] = DM_dens_prof_plot[l,i,j]
            comparelabels_dm[j][l]   = '%s %s'%(simtype[l][2],simtype[l][0])
            comparespecs_dm[j][l]    = DM_dens_prof_plot_specs[l,i,j]

            compareprofiles_temp[j][l] = prof_temp_halo_plot[l,i,j]
            comparelabels_temp[j][l]   = '%s %s'%(simtype[l][2],simtype[l][0])
            comparespecs_temp[j][l]    = prof_temp_halo_plot_specs[l,i,j]

        #pick out the profiles for the current halo (Now we have the profiles for all sims of the same halo)
        compareprofiles_g_list  = compareprofiles_g[j].tolist()
        comparelabels_g_list    = comparelabels_g[j].tolist()
        comparespecs_g_list     = comparespecs_g[j].tolist()

        compareprofiles_dm_list = compareprofiles_dm[j].tolist()
        comparelabels_dm_list   = comparelabels_dm[j].tolist()
        comparespecs_dm_list    = comparespecs_dm[j].tolist()

        compareprofiles_temp_list  = compareprofiles_temp[j].tolist()
        comparelabels_temp_list    = comparelabels_temp[j].tolist()
        comparespecs_temp_list     = comparespecs_temp[j].tolist()

        #Remove any empty objects in case we didn't find any matches for some sims
        if None in compareprofiles_g_list: compareprofiles_g_list.remove(None)
        if None in comparelabels_g_list: comparelabels_g_list.remove(None)
        if None in comparespecs_g_list: comparespecs_g_list.remove(None)

        if None in compareprofiles_dm_list: compareprofiles_dm_list.remove(None)
        if None in comparelabels_dm_list: comparelabels_dm_list.remove(None)
        if None in comparespecs_dm_list: comparespecs_dm_list.remove(None)

        if None in compareprofiles_temp_list: compareprofiles_temp_list.remove(None)
        if None in comparelabels_temp_list: comparelabels_temp_list.remove(None)
        if None in comparespecs_temp_list: comparespecs_temp_list.remove(None)



        currhalo = halolist[0][j] #NOTE THAT WE USE THE NUMBERING OF HALOS FROM SIMULATION 0!

        #gas_prof_dens_all_plot = yt.ProfilePlot.from_profiles(profiles=compareprofiles_g_list,
        #                                                      labels=comparelabels_g_list,
        #                                                      plot_specs=comparespecs_g_list)
        #gas_prof_dens_all_plot.set_unit('radius', 'kpc')
        #gas_prof_dens_all_plot.set_unit('density', 'Msun/kpc**3')
        #gas_prof_dens_all_plot.set_log('radius',False)
        #gas_prof_dens_all_plot.save('yt_%s_%s_haloprofiles_gas.png'%(currsource,currhalo))


        DM_prof_dens_all_plot  = yt.ProfilePlot.from_profiles(profiles=compareprofiles_dm_list,
                                                              labels=comparelabels_dm_list,
                                                              plot_specs=comparespecs_dm_list)
        DM_prof_dens_all_plot.set_unit('radius', 'kpc')
        DM_prof_dens_all_plot.set_unit('io_density', 'Msun/kpc**3')
        DM_prof_dens_all_plot.set_log('radius',False)
        DM_prof_dens_all_plot.save('yt_%s_%s_haloprofiles_DM.png'%(currsource,currhalo))

        gas_prof_temp_all_plot = yt.ProfilePlot.from_profiles(profiles=compareprofiles_temp_list,
                                                              labels=comparelabels_temp_list,
                                                              plot_specs=comparespecs_temp_list)
        gas_prof_temp_all_plot.set_unit('radius', 'kpc')
        gas_prof_temp_all_plot.set_unit('temperature', 'K')
        gas_prof_temp_all_plot.set_log('radius',False)
        gas_prof_temp_all_plot.save('yt_%s_%s_haloprofiles_temp.png'%(currsource,currhalo))
"""

"""
#Compare projections of whole sims in same figure. Maybe also every timestep?
figs = np.zeros(3,dtype=object)
simevo = np.zeros((numsims,numsource),dtype=object)
for k,dim in enumerate(['x','y','z']):
    figs[k] = plt.figure()
    if numsource*numsims < 4:
        grlena = int(np.ceil(numsource/2.))
        grlenb = 2*numsims
    else:
        grlena = 2#numsource
        grlenb = 3#numsims

    grid = AxesGrid(figs[k], (0.075,0.075,0.85,0.85),
                nrows_ncols = (grlena, grlenb),
                axes_pad = 0.05,
                add_all = True,
                label_mode = "1",
                share_all = True,
                aspect = True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="3%",
                cbar_pad="0%")
    for l,name in enumerate(simtype):
        for i,source in enumerate(sourcelist):
            if dens_whole_sim_plot[l,i,k] is not None:
                dens_whole_sim_plot[l,i,k].annotate_clear()
                dens_whole_sim_plot[l,i,k].set_zlim('density',1e5,1e9)
                #dens_whole_sim_plot[l,i,k].annotate_timestamp(corner='upper_left', time=False,redshift=True, draw_inset_box=False)
                simevo[l][i]        = dens_whole_sim_plot[l,i,k].plots['density']
                simevo[l][i].figure = figs[k]
                simevo[l][i].axes   = grid[l+2*i].axes
                simevo[l][i].cax    = grid.cbar_axes[l+2*i]
                dens_whole_sim_plot[l,i,k]._setup_plots()
    figs[k].savefig('multiplot_evolution_projected_density_z0_%s.png'%(dim))
"""
"""
#Compare temp projections of whole sims in same figure. Maybe also every timestep?
figs = np.zeros(3,dtype=object)
simevo = np.zeros((numsims,numsource),dtype=object)
for k,dim in enumerate(['x','y','z']):
    figs[k] = plt.figure()
    if numsource*numsims < 4:
        grlena = int(np.ceil(numsource/2.))
        grlenb = 2*numsims
    else:
        grlena = 2#numsource
        grlenb = 3#numsims

    grid = AxesGrid(figs[k], (0.075,0.075,0.85,0.85),
                nrows_ncols = (grlena, grlenb),
                axes_pad = 0.05,
                add_all = True,
                label_mode = "1",
                share_all = True,
                aspect = True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="3%",
                cbar_pad="0%")
    for l,name in enumerate(simtype):
        for i,source in enumerate(sourcelist):
            if temp_whole_sim_plot[l,i,k] is not None:
                temp_whole_sim_plot[l,i,k].annotate_clear()
                temp_whole_sim_plot[l,i,k].set_zlim('temperature',1e5,1e9)
                #dens_whole_sim_plot[l,i,k].annotate_timestamp(corner='upper_left', time=False,redshift=True, draw_inset_box=False)
                simevo[l][i]        = temp_whole_sim_plot[l,i,k].plots['temperature']
                simevo[l][i].figure = figs[k]
                simevo[l][i].axes   = grid[l+2*i].axes
                simevo[l][i].cax    = grid.cbar_axes[l+2*i]
                temp_whole_sim_plot[l,i,k]._setup_plots()
    figs[k].savefig('multiplot_evolution_projected_temp_z0_%s.png'%(dim))
"""
"""
#Compare DM projections of whole sims in same figure. Maybe also every timestep?
figs = np.zeros(3,dtype=object)
simevo = np.zeros((numsims,numsource),dtype=object)
for k,dim in enumerate(['x','y','z']):
    figs[k] = plt.figure()
    if numsource*numsims < 4:
        grlena = int(np.ceil(numsource/2.))
        grlenb = 2*numsims
    else:
        grlena = 2#numsource
        grlenb = 3#numsims

    grid = AxesGrid(figs[k], (0.075,0.075,0.85,0.85),
                nrows_ncols = (grlena, grlenb),
                axes_pad = 0.05,
                add_all = True,
                label_mode = "1",
                share_all = True,
                aspect = True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="3%",
                cbar_pad="0%")
    for l,name in enumerate(simtype):
        for i,source in enumerate(sourcelist):
            if part_mass_whole_plot[l,i,k] is not None:
                part_mass_whole_plot[l,i,k].annotate_clear()
                part_mass_whole_plot[l,i,k].set_zlim('particle_mass',1e11,1e13)
                #dens_whole_sim_plot[l,i,k].annotate_timestamp(corner='upper_left', time=False,redshift=True, draw_inset_box=False)
                simevo[l][i]        = part_mass_whole_plot[l,i,k].plots['particle_mass']
                simevo[l][i].figure = figs[k]
                simevo[l][i].axes   = grid[l+2*i].axes
                simevo[l][i].cax    = grid.cbar_axes[l+2*i]
                part_mass_whole_plot[l,i,k]._setup_plots()
    figs[k].savefig('multiplot_evolution_projected_DM_z0_%s.png'%(dim))
"""

#############EVERYTHING BELOW THIS POINT IS EXPERIMENTAL AND DOES NOT WORK#############
"""
    if prof_io_density_halo == True:
      field = "io_density"
      print datetime.now().strftime('%H:%M:%S'),'Making DM Density Profile Plots of halo #%d along each axis'%(curr_halo)
      sph = ds.sphere(cen, (5.*wid, 'Mpccm/h'))
      plot = yt.ProfilePlot(sph, "radius", field)#, weight_field="cell_mass")
      plot.set_unit('radius', 'kpc')
      plot.set_unit(field, 'Msun/kpc**3')
      plot.save('yt_%s_%d_prof_%s.png'%(simtype,curr_halo,field))

    if prof_all_density_halo == True:
      field = "all_density"
      print datetime.now().strftime('%H:%M:%S'),'Making DM Density Profile Plots of halo #%d along each axis'%(curr_halo)
      sph = ds.sphere(cen, (5.*wid, 'Mpccm/h'))
      plot = yt.ProfilePlot(sph, "radius", field)#, weight_field="cell_mass")
      plot.set_unit('radius', 'kpc')
      plot.set_unit('all_density', 'Msun/kpc**3')
      plot.save('yt_%s_%d_prof_%s.png'%(simtype,curr_halo,field))

"""

######Volume rendering ######
"""
#Volume rendering is not support by yt for RAMSES simulations. Maybe there will be som change to that in the future.

# Choose a field
field = 'density'
# Do you want the log of the field?
use_log = True
# Find the bounds in log space for your field
dd = ds.all_data()
mi, ma = dd.quantities.extrema(field)

if use_log:
    mi,ma = np.log10(mi), np.log10(ma)

# Instantiate the ColorTransferfunction.
tf = yt.ColorTransferFunction((mi, ma))

# Set up the camera parameters: center, looking direction, width, resolution
c = (ds.domain_right_edge + ds.domain_left_edge)/2.0
L = np.array([1.0, 1.0, 1.0])
W = ds.quan(1., 'Mpccm/h')
N = 256 #Increase for better resolution. Final picture has N*N pixels

# Create a camera object
cam = ds.camera(c, L, W, N, tf, fields = [field], log_fields = [use_log])

# Now let's add some isocontours, and take a snapshot, saving the image
# to a file.
tf.add_layers(10, 0.01, colormap = 'RdBu_r') #


im = cam.snapshot('test_rendering.png')

# To add the domain box to the image:
#nim = cam.draw_domain(im)
#nim.write_png('test_rendering_with_domain.png')

# To add the grid outlines to the image:
#nim = cam.draw_grids(im)
#nim.write_png('test_rendering_with_grids.png')
"""
