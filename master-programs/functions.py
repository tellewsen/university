"""
This is the file containing all the various functions used by main.py
They are all specialized versions of yt plotting routines

Made by Thor Andreas Seiff Ellewsen summer 2016.

Updated 13/12/16
"""

import yt
from datetime import datetime

#################FUNCTIONS START#################
def halo_over_plot(ds,dim,cen,wid,simtype,curr_src,curr_halo):
    print datetime.now().strftime('%H:%M:%S'),'Making overdensity slice plot of halo #%d along %s axis'%(curr_halo,dim)
    field = 'baryon_overdensity'
    s = yt.SlicePlot(ds, dim, field, center=cen, width=wid)
    s.set_axes_unit('Mpccm/h')
    #s.set_zlim('baryon_overdensity',-2,210)    
    s.save('%s_%s_%d_slice_%s_%s.png'%(simtype,curr_src,curr_halo,field,dim))

def halo_over_vel_plot(ds,dim,cen,wid,simtype,curr_src,curr_halo):
    print datetime.now().strftime('%H:%M:%S'),'Making overdensity slice plot of halo #%d along %s axis with velocity annotated'%(curr_halo,dim)
    field = 'baryon_overdensity'
    s = yt.SlicePlot(ds, dim, field, center=cen, width=wid)
    if dim == 'x':
        s.annotate_quiver('velocity_y', 'velocity_z', 16)
    elif dim == 'y':
        s.annotate_quiver('velocity_x', 'velocity_z', 16)
    elif dim == 'z':
        s.annotate_quiver('velocity_x', 'velocity_y', 16)
    else:
        raise ValueError("Function halo_over_vel_plot received an unknown dimension '%s'. Use 'x','y','z'. This is case sensitive"%(dim))
    s.set_axes_unit('Mpccm/h')
    #s.set_zlim('baryon_overdensity',-2,210)    
    s.save('%s_%s_%s_slice_%s_%s_vel.png'%(simtype,curr_src,curr_halo,field,dim))

def halo_over_conv_plot(ds,dim,cen,wid,simtype,curr_src,curr_halo):
    print datetime.now().strftime('%H:%M:%S'),'Making overdensity slice plot of halo #%d along %s axis with convolution line integral annotated'%(curr_halo,dim)
    field = 'baryon_overdensity'
    s = yt.SlicePlot(ds, dim, field, center=cen, width=wid)
    if dim == 'x':
        s.annotate_line_integral_convolution('velocity_y', 'velocity_z', lim=(0.5,0.65))
    elif dim == 'y':
        s.annotate_line_integral_convolution('velocity_x', 'velocity_z', lim=(0.5,0.65))
    elif dim == 'z':
        s.annotate_line_integral_convolution('velocity_x', 'velocity_y', lim=(0.5,0.65))
    else:
        raise ValueError("Function halo_over_conv_plot received an unknown dimension '%s'. Use 'x','y','z'. This is case sensitive"%(dim))
    #s.set_zlim('baryon_overdensity',-2,210)
    s.save('%s_%s_%d_slice_%s_%s_conv.png'%(simtype,curr_src,curr_halo,field,dim))

def halo_over_grid_plot(ds,dim,cen,wid,simtype,curr_src,curr_halo):
    print datetime.now().strftime('%H:%M:%S'),'Making overdensity slice plot of halo #%d along %s axis with grid annotated'%(curr_halo,dim)
    field = 'baryon_overdensity'
    s = yt.SlicePlot(ds, dim, field, center=cen, width=wid)
    s.annotate_grids()
    #s.set_zlim('baryon_overdensity',-2,210)
    s.save('%s_%s_%d_slice_%s_%s_grid.png'%(simtype,curr_src,curr_halo,field,dim))

def halo_over_stream_plot(ds,dim,cen,wid,simtype,curr_src,curr_halo):
    field = 'baryon_overdensity'
    print datetime.now().strftime('%H:%M:%S'),'Making overdensity slice plot of halo #%d along %s axis with streamlines annotated'%(curr_halo,dim)
    s = yt.SlicePlot(ds, dim, field, center=cen, width=wid)
    if dim == 'x':
        s.annotate_streamlines('velocity_y', 'velocity_z')
    elif dim == 'y':
        s.annotate_streamlines('velocity_x', 'velocity_z')
    elif dim == 'z':
        s.annotate_streamlines('velocity_x', 'velocity_y')
    else:
        raise ValueError("Function halo_over_stram_plot received an unknown dimension '%s'. Use 'x','y','z'. This is case sensitive"%(dim))    
    #s.set_zlim('baryon_overdensity',-2,210)
    s.save('%s_%s_%d_slice_%s_%s_stream.png'%(simtype,curr_src,curr_halo,field,dim))




###WHOLE SIMS PLOTS###

#Density Slice plot whole sim through center along each axis. Not particularly useful.
def slic_whole(ds,dim,gas_density_extrema,simtype,curr_src):
    """
    ds - dataset
    dim  - 'x','y','z', or array [1,0.5,0.2] for off axis slice
    gas_density_extrema - array [min,max]
    simtype - string 'LCDM' /'Sym D' / etc
    curr_src - current output folder as string, e.g '00009'
    Function assumes that the gas_density_extrema are in cgs units (g/cm**3)
    """
    print datetime.now().strftime('%H:%M:%S'),'Making Slice plots of whole sim along %s axis'%dim

    #Apply units to the extremas and convert to solar masses per kilparsec cubed
    gas_density_extrema = gas_density_extrema*yt.units.g/yt.units.cm**3
    gas_density_extrema = gas_density_extrema.in_units('Msun/kpc**3')

    field = 'density'

    output = yt.SlicePlot(ds, dim, field)
    #save plot object to file
    #output.frb.save_as_dataset('yt_%s_%s_%s_slice_%s'%(simtype,curr_src,field,dim))

    output.annotate_timestamp(corner='upper_left', time=False,redshift=True, draw_inset_box=False)
    output.annotate_scale(corner='upper_right')
    output.set_axes_unit('Mpccm/h')
    output.set_zlim(('gas','density'), gas_density_extrema[0],gas_density_extrema[1])
    output.set_unit(field, 'Msun/kpc**3')
    output.save('%s_%s_%s_slice_%s.png'%(simtype,curr_src,field,dim))
    return output

#Particle projection with mass weights plot
def part_mass_whole(ds,dim,simtype,curr_src):
    print datetime.now().strftime('%H:%M:%S'),'Making Particle Projection plots of whole sim along %s axis'%dim
    field = 'particle_mass'

    #Make plot
    if dim == 'x':
        output = yt.ParticlePlot(ds, 'particle_position_y', 'particle_position_z', field)
    elif dim == 'y':
        output = yt.ParticlePlot(ds, 'particle_position_z', 'particle_position_x', field)
    elif dim == 'z':
        output = yt.ParticlePlot(ds, 'particle_position_x', 'particle_position_y',field)
    else:
        raise ValueError("Function part_mass_whole received an unknown dimension '%s'. Use 'x','y','z'. This is case sensitive"%(dim))
    #Save object to disk
    #output.frb.save_as_dataset('yt_%s_%s_part_%s_%s'%(simtype,curr_src,field,dim))

    #Format color bar and axis units.
    #output.set_zlim(('all','particle_mass'), all_particle_mass_extrema[0],all_particle_mass_extrema[1])
    output.annotate_timestamp(corner='upper_left', time=False,redshift=True, draw_inset_box=False)
    output.annotate_scale(corner='upper_right')
    output.set_axes_unit('Mpccm/h')
    output.set_unit(field, 'Msun')

    #Save figure to disk
    output.save('%s_%s_part_%s_%s.png'%(simtype,curr_src,field,dim))
 
    
    #Return the plot object for future use.
    return output

#Density projection whole simulation
def dens_whole_sim(ds,dim,simtype,curr_src):
    print datetime.now().strftime('%H:%M:%S'),'Making density projection plots of the whole simulation along %s axis'%dim
    fieldtype = 'gas'
    field = 'density'

    #make projection and save to disk
    #proj = ds.proj(field, dim, weight_field=field)
    #proj.save_as_dataset('yt_%s_%s_%s_proj_%s'%(simtype,curr_src,field,dim))
    #proj = yt.load('yt_%s_%s_%s_proj_%s.h5'%(simtype,curr_src,field,dim))
    #make a plot with appropriate units and save to disk
    output = yt.ProjectionPlot(ds, dim, field)
    output.annotate_timestamp(corner='upper_left', time=False,redshift=True, draw_inset_box=False)
    output.annotate_scale(corner='upper_right')
    output.set_axes_unit('Mpccm/h')
    output.set_unit(field, 'Msun/kpc**2')

    #save figure to disk
    output.save('%s_%s_%s_proj_%s.png'%(simtype,curr_src,(fieldtype,field),dim))
    return output

#Temperature projection whole sim
def temp_whole_sim(ds,dim,simtype,curr_src):
    print datetime.now().strftime('%H:%M:%S'),'Making temperature Projection plots of the whole simulation along %s axis'%dim
    field = 'temperature'

    #make projection and save to disk
    #proj = ds.proj(field, dim, weight_field=field)
    #proj.save_as_dataset('yt_%s_%s_%s_proj_%s'%(simtype,curr_src,field,dim))

    #make a plot with appropriate units and save to disk
    output = yt.ProjectionPlot(ds, dim, field)
    output.annotate_timestamp(corner='upper_left', time=False,redshift=True, draw_inset_box=False)
    output.annotate_scale(corner='upper_right')
    output.set_axes_unit('Mpccm/h')
    output.set_unit(field, 'K*Mpc')
    #output.set_zlim(('gas','density'), gas_density_extrema[0],gas_density_extrema[1])
    output.save('%s_%s_%s_proj_%s.png'%(simtype,curr_src,field,dim))
    return output


###HALO PLOTS###

#Temperature Projection through a halo. Stops at the edge of the halo
def temp_plot_halo(ds,dim,cen,wid,simtype,curr_src,curr_halo):
    print datetime.now().strftime('%H:%M:%S'),'Making temp projection plots of halo #%d along %s axis'%(curr_halo,dim)
    fieldtype = 'gas'
    field = 'temperature'

    #make projection and save to disk
    #proj = ds.proj(field, dim, weight_field=field,center=cen,width=wid)
    #proj.save_as_dataset('yt_%s_%s_%d_proj_%s_%s.png'%(simtype,curr_src,curr_halo,field,dim))

    #make a plot with appropriate units and save to disk
    output = yt.ProjectionPlot(ds, dim, (fieldtype,field), center=cen, width=wid)

    output.annotate_timestamp(corner='upper_left', time=False,redshift=True, draw_inset_box=False)
    #output.annotate_scale(corner='upper_right')
    output.set_axes_unit('Mpccm/h')
    output.set_unit(field, 'K*Mpc')
    output.save('%s_%s_%d_proj_%s_%s.png'%(simtype,curr_src,curr_halo,field,dim))
    return output

#Presssure Projection through a halo. Stops at the edge of the halo
def pres_plot_halo(ds,dim,cen,wid,simtype,curr_src,curr_halo):
    print datetime.now().strftime('%H:%M:%S'),'Making pressure projection plots of halo #%d along %s axis'%(curr_halo,dim)
    fieldtype = 'gas'
    field = 'pressure'

    #make projection and save to disk
    #proj = ds.proj(field, dim, weight_field=field,center=cen,width=wid)
    #proj.save_as_dataset('yt_%s_%s_%d_proj_%s_%s.png'%(simtype,curr_src,curr_halo,field,dim))

    #make a plot with appropriate units and save to disk
    output = yt.ProjectionPlot(ds, dim, (fieldtype,field), center=cen, width=wid)

    output.annotate_timestamp(corner='upper_left', time=False,redshift=True, draw_inset_box=False)
    #output.annotate_scale(corner='upper_right')
    output.set_axes_unit('Mpccm/h')
    #output.set_unit(field, 'K*Mpc')
    output.save('%s_%s_%d_proj_%s_%s.png'%(simtype,curr_src,curr_halo,field,dim))
    return output



#Density Projection through a halo.
def proj_plot_halo(ds,dim,cen,wid,simtype,curr_src,curr_halo):
    print datetime.now().strftime('%H:%M:%S'),'Making density projection plots of halo #%d along %s axis'%(curr_halo,dim)
    fieldtype = 'gas'
    field = 'density'
    output = yt.ProjectionPlot(ds, dim, (fieldtype,field), center=cen, width=wid)

    output.annotate_timestamp(corner='upper_left', time=False,redshift=True, draw_inset_box=False)
    #output.annotate_scale(corner='upper_right')
    output.set_axes_unit('Mpccm/h')
    output.set_unit(field, 'Msun/kpc**2')
    output.save('%s_%s_%d_proj_%s_%s.png'%(simtype,curr_src,curr_halo,field,dim))
    return output


#Cube region around a halo (test)
def cube_proj_halo(ds,dim,cen,wid,simtype,curr_sim,curr_halo,cube):
    print datetime.now().strftime('%H:%M:%S'),'Making density Projection plots of halo #%d along %s axis with cube source'%(curr_halo,dim)
    field = 'density'
    output = yt.ProjectionPlot(ds, dim, field, center=cen, width=wid, data_source=cube)

    output.annotate_timestamp(corner='upper_left', time=False,redshift=True, draw_inset_box=False)
    #output.annotate_scale(corner='upper_right')
    output.set_axes_unit('Mpccm/h')
    output.set_unit(field, 'Msun/kpc**2')
    output.save('%s_%s_%d_proj_%s_%s_cube.png'%(simtype,curr_sim,curr_halo,field,dim))
    return output
  
#Density Slice plots of each halo
def slic_halo(ds,dim,cen,wid,simtype,curr_sim,curr_halo):
    print datetime.now().strftime('%H:%M:%S'),'Making density Slice plots of halo #%d along %s axis'%(curr_halo,dim)
    field = 'density'
    output = yt.SlicePlot(ds, dim, field, center=cen, width=wid)

    #output.frb.save_as_dataset('yt_%s_%s_%d_slice_%s_%s.png'%(simtype,curr_sim,curr_halo,field,dim))

    output.annotate_timestamp(corner='upper_left', time=False,redshift=True, draw_inset_box=False)
    #output.annotate_scale(corner='upper_right')
    output.set_axes_unit('Mpccm/h')
    output.set_unit(field, 'Msun/kpc**3')
    output.save('%s_%s_%d_slice_%s_%s.png'%(simtype,curr_sim,curr_halo,field,dim))
    return output

#Particle mass plot of each halo in each dim
def part_mass_halo(ds,dim,cen,wid,simtype,curr_sim,curr_halo):
    """
    usage: part_mass_halo(ds,wid,cen,simtype,curr_sim,curr_halo,field)\
           ds: dataset\
           wid: width of halo\
           cen: center of halo (x,y,z) in code units\
           simtype: simulation type (used for file formating)"
    """
    print datetime.now().strftime('%H:%M:%S'),'Making Particle mass Projection plots of halo #%d along %s axis'%(curr_halo,dim)
    weightfield = 'particle_mass'

    if dim == 'x': 
        output = yt.ParticlePlot(ds, 'particle_position_y', 'particle_position_z', weightfield, width=wid,center=cen)
    elif dim == 'y': 
        output = yt.ParticlePlot(ds, 'particle_position_z', 'particle_position_x', weightfield, width=wid,center=cen)
    elif dim == 'z': 
        output = yt.ParticlePlot(ds, 'particle_position_x', 'particle_position_y', weightfield, width=wid,center=cen)
    else:
        raise ValueError("Function part_mass_halo received an unknown dimension '%s'. Use 'x','y','z'. This is case sensitive"%(dim))

    #Format pictures with nice units and color map.
    output.annotate_timestamp(corner='upper_left', time=False,redshift=True, draw_inset_box=False)
    #output.annotate_scale(corner='upper_right')
    output.set_axes_unit('Mpccm/h')
    output.set_unit(weightfield, 'Msun')
    #output.set_cmap('particle_mass','RdBu_r')

    #Save with filename corresponding to dimension
    output.save('%s_%s_%d_%s_%s.png'%(simtype,curr_sim,curr_halo,weightfield,dim))

    return output



#Radial gas density profiles of halo
def prof_density_halo(ds,sph,simtype,curr_sim,curr_halo):
    "This function makes a profile of the halo chosen, saves it, and returns the profile as an object for future use"
    print datetime.now().strftime('%H:%M:%S'),'Making gas Density Profile Plot of halo #%d'%(curr_halo)
    fieldtype = 'gas'
    field = 'density'
    weight = None
    
    rp    = yt.create_profile(sph,'radius',(fieldtype,field),weight_field=weight)
    plot  = yt.ProfilePlot(sph,'radius',(fieldtype,field),weight_field=weight)

    plot.set_unit('radius', 'Mpccm/h')
    plot.set_unit(field, 'Msun/kpc**3')
    plot.set_log('radius',False)
    plot.save('%s_%s_%d_prof_%s.png'%(simtype,curr_sim,curr_halo,field))
    #Return the profile object for future use.
    return rp

#Radial DM density profiles of halo
def DM_dens_prof(ds,sph,simtype,curr_sim,curr_halo):
    print datetime.now().strftime('%H:%M:%S'),'Making DM Density Profile Plot of halo #%d'%(curr_halo)
    field_type = 'deposit'
    field = 'io_density'
    weight= None


    rp  = yt.create_profile(sph, 'radius', (field_type,field),weight_field=weight)
    output = yt.ProfilePlot(sph, "radius", (field_type,field),x_log=None,y_log=None)

    output.set_unit('radius', 'Mpccm/h')
    output.set_unit(field, 'Msun/kpc**3')
    output.set_log('radius',False)
    output.save('%s_%s_%d_prof_%s.png'%(simtype,curr_sim,curr_halo,field))
    return rp

#Radial gas density profiles of halo
def prof_temp_halo(ds,sph,simtype,curr_sim,curr_halo):
    "This function makes a profile of the halo chosen, saves it, and returns the profile as an object for future use"
    print datetime.now().strftime('%H:%M:%S'),'Making gas temperature Profile Plot of halo #%d'%(curr_halo)
    fieldtype = 'gas'
    field = 'temperature'
    weight= None

    rp    = yt.create_profile(sph,'radius',(fieldtype,field),weight_field=weight)
    plot  = yt.ProfilePlot(sph,'radius',(fieldtype,field),weight_field=weight)

    plot.set_unit('radius', 'Mpccm/h')
    plot.set_unit(field, 'K')
    plot.set_log('radius',False)
    plot.save('%s_%s_%d_prof_%s.png'%(simtype,curr_sim,curr_halo,field))
    #Return the profile object for future use.
    return rp

#########################FUNCTIONS END###########################################

