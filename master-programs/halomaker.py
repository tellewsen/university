#This file converts the AHF halocatalog to a yt compatible format loadable by using
#ds_halo = yt.load('test.0.h5')

import numpy as np

#Load AHF halo catalog
halodata = np.genfromtxt('../halomatch/Runs/LCDM_vanilla.txt',skip_header=1,unpack=True,max_rows=4) #remove [:,0] when this works
#Load the simulation the halo belongs to
ds = load('/mn/stornext/d5/astcosim/AmirBridgetDavid/Data_From_Hexagon/Runs/LCDM/vanilla/output_00011/info_00011.txt')

#num,x,y,z,vx,vy,vz,rvir,mvir,mgas,....
npart = halodata[0]
x = halodata[1]
y = halodata[2]
z = halodata[3]
vx= halodata[4]
vy= halodata[5]
vz= halodata[6]
rvir  = halodata[7]
mvir  = halodata[8]

#assign units to everything
x = ds.arr(x,'Mpc/h')
y = ds.arr(y,'Mpc/h')
z = ds.arr(z,'Mpc/h')
vx = ds.arr(vx,'km/s')
vy = ds.arr(vy,'km/s')
vz = ds.arr(vz,'km/s')
rvir = ds.arr(rvir, 'Mpc/h')
mvir = ds.arr(mvir,'Msun/h')
n_halos = rvir.size
partid = np.linspace(1,n_halos,n_halos,dtype=int)

extra_attrs = {"data_type": "halo_catalog",
               "num_halos": n_halos}

data = {"particle_identifier": partid,
		"particle_number": npart,
        "particle_position_x": x,
        "particle_position_y": y,
        "particle_position_z": z,
        "particle_velocity_x": vx,
        "particle_velocity_y": vy,
        "particle_velocity_z": vz,
        "virial_radius": rvir,
        "particle_mass": mvir}
ftypes = dict([(field, ".") for field in data])
save_as_dataset(ds,'test.0.h5', data, field_types=ftypes, extra_attrs=extra_attrs)
