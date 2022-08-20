import yt
import numpy as np


#Defining particle filters for RAMSES
def Stars(pfilter,data):
    age  = data[(pfilter.filtered_type,"particle_age")]
    filter = np.logical_or(age < 0,age > 0)
    return filter

def DM(pfilter,data):
    filter = data[(pfilter.filtered_type,"particle_age")] == 0
    return filter

yt.add_particle_filter("stars",function=Stars,filtered_type='all', 
                    requires=['particle_age'])

yt.add_particle_filter("dm",function=DM,filtered_type='all', 
                    requires=['particle_age'])



#datapath = 'output_00080/info_00080.txt'
datapath = 'output_00016/info_00016.txt'
ds = yt.load(datapath)






ds.add_particle_filter("stars")
ds.add_particle_filter("dm")
ds.derived_field_list

ad = ds.all_data()
#ages = ad['particle_age']
#ids  = ad['particle_identifier']
#original_particles_ids = ids[ages == 0]

age = ad['age']
minage = np.min(age)
print sum(age<0)
print sum(age>0)
print sum(age==0)
print sum(age!=0)
print len(age)
