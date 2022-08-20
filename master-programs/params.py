"""
This is the parameter file used by main.py
It contains all the various parameters specifying which simulations are placed where,and the paths to them.

Made by Thor Andreas Seiff Ellewsen summer 2016.

Updated 13/12/16
"""

datapath = []

#NEWEST SIMS
#            [['folder_name',halomatchfile_colnumber,parent_folder_name,Halofile_file_name_to_use],[],[]]
simtype    = [['vanilla',0,'LCDM','LCDM_vanilla'],['vanilla',1,'SymA','SymA_vanilla'],
              ['vanilla',2,'SymB','SymB_vanilla'],['vanilla',3,'SymC','SymC_vanilla'] ,['vanilla',4,'SymD','SymD_vanilla'],
              ['cool_star',5,'LCDM','LCDM_cool_star'],['cool_star',6,'SymA','SymA_cool_star'],
              ['cool_star',7,'SymB','SymB_cool_star'],['cool_star',8,'SymC','SymC_cool_star'],['cool_star',9,'SymD','SymD_cool_star'],
              ['cool_star_SN',10,'LCDM','LCDM_cool_star_SN'],['cool_star_SN',11,'SymA','SymA_cool_star_SN'],['cool_star_SN',12,'SymB','SymB_cool_star_SN'],
              ['cool_star_SN',13,'SymC','SymC_cool_star_SN'],['cool_star_SN',14,'SymD','SymD_cool_star_SN'],
              ['DM',15,'LCDM','LCDM_DM'],['DM',16,'SymA','SymA_DM'],['DM',17,'SymB','SymB_DM'],
              ['DM',18,'SymC','SymC_DM'],['DM',19,'SymD','SymD_DM']]

sourcelist = [11]#None

#Set path to sim data.
for src in sourcelist:
    for sim in simtype:
        #datapath.append('/mn/stornext/d5/astcosim/AmirBridgetDavid/Data_From_Hexagon/Runs/%s/%s/output_000??/info_000??.txt'%(sim[2],sim[0]))
        datapath.append('/mn/stornext/d5/astcosim/AmirBridgetDavid/Data_From_Hexagon/Runs/%s/%s/output_000%s/info_000%s.txt'%(sim[2],sim[0],src,src))


#halodir    = 'Halos/'
#halomatchdir    = 'halomatch/Runs/'

################Load matching list########################
#match = np.empty(len(simtype),dtype=object)
#
#
#itersims = iter(simtype)
#skip first element
#next(itersims)
#for i,element in enumerate(itersims):
#    if i == 0:
#        match[i] = np.loadtxt('%s''match_LCDM_vanilla_%s.txt'%(halomatchdir,element[3]),dtype=object,skiprows=3,usecols=0,unpack=True)
#        match[i+1] = np.loadtxt('%s''match_LCDM_vanilla_%s.txt'%(halomatchdir,element[3]),dtype=object,skiprows=3,usecols=1,unpack=True)
#    else:
#        match[i+1] = np.loadtxt('%s''match_LCDM_vanilla_%s.txt'%(halomatchdir,element[3]),dtype=object,skiprows=3,usecols=1,unpack=True)
#

#######Load halocatalog##############
    
#halos = np.zeros(len(simtype),dtype=object)
#haloradius = np.zeros(len(simtype),dtype=object)
#for i,sim in enumerate(simtype):
#    halos[i] = np.loadtxt('%s''%s.txt'%(halodir,sim[3]),skiprows=1,unpack=False,usecols=(1,2,3))
#    haloradius[i] = np.loadtxt('%s''%s.txt'%(halodir,sim[3]),skiprows=1,usecols=(7))

# #Scale halos to code units
# halos      = halos/64.      #In Mpc/h
# haloradius = haloradius/64. #In Mpc/h

#Here you insert the INDEXES of halos in halos/haloradius corresponding to the simtypes chosen.
#This must either be a an array with each row corresponding to the halos in each sim like the example commented out below or None.
#If left as None, the main program will fill the array on its own with the 2 largest halos in one of the sims and the corresponding halos in the others.
#(Assuming of course that the match.txt file is formatted the right way.
#halolist = None

#halolist = np.zeros((numsims,numhalos),dtype=int)
#for k,simname in enumerate(simtype):
#    for i in range(2):
#        halolist[k][i] = match[i][simname[1]] #load halos for right sim
