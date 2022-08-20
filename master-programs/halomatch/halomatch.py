"""
This program matches halos from two different simulations together.
It accepts the output format of the AMIGA HALO FINDER.

This is a slightly modified version based on the halo matcher written in C by Bridget Falck.
It first finds the halos within the mass threshold, and then finds the closest of those to the original halo.

example run:
python halomatch.py all_LCDM_halos_z0 /mn/stornext/d5/astcosim/AmirBridgetDavid/HaloFinder/HaloList/ all_Sym_D_halos_z0 64 .5 match

Latest edit: 24th october 2016 by Andreas Ellewsen
"""

import sys
import numpy as np
if (len(sys.argv) != 7):
    print "Wrong number of arguments.\n"
    print "arg1: directory of both halo files\n"
    print "arg2: filename of reference halo file\n"
    print "arg3: filename of alternate halo file\n"
    print "arg4: box size in Mpc/h\n"
    print "arg5: mass difference threshold\n"
    print "arg6: label for output match file\n"
    sys.exit(0)

outdir   = sys.argv[1]
label    = sys.argv[2]
altlabel = sys.argv[3]

catfile  = '%s%s.txt'%(outdir,label)
altfile  = '%s%s.txt'%(outdir,altlabel)


try:
    boxsize = float(sys.argv[4])
except ValueError:
    print "That's no boxsize; try again.\n"
    sys.exit(0)

b2    =  boxsize/2.
negb2 = -boxsize/2.

try:
    mthresh = float(sys.argv[5])
except ValueError:
    print "That's no mass threshold; try again.\n"
    sys.exit(0)    

mlabel    = sys.argv[6];
matchfile = '%s.txt'%(mlabel)


"""
#Need: halo centers, halo masses, npart, make sure centers in Mpc
#   Possibly hard code: mass threshold, min and max mass, min npart 
#   Output ORG index, alt index, and webenv of matched halos
"""

# Define mass and npart thresholds (if not reading from stdin)

massmin = 1.e11
massmax = 1.e15
#mthresh = 0.1; // abs(mass_org - mass_alt)/mass_alt < mthresh
npmin = 100


#Read both catalogs

try:
    nmem1,x1,y1,z1,mvir1 = np.loadtxt(catfile,unpack=True,skiprows=1,usecols=(0,1,2,3,8))
    pos1 = np.array([x1,y1,z1])
except IOError:
    print 'Unable to open LCDM catalog file %s.\n'%(catfile)
    sys.exit(0)


try:
    nmem2,x2,y2,z2,mvir2 = np.loadtxt(altfile,unpack=True,skiprows=1,usecols=(0,1,2,3,8))
    pos2 = np.array([x2,y2,z2])
except IOError:
    print 'Unable to open Sym D catalog file %s.\n'%(altfile)
    sys.exit(0)

"""
#For each ALT halo:
#   - Check that mass and nmem within acceptible range
#   - For each ORG halo:
#     -- If m200 == 0, break
#     -- calculate PBC distance between halos
#     -- save minimum distance and ORG index
#   - Check that minimum distance reasonable (< 1 Mpc)
#   - Check that mass diff. within threshold
#   - Check that nmem of ORG index halo acceptable
#   - If checks pass:
#     -- save ALT index, ORG index, and webenv
#     -- update MATCH index and nhm
"""


#Initialize final arrays to size nha (max possible) since don't know 
#nhm yet, the final number of halos in the matched catalog
nho = len(nmem1) #number of halos original list
nha = len(nmem2) #number of halos alternate list

#iorg = np.empty(nho,dtype=object)#(int *)malloc(nha*sizeof(int));
iorg = np.arange(nho)
ialt = np.empty(nho,dtype=object)#(int *)malloc(nha*sizeof(int));
ialt[:] = None
#envmatch = (unsigned char *)malloc(nha*sizeof(unsigned char));

nhm = 0;
svec = np.zeros(3)
can_num = 0
can_s   = np.zeros(nha)
can_j   = np.zeros(nha)
s_min = 1.

for i in range(nho):
  print 'i =%d/%d'%(i,nho-1)
  if ((mvir1[i] >= massmin) and (mvir1[i] <= massmax)): #check halo in mass range
    #print mvir2[i]
    if (nmem1[i] > npmin): #check more particles than minimum
      for j in range(nha): #run through alternate halos
        if (mvir2[j] != 0): #Check non-zero mass
          if (abs(mvir2[j]/mvir1[i] - 1.) <= mthresh): #Check masses within threshold
            #Make a list of all candidates within mass threshold
            for k in range(3):
              svec[k] = pos2[k][j] - pos1[k][i] #distance between org and alt
              if(svec[k] <= negb2): #Consider PBCs
                svec[k] += boxsize
              if(svec[k] > b2): #Consider PBCs
                svec[k] -= boxsize
            can_s[can_num] = np.sqrt(svec[0]*svec[0]+svec[1]*svec[1]+svec[2]*svec[2]) #save candidate distance from original
            can_j[can_num] = j #save candidate index
            can_num +=1 

  #Run through candidates and pick the one closest to original halo
  #save index of org
  iorg[i] = i
  for element in range(np.count_nonzero(can_s)): 
    #Find closest to org
    if(can_s[element] < s_min and can_s[element] < 1.): 

      #save index of alt
      ialt[i] = can_j[element] 
      #reduce min to find closest
      s_min = can_s[element] 
      match = 1 #register that a match was found

  #Keep count of the total number of matches
  if match > 0:
      #Increase number of matches by 1
      nhm += 1
      print nhm,"matches found"

  #return to start for next halo
  match    = 0
  s_min    = 1.
  can_j[:] = 0
  can_s[:] = 0
  can_num  = 0
  

print "%d matches total\n"%nhm

# Output
match = open(matchfile,"w")
match.write('%d matches\n'%nhm)
match.write('%s\n'%(sys.argv))
match.write('%s\t%s\n'%(label,altlabel))
for i in range(nho):
  if ialt[i] == None:
    match.write('%d\tNone\n'%(iorg[i]))
  else:
    match.write('%d\t%d\n'%(iorg[i],ialt[i]))
match.close()
#for (i=0; i<nhm; i++) fwrite(&envmatch[i],1,sizeof(unsigned char),match);





#check that the x coordinates of each of the halos are roughly the same.
