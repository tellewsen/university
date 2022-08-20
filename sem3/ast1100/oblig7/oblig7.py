from scitools.all import *

x_arc,y_arc,r,wl=loadtxt("galaxies.txt",unpack=True)

c      = 3.0  *10**8
wl0    = 21.2 *10**-2
G      = 6.67 *10**-11
Mpctom = 3.08567758 * 10**22 	#Megaparsec to meters 

#Doppler function
def Doppler(wl):
    return c*((wl - wl0)/wl0)

Vp = Doppler(wl)	#Array of all the velocities
V_cluster = mean(Vp)	#Speed of cluster

print "Speed of Cluster: ",V_cluster, "m/s"

#Plot the position of the galaxies in the cluster
plot(x_arc,y_arc,'o')
xlabel("X-axis [arcminutes]")
ylabel("Y-axis [arcminutes]")
hardcopy("oblig7.png")

#Makes a vector r_vector with x,y,z values instead of just length from a vector r
r_vector = zeros((len(r),3))
for i in range(0,len(r)):
    r_vector[i,0] = r[i]*sin((x_arc[i]*pi)/(60*180))*Mpctom 
    r_vector[i,1] = r[i]*sin((y_arc[i]*pi)/(60*180))*Mpctom
    r_vector[i,2] = r[i]*cos((x_arc[i]*pi)/(60*180))*Mpctom

#Calculates length of a vector
def length(r):
    return sqrt(r[0]**2 + r[1]**2 + r[2]**2)

#Calculates the mass of one galaxy in the cluster (This is the formula we found analytically)
v_sum = 0
r_sum = 0
for i in range(0,len(r)):
    v_sum += (Vp[i]-V_cluster)**2
    for j in range(0,len(r)):
        if i<j:
            r_ij = length(r_vector[i,:]-r_vector[j,:])
            r_sum += 1/r_ij
m_galaxy = v_sum/(G*r_sum)		#Mass of one galaxy in the cluster


m_cluster = m_galaxy*100 		#There were 100 galaxies in the cluster
m_luminous = 4E43 			#This was found in the lecture notes
Ratio = m_cluster/m_luminous 		#Calculates the percentage of luminous matter in the cluster
print "Mass of cluster: ", m_cluster
print "Mass of cluster divides by luminous mass: ",Ratio,
