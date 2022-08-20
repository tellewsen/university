#!/usr/bin/env python
# Tools for sparse matrices
import scipy.sparse as sparse
import scipy.sparse.linalg

# Numerical tools
from numpy import *

# Plotting library
from matplotlib.pyplot import *

"""Physical constants"""
_u  	= 931.494061			# Atomic mass unit [MeV/c^2]
_c 	= 3e5          			# Speed of light [nm / ps]
_p 	= 938.27        		# Rest energy for a proton [MeV]
_E0Cl35 = 34.97*_u*_c**2		# Rest energy for CL35 [MeV]
_E0p    = (_p*_E0Cl35)/(_p+_E0Cl35)	# Rest energy for reduced mass [MeV]
_hbarc 	= 197.32e-6   			# [MeV nm]
_omega 	= 564 				# [1/ps]
_x0	= 0.127				# [nm]
def Psi0( x ):
    """
    Initial state for a stationairy gaussian wave packet.
    """
    x1 	= 0.240 # [nm] Starts at center.
    a 	= 0.004 # [nm]

    A = ( 1. / ( 2 * pi * a**2 ) )**0.25
    K1 = exp( - ( x - x1 )**2 / ( 4. * a**2 ) )

    return A * K1

def harmonicOscillator( x ):
    """
    The potential for a quantum harmonic oscillator for a proton.

    @param k The wave number. Default value chosen by eye.
    """
    potential = 0.5*(_E0p/(_c*_c))*_omega*_omega*(x-_x0)*(x-_x0)

    return potential

def Morse( x ):
    """
    The Morse potential for a quantum harmonic oscillator for a proton.

    @param k The wave number. Default value chosen by eye.
    """
    V0 = 10.2e-6 #[MeV]
    b  = 12.6    #[nm^-1]
    potential = V0*(exp(-2*b*(x-_x0)) -2*exp(-b*(x-_x0)))

    return potential

if __name__ == '__main__':
    nx = 1001  # Number of points in x direction
    dx = 0.001 # Distance between x points [pm]

    # Use zero as center, same amount of points each side
    a = - 0.5 * nx * dx
    b =   0.5 * nx * dx
    x = linspace( a, b, nx )
    
    # Time parameters
    T  = 1.4 	# How long to run simulation [ps]
    dt = 1e-4 	# The time step [ps]
    t  = 0
    time_steps = int( T / dt ) # Number of time steps
    # Constants - save time by calculating outside of loop
    k1 = - ( 1j * _hbarc * _c) / (2. * _E0p )
    k2 =   ( 1j * _c ) / _hbarc

    # Create the initial state Psi
    Psi = Psi0(x)

    # Create the matrix containing central differences. It it used to
    # approximate the second derivative.
    data = ones((3, nx))
    data[1] = -2*data[1]
    diags = [-1,0,1]
    D2 = k1 / dx**2 * sparse.spdiags(data,diags,nx,nx)

    # Identity Matrix
    I = sparse.identity(nx)

    # Create the diagonal matrix containing the potential.
    V_data = Morse(x)
    V_diags = [0]
    V = k2 * sparse.spdiags(V_data, V_diags, nx, nx)
    
    # Put mmatplotlib in interactive mode for animation
    ion()
    
    # Setup the figure before starting animation
    fig = figure() # Create window
    ax = fig.add_subplot(111) # Add axes
    line, = ax.plot( x, abs(Psi)**2, label='$|\Psi(x,t)|^2$' ) # Fetch the line object

    # Also draw a green line illustrating the potential
    ax.plot( x, V_data, label='$V(x)$' )

    # Add other properties to the plot to make it elegant
    fig.suptitle("Solution of Schrodinger's equation with harmonic oscillator") # Title of plot
    ax.grid('on') # Square grid lines in plot
    ax.set_xlabel('$x$ [nm]') # X label of axes
    ax.set_ylabel('$|\Psi(x, t)|^2$ [1/nm] and $V(x)$ [MeV]') # Y label of axes
    ax.legend(loc='best') # Adds labels of the lines to the window
    axis([-0.01, 0.5, -10, 100])
    #draw() # Draws first window
    
    #Extra stuff to calculate expectation value for position and moment
    tlist        = []
    expectX      = []
    expectX2     = []
    expectP      = []

    counter      = 0

    # Time loop
    while t < T:
        """
        For each iteration: Solve the system of linear equations:
        (I - k/2*D2) u_new = (I + k/2*D2)*u_old
        """
        # Set the elements of the equation
    	A = (I - dt/2. * (D2 + V))
    	b = (I + dt/2. * (D2 + V)) * Psi

        # Calculate the new Psi
    	Psi = sparse.linalg.spsolve(A,b)

	# Calculate expectation values
	expectX.append(trapz(x*abs(Psi)**2,dx=dx))
	expectX2.append(trapz((x**2)*abs(Psi)**2,dx=dx))
	tlist.append(t)

        # Update time
    	t += dt
	counter +=1
        print float(counter)/float(time_steps)*100,'% done'
        #savefig('tmp%04d.png'%counter)

    	# Plot this new state
        #line.set_ydata( abs(Psi)**2 ) # Update the y values of the Psi line
        #draw() # Update the plot


    # Turn off interactive mode
    ioff()
    
    # Add show so that windows do not automatically close
    #show()
    
    # Calculate expectation values and deviation
    tlist     = array(tlist)
    expectX   = array(expectX)
    expectX2  = array(expectX2)
    expectP   = []
    sigmaX    = sqrt(expectX2 - expectX**2)
    
    for i in range(len(expectX)-1):
        expectP.append(_E0p/(_c*_c)*(expectX[i+1] - expectX[i]) /  dt)
    
    
    figure()
    plot(tlist,expectX)
    xlabel('Time[ps]')
    ylabel('Expectation value of position [nm]')
    title('Expectation value of position vs time')
    legend(['Quantum'])
    
    figure()
    plot(tlist,sigmaX)
    xlabel('Time[ps]')
    ylabel('Standard deviation [nm]')
    title('Standard deviation of position vs time')
    
    figure()
    plot(tlist[:-1],expectP)
    xlabel('Time [ps]')
    ylabel('Expectation value of momentum [MeV/c]')
    title('Expectation value of momentum vs time')
    show()
    """
    from scitools.easyviz import *
    movie('tmp*.png')
    """
