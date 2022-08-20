#Exercise 5.21 p.247
from scitools.std import *
import time

def f(x,t):
    return exp(-(x-3*t)**2)*sin(3*pi*(x-t))


t_start = -1
t_stop = 1
t_values = linspace(t_start, t_stop, 10)
x = linspace(-6,6, 1000)

# Show the movie on the screen
# and make hardcopies of frames simultaneously.

counter = 0
for t in t_values:
    y = f(x,t)
    plot(x, y, axis=[x[0], x[-1], -1.1, 1.1],
        xlabel='x', ylabel='f(x,s)', legend='t=%4.2f' % t,
        savefig='tmp%04d.png' % counter)
    counter += 1
    time.sleep(0.2) #pause to control movie speed

# Make movie file the simplest possible way:
movie('tmp*.png')

"""
Terminal> python plot_wavepacket_movie.py 
Could not find/open font when opening font "arial", using internal non-scalable font
Could not find/open font when opening font "arial", using internal non-scalable font
Could not find/open font when opening font "arial", using internal non-scalable font
Could not find/open font when opening font "arial", using internal non-scalable font
Could not find/open font when opening font "arial", using internal non-scalable font
Could not find/open font when opening font "arial", using internal non-scalable font
Could not find/open font when opening font "arial", using internal non-scalable font
Could not find/open font when opening font "arial", using internal non-scalable font
Could not find/open font when opening font "arial", using internal non-scalable font
Could not find/open font when opening font "arial", using internal non-scalable font
Could not find/open font when opening font "arial", using internal non-scalable font
Could not find/open font when opening font "arial", using internal non-scalable font
Could not find/open font when opening font "arial", using internal non-scalable font
Could not find/open font when opening font "arial", using internal non-scalable font
Could not find/open font when opening font "arial", using internal non-scalable font
Could not find/open font when opening font "arial", using internal non-scalable font
Could not find/open font when opening font "arial", using internal non-scalable font
Could not find/open font when opening font "arial", using internal non-scalable font
Could not find/open font when opening font "arial", using internal non-scalable font
Could not find/open font when opening font "arial", using internal non-scalable font



Found 10 files of the format tmp*.png.

scitools.easyviz.movie function runs the command: 

convert -delay 4 tmp0000.png tmp0001.png tmp0002.png tmp0003.png tmp0004.png tmp0005.png tmp0006.png tmp0007.png tmp0008.png tmp0009.png movie.gif



movie in output file movie.gif
"""

