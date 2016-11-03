# Lissajous curves. A family of parametric curves of the type
#   x = A sin(at+d)
#   y = B sin(bt)
# More info: https://en.wikipedia.org/wiki/Lissajous_curve
#
# Author:    Kjetil Andre Johannessen
# Institute: Norwegian Univeristy of Science and Technology (NTNU)
# Date:      October 2016
#


from sys import path
path.append('../')
from splipy import *
from fractions import gcd
import splipy.curve_factory as curves
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from math import pi

def lissajous(a, b, d):
  # request a,b integers, so we have closed, periodic curves
  n = gcd(a,b)
  N = (a/n) * (b/n) # number of periods before looping

  # error test input
  if N > 1e4:       # non-integer (a,b) or otherwise too irregular
    raise Exception('Non-periodic', 'a,b must be integers (of moderate size)')

  # compute a set of interpolation points
  t = np.linspace(0,2*pi/n, 3*N) # using 3N interpolation points is decent enough
  x = np.array([np.sin(a*t + d), np.sin(b*t)])

  # do a cubic curve interpolation with periodic boundary conditions
  return curves.cubic_curve(x.T, curves.Boundary.PERIODIC)


### main program ###

# create the curve
crv = lissajous(60, 44, pi/2);
# crv = lissajous(3, 4, pi/2);

# evalute the curve at n points
n = 3000
t = np.linspace(crv.start(0), crv.end(0), n);
x = crv(t)



### do the  plotting animation
fig = plt.figure()
lines = []
frames = 100        # number of frames in animation
fps    = n / frames # logical error if fps is not an integer
for i in range(frames):
    j = range(i*fps, min((i+1)*fps+1, n-1))
    l, = plt.plot(x[j,0], x[j,1], color=[1,1,1])
    lines.append(l)
plt.axis('equal')
plt.axis([-2, 2, -2, 2])

def animate(i):
    m = len(lines)
    # first sweep, make all the lines black (and thus visible)
    if i<m:
        lines[i].set_color([0.0,0.0,0.0])
        lines[i].set_zorder(2*m-i);
    # second sweep, make all lines white again
    else:
        lines[i-m].set_color([1.0,1.0,1.0])
        lines[i-m].set_zorder(i-m);

    
# create and show the animation
ani = animation.FuncAnimation(fig, animate, np.arange(0,2*n/fps), interval=10)
plt.show()

# save results as an animated gif for web display (PS: this function call is slow)
# ani.save('lissajous34.gif', writer='imagemagick', fps=30);
