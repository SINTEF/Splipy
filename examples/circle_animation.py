# Animation to show the NURBS circle parametrization. It is not an arc-length
# parametrization, and the acceleration vector is shown to be discontinuous.
#
# Author:    Kjetil Andre Johannessen
# Institute: Norwegian University of Science and Technology (NTNU)
# Date:      March 2016
#
from sys import argv
from splipy import *
import splipy.curve_factory as curves
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


INTERACTIVE = "--ci" not in argv[1:]


n = 250                                  # number of evaluation points
c = curves.circle()                      # create the NURBS circle
t = np.linspace(c.start(0), c.end(0), n) # parametric evaluation points
x = c(t)                                 # physical (x,y)-coordinates, size (n,2)
v = c.derivative(t, 1)                   # velocity at all points
a = c.derivative(t, 2)                   # acceleration at all points

# plot the circle and get reference to the acceleration/velocity lines which we
# will update during animation
fig = plt.figure(figsize=[10, 10])
plt.plot(x[:,0], x[:,1], 'k-')
velocity,     = plt.plot([x[0,0], x[0,0]+v[0,0]], [x[0,1], x[0,1]+v[0,1]], 'r-', linewidth=2)
acceleration, = plt.plot([x[0,0], x[0,0]+a[0,0]], [x[0,1], x[0,1]+a[0,1]], 'b-', linewidth=3)
plt.legend(('NURBS Circle', 'Velocity', 'Acceleration'))

# update the velocity/acceleration lines for frame *i* in the animation
def animate(i):
    velocity.set_data(    [x[i,0], x[i,0]+v[i,0]], [x[i,1], x[i,1]+v[i,1]])
    acceleration.set_data([x[i,0], x[i,0]+a[i,0]], [x[i,1], x[i,1]+a[i,1]])
    plt.axis([-3, 3, -3, 3])
    
# create and show the animation
ani = animation.FuncAnimation(fig, animate, np.arange(1,n), interval=24)

if INTERACTIVE:
    plt.show()
else:
    # save results as an animated gif for web display (PS: this function call is slow)
    ani.save('circle.gif', writer='imagemagick', fps=30);
