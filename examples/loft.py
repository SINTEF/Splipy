# Lofting example.
#
# Author:    Kjetil Andre Johannessen
# Institute: SINTEF Digital, Mathematics & Cybernetics
# Date:      February 2021
#

# Slightly more interesting input geometries. Read more on NACA wing 
# profiles here: https://en.wikipedia.org/wiki/NACA_airfoil
from splipy.utils.NACA import NACA

# Standard splipy imports
from splipy import curve_factory
from splipy import surface_factory
from splipy.io.stl import STL

crv1 = curve_factory.circle(r=1) # base circle
crv2 = 1.2 * crv1 + [0,0,1]      # increase in size and offset in height
crv3 = 0.9 * crv1 + [0,0,2]
crv4 = 1.3*NACA(4,2,35) + [-.5,0,4] # sequence of NACA wing profiles
crv5 = 1.2*NACA(5,3,25) + [-.5,0,5]
crv6 = 0.8*NACA(8,3,15) + [-.5,0,7]
crv7 = 0.2*NACA(8,3,15) + [-.5,0,9]

wind_turbine_blade = surface_factory.loft(crv1, crv2, crv3, crv4, crv5, crv6, crv7)

# Dump result as an stl file which can be viewed in for instance Meshlab
with STL('blade.stl') as myfile:
    myfile.write(wind_turbine_blade.swap()) # swap() is to make sure normals are pointing out of the object
