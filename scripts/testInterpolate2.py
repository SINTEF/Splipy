from math import *
from GoTools import *
from GeoUtils.Interpolate import *

n = 20
t = [float(i)/(n-1)*2*pi for i in range(n)]
y = []
dx = []
dy = []
for x in t:
	y.append(sin(x))
	dx.append(1)
	dy.append(cos(x))
crvs = []
crvs.append(CubicCurve(t,y, [], HERMITE, dx, dy))
crvs.append(CubicCurve(t,y, [], TANGENT, [1,1], [1,1]))
crvs.append(CubicCurve(t,y, [], NATURAL))
crvs.append(CubicCurve(t,y, [], PERIODIC))
crvs.append(CubicCurve(t,y, [], FREE))


# output file name
SetFinalOutput("curves.g2")
FinalOutput( crvs, True);
