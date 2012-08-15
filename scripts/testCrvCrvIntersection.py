from math import *
from GoTools import *
from GoTools.CurveFactory import *
from GoTools.SurfaceFactory import *
from GoTools.VolumeFactory import *
from GeoUtils.Interpolate import *

# standard nice things
origin = Point(0,0,0)
x_axis = Point(1,0,0)
y_axis = Point(0,1,0)
z_axis = Point(0,0,1)

N = 7;
M = 30;
x = [];
y = [];
z = [];

for i in range(N):
	theta = pi/4*i/(N-1) + 2*pi/2;
	x.append(  cos(theta) - 4 );
	y.append(  sin(theta) + 1 );
	z.append(      0          );

for i in range(M):
	t = 8.0* i/(M-1);
	x.append( -4 + t);
	y.append( 0 );
	z.append( 0 );

for i in range(N):
	theta = - pi/4*i/(N-1) + pi/2;
	x.append(  cos(theta) + 4 );
	y.append(  sin(theta) - 1 );
	z.append(      0          );

wigglyLine = Cubic(x,y,z);
shortLine  = LineSegment(origin, x_axis);

circ1 = Circle(x_axis, 1.5, z_axis);
circ2 = Circle(y_axis, 0.5, z_axis);

# intersect curve returns a Curve
oneLine, noPoints  = IntersectCurve(shortLine, wigglyLine);
# intersect curve returns 2 Points
noLine, twoPoints = IntersectCurve(circ1, circ2);

# have to make a line out of the two points since we can't store points
seg = LineSegment(twoPoints[0], twoPoints[1]);

FinalOutput([oneLine[0], seg], True)

