from math import *
from GoTools import *
from GoTools.CurveFactory import *

# standard nice things
origin = Point(0,0,0)
x_axis = Point(1,0,0)
y_axis = Point(0,1,0)
z_axis = Point(0,0,1)

SetFinalOutput('circleSegments.g2');

# test different normal axis
c1 = CircleSegment(origin, x_axis, pi/2, z_axis);
c2 = CircleSegment(origin, y_axis, pi/2, x_axis);
c3 = CircleSegment(origin, x_axis, -pi/2, y_axis);

result = [c1,c2,c3];

# test non-90 degree angles
N  = 10;
dt = 2.0*pi/N;
for i in range(N):
	t = i*dt;
	x = Point(cos(t), sin(t), 0);
	result.append(CircleSegment(origin, 2*x, dt, z_axis));

# test large circles above 120 degrees
result.append(CircleSegment(origin, [3*cos(pi/10), 3*sin(pi/10), 0], 5*pi/6, z_axis)); # two knotspans
result.append(CircleSegment(origin, [4*cos(pi/8),  4*sin(pi/8),  0], 4*pi/3, z_axis)); # two knotspans
result.append(CircleSegment(origin, [5*cos(pi/6),  5*sin(pi/6),  0], 5*pi/3, z_axis)); # 3   knotspans
result.append(CircleSegment(origin, [6*cos(pi/4),  6*sin(pi/4),  0], 6*pi/3, z_axis)); # 3   knotspans

# test non-origo circles
result.append( CircleSegment(4*z_axis, [4,0,4], pi, z_axis) );
result.append( CircleSegment([4,0,4], [4,3,4], pi, [1, 0, 1]) );

FinalOutput(result)
