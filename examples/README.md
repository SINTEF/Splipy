# Examples

This folder contains several examples aimed as a getting-started guide.


## Circle_animation.py

This is a simple code to generate a parametric NURBS circle. Note that it is
not arclength-parametrized and the acceleration vector is not continuous.
Disregarding the include-part of the code, as well as the plotting itself, the
code looks like this:

``` python

  n = 250                                  # number of evaluation points
  c = curve_factory.circle()               # create the NURBS circle (r=1)
  t = np.linspace(c.start(0), c.end(0), n) # parametric evaluation points
  x = c(t)                                 # physical (x,y)-coordinates, size (n,2)
  v = c.derivative(t, 1)                   # velocity at all points
  a = c.derivative(t, 2)                   # acceleration at all points
```

![Missing circle animation](http://i.imgur.com/8MaBiTW.gif "Circle animation")

## Read.py

Displays the interface to write geometries to file

``` python

  # Read multiple NURBS patches from the file 'teapot.g2'
  with G2('teapot.g2') as my_file:
      my_teapot = my_file.read()
```

## Write.py

Displays the interface to write geometries to file

``` python

  # create a NURBS torus
  torus = surface_factory.torus(minor_r=1, major_r=4)

  # G2 files are native GoTools (http://www.sintef.no/projectweb/geometry-toolkits/gotools/)
  with G2('torus.g2') as my_file:
      my_file.write(torus)
```

## Lissajous.py
Lissajous curves. A family of parametric curves of the type

```
x = A sin(at+d) 
y = B sin(bt)
```

More info: [https://en.wikipedia.org/wiki/Lissajous_curve](https://en.wikipedia.org/wiki/Lissajous_curve). Again, stripping all inclusion, and animation parts of the code, one can generate these curves in the following way


``` python 

def lissajous(a, b, d):
  # request a,b integers, so we have closed, periodic curves
  n = gcd(a,b)
  N = (a/n) * (b/n) # number of periods before looping

  # compute a set of interpolation points
  numb_pts = max(3*N, 100) # using 3N interpolation points is decent enough
  t = np.linspace(0,2*pi/n, numb_pts)
  x = np.array([np.sin(a*t + d), np.sin(b*t)])

  # do a cubic curve interpolation with periodic boundary conditions
  return curve_factory.cubic_curve(x.T, curve_factory.Boundary.PERIODIC)
```

![Missing Lissajous curve animation](http://i.imgur.com/HKr59BT.gif "lissajous(3,4,pi/2)")

Animation of the lissajous curve with a=3, b=4 and d=pi/2

![Missing Lissajous curve animation](http://i.imgur.com/6q7aAUM.gif "lissajous(60,44,pi/2)")

Animation of the lissajous curve with a=60, b=44 and d=pi/2
