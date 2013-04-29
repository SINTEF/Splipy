__doc__ = 'Implementation of various interpolation schemes.'

import GoTools 
import GoTools.CurveFactory 

def Linear(x,y=[],z=[]):
        """Linear interpolate a list of points (arclength parametrization) in an approximative sense
        @param x: The x-coordinate of the points to interpolate
        @type x: List of floats
        @param y: y-coordinates
        @type y: List of floats
        @param z: z-coordinates
        @type z: List of floats
        @return: Linear interpolated curve in 3 dimensions
        @rtype: Curve
        """
        if(len(y) == 0):
                for i in range(len(x)):
                        y.append(0)
        if(len(z) == 0):
                for i in range(len(x)):
                        z.append(0)
        if(len(x) != len(y) or len(x) != len(y) or len(y) != len(z)):
                print "Error: coordinate lists needs to be of same length"
                return None 

        N = len(x)
        pts = []
        t   = [0]
        for i in range(N):
                pts.append(GoTools.Point(x[i],y[i],z[i]))
                if(i>0):
                        t.append(t[i-1] + abs(pts[i]-pts[i-1])) # eucledian distance between points i and (i-1)
        return GoTools.CurveFactory.ApproximateCurve(pts, t, 2)

def LinearP(pts):
        """Linear interpolate a list of points (arclength parametrization) in an approximative sense
        @param pts: The coordinate of the points to interpolate
        @type pts: List of Points
        @return: Linear interpolated curve
        @rtype: Curve
        """
        N = len(pts)
        t   = [0]
        for i in range(1,N):
                t.append(t[i-1] + abs(pts[i]-pts[i-1])) # eucledian distance between points i and (i-1)
        return GoTools.CurveFactory.ApproximateCurve(pts, t, 2)

def Cubic(x, y=[], z=[]):
        """Cubic spline interpolation a list of points by arclength parametrization  in an approximative sense
        @param x: The x-coordinate of the points to interpolate
        @type x: List of floats
        @param y: y-coordinates
        @type y: List of floats
        @param z: z-coordinates
        @type z: List of floats
        @return: Cubic interpolated curve in 3 dimensions
        @rtype: Curve
        """
        if(len(y) == 0):
                for i in range(len(x)):
                        y.append(0)
        if(len(z) == 0):
                for i in range(len(x)):
                        z.append(0)
        if(len(x) != len(y) or len(x) != len(y) or len(y) != len(z)):
                print "Error: coordinate lists needs to be of same length"
                return None 

        N   = len(x)
        pts = []
        t   = [0]
        for i in range(N):
                pts.append(GoTools.Point(x[i],y[i],z[i]))
                if(i>0):
                        t.append(t[i-1] + abs(pts[i]-pts[i-1])) # eucledian distance between points i and (i-1)
        return GoTools.CurveFactory.ApproximateCurve(pts, t, 4)

def CubicP(pts):
        """Cubic spline interpolation a list of points by arclength parametrization in an approximative sense
        @param pts: The coordinate of the points to interpolate
        @type pts: List of Points
        @return: Cubic interpolated curve 
        @rtype: Curve
        """
        N = len(pts)
        t   = [0]
        for i in range(1,N):
                t.append(t[i-1] + abs(pts[i]-pts[i-1])) # eucledian distance between points i and (i-1)
        return GoTools.CurveFactory.ApproximateCurve(pts, t, 4)
        
def UniformCubic(x, y=[], z=[]):
        """Cubic spline interpolation a list of points by uniform parametrization
        @param x: The x-coordinate of the points to interpolate
        @type x: List of floats
        @param y: y-coordinates
        @type y: List of floats
        @param z: z-coordinates
        @type z: List of floats
        @return: Cubic interpolated curve in 3 dimensions
        @rtype: Curve
        """
        if(len(y) == 0):
                for i in range(len(x)):
                        y.append(0)
        if(len(z) == 0):
                for i in range(len(x)):
                        z.append(0)
        if(len(x) != len(y) or len(x) != len(y) or len(y) != len(z)):
                print "Error: coordinate lists needs to be of same length"
                return None 

        N   = len(x)
        pts = []
        for i in range(N):
                pts.append(GoTools.Point(x[i],y[i],z[i]))
        return GoTools.CurveFactory.ApproximateCurve(pts, range(N), 4)
        
def UniformCubicP(pts):
        """Cubic spline interpolation a list of points by uniform parametrization
        @param pts: The coordinate of the points to interpolate
        @type pts: List of Points
        @return: Cubic interpolated curve 
        @rtype: Curve
        """
        return GoTools.CurveFactory.ApproximateCurve(pts, range(len(pts)), 4)

