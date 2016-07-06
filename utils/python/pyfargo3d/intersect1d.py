import numpy as np
from scipy.interpolate import spline,interp1d
from scipy.optimize import bisect

class Curve():
    """
    Curve class.
    Input x,y arrays of the same lenght
    """
    def __init__(self,x,y):
        assert len(x) == len(y), \
            "x,y must to have the same lenght"
        self.x = x
        self.y = y

class Intersect():
    """
    Class for computing intersections between curves...

    Parameters:
    -----------
    "n": integer
         is the number of neighbors
    "rseach": integer 
              in the number of oversampling in the common domain
              for building F(x)=(f2-f1)(x) --> F(x0)=0
              If zeros the curves have a strong curvature, try to increase this 
              value to 5000 or 10000. The performance of the algorithm is governed
              by rsearch
    """
    def __init__(self,Curve1,Curve2,n=1,rsearch=2000):
        self.Curve1 = Curve1
        self.Curve2 = Curve2
        self.__f1 = interp1d(Curve1.x,Curve1.y)
        self.__f2 = interp1d(Curve2.x,Curve2.y)
        self.x1max = np.max(self.Curve1.x)
        self.x2max = np.max(self.Curve2.x)
        self.x1min = np.min(self.Curve1.x)
        self.x2min = np.min(self.Curve2.x)
        self.x,self.f  = self.build_f(rsearch)
        if self.x != None:
            self.xinit = self.scan_domain(self.x,self.f,n=n)
        else:
            self.xinit = None

#        self.commons = self.__check_common_points()
#        print self.commons

#    def __check_common_points(self):
#        x1 = self.Curve1[0]
#        x2 = self.Curve2[0]
#        y1 = self.Curve1[1]
#        y2 = self.Curve2[1]
#        common = np.where(x1 == x2)[0]
#        dif = y1[common]-y2[common]
#        zero = np.where(dif == 0)[0]
#        if len(common[zero])>1:
#            print "hola"
#        print x1[common[zero]],y1[common[zero]]

    def intersect(self):
        intersections = []
        if self.xinit != None:
            for ab in self.xinit:
                if(self.f(ab[0])*self.f(ab[1])<0):
                    intersection = bisect(self.f,ab[0],ab[1])
                    intersections.append((intersection,
                                          self.__f1(intersection).tolist()))
            if len(intersections)==0:
                intersections = None
            return intersections
        else:
            return None

    def build_f(self,rsearch):
        x = self.common_domain(rsearch)
        if x != None:
            y1 = self.__f1(x)
            y2 = self.__f2(x)
            f = interp1d(x,y2-y1)
            return x,f
        else:
            return None, None
        
    def common_domain(self,rsearch):

        #Testing if the common domain is empty:
        if self.x1min > self.x2max or self.x2min>self.x1max:
            return None 

        if self.x1max > self.x2max:
            xmax = self.x2max
        else:
            xmax = self.x1max
        if self.x1min < self.x2min:
            xmin = self.x2min
        else:
            xmin = self.x1min
        return np.linspace(xmin,xmax,rsearch)

    def scan_domain(self,x,f,n=5):
        """
        n is the security offset in points for to define the initial 
        point for bisection.
        """
        sign = np.sign(f(x[:-1])*f(x[1:]))
        indices = np.where(sign == -1)[0]
        xinit = []
        for i in indices:
            #xinit.append((x[i-n],x[i+n])) #Old version, tested and working...
            xinit.append((x[i],x[i+1])) #plus  one for problems in the sign when the change is in the first cell... Be careful with it. Experimental!!!
        if len(xinit) == 0:
            xinit = None
        return xinit

def intersect1d(x1,y1,x2,y2,n=1):
    """
    (C) Pablo Benitez Llambay 2013

    Find the crossing points between two functions.
    
    Parameters:
    ----------
    
    xi : array_like
         A 1-D array of "monotonically increasing" real domain values.
    yi : array_like
         A 1-D array of real image values.
    n  : integer value,
         offset for the zero-search. If the algorithm fails, try to change it
         
    Output:
    ------
    list of intersection points ([x1,y1],[x2,y2],...,[xn,yn]). If there are
    not common points or the algorithm fails, None is returned

    Example:
    ---------
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np

    >>> x1 = np.linspace(-3*np.pi,3*np.pi,100)
    >>> y1 = np.sin(x1)

    >>> x2 = np.linspace(-np.pi,np.pi,50)
    >>> y2 = np.cos(3*x2)

    >>> intersections = intersect1d(x1,y1,x2,y2)

    >>> print intersections
    [(-2.7490875470774356, -0.3820514515610326), (-1.1697127956833637, -0.9168432341132113), (-0.7913747419692514, -0.7084720265124655), (0.3911822841670822, 0.37960439194115453), (1.9676368248483191, 0.9199381857557692), (2.351982808803209, 0.7083950296504822)]

    for plotting:
    -------------

    plt.plot(x1,y1,'k')
    plt.plot(x2,y2,'r')
    for point in intersections:
      plt.plot(point[0],point[1],'go',markersize=10)
    plt.show()

    Note: Your curve must to be x-monotonically increasing!

    """

    curve1 = Curve(x1,y1)
    curve2 = Curve(x2,y2)
    
    intersection = Intersect(curve1,curve2,n=n).intersect()

    return intersection

if __name__ == "__main__":
    import doctest
    doctest.testmod()

    import matplotlib.pyplot as plt
    import numpy as np

    x1 = np.linspace(-3*np.pi,3*np.pi,100)
    y1 = np.sin(x1)
    x2 = np.linspace(-np.pi,np.pi,50)
    y2 = np.cos(3*x2)

    intersections = intersect1d(x1,y1,x2,y2)

    plt.plot(x1,y1,'k')
    plt.plot(x2,y2,'r')
    for point in intersections:
      plt.plot(point[0],point[1],'go',markersize=10)
    plt.show()
