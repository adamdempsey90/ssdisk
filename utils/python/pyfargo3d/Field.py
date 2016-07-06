from pylab import *
import re

class Parameters(object):
    """
    Class for reading the simulation parameters.
    """
    def __init__(self, paramfile):

        try:
            params = open(paramfile,'r')
        except IOError:
            print  paramfile + " not found."
            return

        lines = params.readlines()
        params.close()
        par = {}

        for line in lines:
            name, value = line.split()
            try:
                float(value)
            except ValueError:
                try:
                    int(value)
                except ValueError:
                    value = '"' + value + '"'
            par[name] = value

        self._params = par

        for name in par:
            exec("self."+name.lower()+"="+par[name])

class Domain(object):
    """
    Class for reading the domain of the simulation. The inputs are the
    name of the domain files.
    """
    def __init__(self, xfile=None, yfile=None, zfile=None):

        if xfile != None:
            try:
                domain_x = np.loadtxt(xfile)
                self.xmin = domain_x
                self.xmed = 0.5*(self.xmin[:-1]+self.xmin[1:])
            except IOError:
                print "IOError with ", xfile
                pass

        if yfile != None:
            try:
                domain_y = np.loadtxt(yfile)
                self.ymin = domain_y[3:-3]
                self.ymed = 0.5*(self.ymin[:-1]+self.ymin[1:])
            except IOError:
                print "IOError with ", yfile
                pass

        if zfile != None:
            try:
                domain_z = np.loadtxt(zfile)
                if domain_z.shape[0] > 3:
                    self.zmin = domain_z[3:-3]
                    self.zmed = 0.5*(self.zmin[:-1]+self.zmin[1:])
                else:
                    self.zmin = domain_z
                    self.zmed = 0.5*(self.zmin[:-1]+self.zmin[1:])
            except IOError:
                print "IOError with ", zfile
                pass

class Planet(object):
    """
    Planet class. It stores the position of the planet at a certain
    time. The time is given by the integer number n (output number), and the first occurrence is
    taken.
    """
    
    def __init__(self,planetfile,n):
        p = np.loadtxt(planetfile)
        index = np.where(p[:,0] == n)[0]
        index = index[0]
        xp = p[index,1];           yp  = p[index,2]
        
        r = np.sqrt(xp**2+yp**2);  phi = np.arctan2(yp,xp)

        self.x = float(xp);        self.y = float(yp)
        self.phi = float(phi);     self.r = float(r)
        self.z  = float(p[index,3])
        self.vx = float(p[index,4]) 
        self.vy = float(p[index,5]) 
        self.vz = float(p[index,6])
        self.m  = float(p[index,7])
        self.fileline = int(index)

class Field(object):
    def __init__(self, field, staggered = None, varname=None, domname=None, dtype='float64'):

        """
        A field class. It stores the scalar, parameters and the domain
        data, all in one place.  For this reason, always is needed the
        existence of the variables.par file. 

        staggered can be 'x', 'y', or 'yx' or 'xy', or 'xyz' and all
        the allowed permutations. If y is in the string, then the y
        direction will be filled with ymin, etc...
        """

        #Getting the parameters.
        if varname == None:
            var = "variables.par"
        else:
            var = varname
        directory = re.search("(.*)/",field).group(0)
        varfile = directory + var
        self.Params = Parameters(varfile)
        #Getting the domain
        if domname == None:
            dom = "domain"
        else:
            dom = domname
        domfile = directory + dom
        domain = Domain(domfile+"_x.dat",domfile+"_y.dat",domfile+"_z.dat")

        if staggered == None:
            self.x = domain.xmed
            self.y = domain.ymed
            self.z = domain.zmed
            self.xm = domain.xmin
            self.ym = domain.ymin
            self.zm = domain.zmin
            staggered = 'c'
            staggered = 'c'
        else:
            if staggered.count('x')>0:
                self.x = domain.xmin[:-1]
            else:
                self.x = domain.xmed
            if staggered.count('y')>0:
                self.y = domain.ymin[:-1]
            else:
                self.y = domain.ymed
            if staggered.count('z')>0:
                self.z = domain.zmin[:-1]
            else:
                self.z = domain.zmed

        self.staggered = staggered

        self.__dtype = dtype
        self.field = field
        self.data = self.__open_field()

        #Helpers
        self.xmin = self.x.min()
        self.xmax = self.x.max()
        self.ymin = self.y.min()
        self.ymax = self.y.max()
        self.zmin = self.z.min()
        self.zmax = self.z.max()
            
    def __open_field(self):
        """
        Reading and reshaping a field.
        """
        f = np.fromfile(self.field,dtype=self.__dtype)
        if self.Params.nz == 1:
            return f.reshape(self.Params.ny,self.Params.nx)
        else:
            return f.reshape(self.Params.nz,self.Params.ny,self.Params.nx)

    def plot(self,cmap='cubehelix',**karg):
        plt.imshow(self.data,interpolation='nearest',
                   origin='lower',cmap=cmap,
                   extent=[self.x.min(),self.x.max(),self.y.min(),self.y.max()],
                   aspect='auto',
                   **karg)

    def contour(self,**karg):
        plt.contour(self.data,interpolation='nearest',
                    origin='lower',
                    extent=[self.x.min(),self.x.max(),self.y.min(),self.y.max()],
                    aspect='auto',
                    **karg)

    def __bilinear(self,x,y,f,p):
        """
        Computing the bilinear interpolation.
        Parameters
        ----------
	
        x = (x1,x2); y = (y1,y2)
        f = (f11,f12,f21,f22)
        p = (x,y)
	
        where x,y is the interpolated point and
        fij is the value of the function at the
        point (xi,yj).
	
        Output
        ------
	
        f(p): Float. 
        The interpolated value of the function f(p) = f(x,y)
        """
	
        xp  = p[0]; yp   = p[1]; x1  = x[0]; x2  = x[1]
        y1  = y[0]; y2  = y[1];  f11 = f[0]; f12 = f[1]
        f21 = f[2]; f22 = f[3]
        t = (xp-x1)/(x2-x1);    u = (yp-y1)/(y2-y1)
	
        return (1.0-t)*(1.0-u)*f11 + t*(1.0-u)*f12 + t*u*f22 + u*(1-t)*f21
    
    def get_value(self, x, y):
        """
        For a real set of coordinates (x,y), returns the bilinear
        interpolated value of a Field class.
        """
	
        i = int((x-self.xmin)/(self.xmax-self.xmin)*self.Params.nx)
        j = int((y-self.ymin)/(self.ymax-self.ymin)*self.Params.ny)
	
        if i<0 or j<0 or i>self.data.shape[1]-2 or j>self.data.shape[0]-2:
            return None
	
        f11 = self.data[j,i]
        f12 = self.data[j,i+1]
        f21 = self.data[j+1,i]
        f22 = self.data[j+1,i+1]
        try:
            x1  = self.x[i]
            x2  = self.x[i+1]
            y1  = self.y[j]
            y2  = self.y[j+1]
            return self.__bilinear((x1,x2),(y1,y2),(f11,f12,f21,f22),(x,y))
        except IndexError:
            return None

def open_mscalar(monitor_name):
    """
    Reading a scalar monitor.
    It returns (t,monitor), where t was nomalized by 2pi.
    
    Use: t,torque = open_mscalar(torque.dat)
    """

    data = np.loadtxt(monitor_name)
    return data[:,0]/(2.0*np.pi),data[:,1]
