import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d

"""
For an efficiently developing, add these two lines in ipython shell:

In [1]: %load_ext autoreload
In [2]: %autoreload 2

"""

class Snapshot():
    def __init__(self, n, fields=['Density'], reader = 'Fargo3d', directory=None):
        """
        fields is a list with the name of the fields to be loaded.
        """
        for field in fields:
            """
            Here, I'm creating a list of variables with the name of the root fields.
            Note the interesting way to do that. I'm declaring dynamicaly the name of
            variables in a secondary interpreter, that shares the same namespace of
            the main program. The function exec is amazing for interactive sessions!
            """
#           exec('self.' + i + "=Reader.Fargo3d(input_file = (i+'{:06d}.dat').format(n))")
#           exec('self.' + field + "=Field((field+'{:06d}.dat').format(n))")
            exec('self.' + field + "=Field(n=n,name=field,directory=directory)")

        try:
            vx = self.Vx
            vy = self.Vy
        except AttributeError:
#            print "streams method is avalaible only for snapshots with Vx & Vy fields."
#            print "If you want to use it, please add Vx & Vy in fields. ej: Snapshot(n,fields=['Vx',Vy])"
            pass

class Streams():
    def __init__(self,vx,vy,planet):
        """
        vx & vy are Field classes.
        """

        self.Vx = vx
        self.Vy = vy
        self.xmin = vx.parameters['xmin']
        self.xmax = vx.parameters['xmax']
        self.ymin = vx.parameters['ymin']
        self.ymax = vx.parameters['ymax']
        self.domain_x = vx.domain['x'] #xmin(i)
        self.domain_y = vx.domain['y'] #ymin(j)
        self.nx   = vx.parameters['nx']
        self.ny   = vx.parameters['ny']
        self.planet = planet
        
    def bilinear(self,x,y,f,p):
        """
        x = (x1,x2); y = (y1,y2)
        f = (f11,f12,f21,f22)
        p = (x,y)
        where x,y is the interpolated point and
        fij is the value of the function at the 
        point (xi,yj).
        """
        xp  = p[0]; yp   = p[1]; x1  = x[0]; x2  = x[1]
        y1  = y[0]; y2  = y[1];  f11 = f[0]; f12 = f[1]
        f21 = f[2]; f22 = f[3]
        t = (xp-x1)/(x2-x1);    u = (yp-y1)/(y2-y1)

        return (1.0-t)*(1.0-u)*f11 + t*(1.0-u)*f12 + t*u*f22 + u*(1-t)*f21

    def get_vphi(self,x,y):
        i = int((x-self.xmin)/(self.xmax-self.xmin)*self.nx)
        j = int((y-self.ymin)/(self.ymax-self.ymin)*self.ny-0.5)
        f11 = self.Vx.data[j,i]
        f12 = self.Vx.data[j,i+1]
        f21 = self.Vx.data[j+1,i]
        f22 = self.Vx.data[j+1,i+1]
        x1 = self.domain_x[i]
        x2 = self.domain_x[i+1]
        y1 = 0.5*(self.domain_y[j]+self.domain_y[j+1])
        y2 = 0.5*(self.domain_y[j+1]+self.domain_y[j+2])
        vphi = self.bilinear((x1,x2),(y1,y2),
                             (f11,f12,f21,f22),
                             (x,y))
        return vphi

    def get_vrad(self,x,y):
        i = int((x-self.xmin)/(self.xmax-self.xmin)*self.nx+0.5)
        j = int((y-self.ymin)/(self.ymax-self.ymin)*self.ny)
        f11 = self.Vy.data[j,i]
        f12 = self.Vy.data[j,i+1]
        f21 = self.Vy.data[j+1,i]
        f22 = self.Vy.data[j+1,i+1]
        x1 = 0.5*(self.domain_x[i]+self.domain_x[i+1])
        x2 = 0.5*(self.domain_x[i+1]+self.domain_x[i+2])
        y1 = self.domain_y[j]
        y2 = self.domain_y[j+1]
        vrad = self.bilinear((x1,x2),(y1,y2),
                             (f11,f12,f21,f22),
                             (x,y))
        return vrad

    def __euler(self,x,y,frac=0.4,reverse=False):
        """
        Reverse inverts the sign of velocity
        """
        sign = 1.0
        if reverse:
            sign = -1
        vphi = self.get_vphi(x,y)
        vrad = self.get_vrad(x,y)
        l = ((self.xmax-self.xmin)/self.nx)**2 + ((self.ymax-self.ymin)/self.ny)**2
        h = np.sqrt(l/(vphi**2+vrad**2))
        
#        h = np.min(((self.xmax-self.xmin)/self.nx/np.abs(vphi/y),
#                    (self.ymax-self.ymin)/self.ny/np.abs(vrad)))
        h *= frac
        return sign*h*np.array([vphi/y,vrad])
        
    def __get_stream(self,x0,y0,reverse=False,frac=0.4,nmax=10**6,bidirectional=False):

        if(bidirectional):
            reverse = False
            s0 = self.__get_stream(x0,y0,reverse=False,bidirectional=False,nmax=nmax,frac=frac)
            print "s0 lista"
            s1 = self.__get_stream(x0,y0,reverse=True,bidirectional=False,nmax=nmax,frac=frac)
            print "s1 lista"
            return (s0,s1)

        print 'Computing streamline...'

        x = [];  y = []
        x.append(x0)
        y.append(y0)

        for i in xrange(nmax):
            ds = self.__euler(x0,y0,frac=frac,reverse=reverse)
            dx = ds[0];   dy = ds[1]
            if(np.sqrt(dx**2+dy**2)<1e-10):
                print "Warning: dt is very small, maybe you're in a stagnation point!" +\
                    "Please, select another initial point."
                break
            x0 += dx
            y0 += dy
            
            if (x0 > self.xmax-2*(self.xmax-self.xmin)/self.nx) or \
                    (y0 > self.ymax-2*(self.ymax-self.ymin)/self.ny) \
                    or x0<self.xmin or y0<self.ymin:
                print "Warning: It was reached the limit of the box, breaking..."
                break
            x.append(x0)
            y.append(y0)
        print "streamline was done in",i, "steps"
        stream = np.array([x,y])
        return stream

    def get_streams(self,x0,y1,y2,n=30,frac=0.5,nmax=10**6):
        values = np.arange(y1,y2,(y2-y1)/float(n))
        streams = []
        i = 0
        for y in values:
            s = self.__get_stream(x0,y,frac=frac,nmax=nmax)
            if(len(s[0])==1):
                print "Warning! Sreamline was null. Recomputing it in reverse mode."
                s = self.__get_stream(x0,y,frac=frac,nmax=nmax,reverse=True)
                if(len(s[0])==1):
                    print "The streamline",y,"have a problem, the lenght is Null... Skipped"
                    continue
            print "Streamline",i, "OK"
            streams.append(s)
            i += 1
        return streams


    def get_stagnation(self,w=0.1,tol=None):
        """
        Computes the stagnation points. Tol is the tolerance (in pixels)
        for the stagnation points.
        if local is true, the stagnation point is searched near to (xp,yp).
        w controls the width of the local area. w is the fraction of the size of the box that represents the local part. 
        """

        xmin = self.xmin
        xmax = self.xmax
        ymin = self.ymin
        ymax = self.ymax
        nx = self.nx
        ny = self.ny
        vx = self.Vx.data
        vy = self.Vy.data
        xp = self.planet[0]
        yp = self.planet[1]

        w = 0.5*w      #redefining w!!

        if tol == None:
            l = np.sqrt(((self.xmax-self.xmin)/self.nx)**2 +
                        ((self.ymax-self.ymin)/self.ny)**2)
            tol = 0.5*l # half width of a cell.

        lx = (xmax-xmin)
        ly = (ymax-ymin)
        ip = int((xp-xmin)/(xmax-xmin)*nx)
        jp = int((yp-ymin)/(ymax-ymin)*ny)

        wx = int(w*nx)
        wy = int(w*ny)
        x = self.domain_x
        y = self.domain_y + 0.5*ly/ny #Be careful with centering!!

        xx,yy = np.meshgrid(x,y)

        cx = plt.contour(xx[jp-wy:jp+wy,ip-wx:ip+wx],yy[jp-wy:jp+wy,ip-wx:ip+wx],vx[jp-wy:jp+wy,ip-wx:ip+wx],
                         levels=(0,)) #getting 0 contours

        x = self.domain_x + 0.5*lx/nx
        y = self.domain_y

        xx,yy = np.meshgrid(x,y)

        cy = plt.contour(xx[jp-wy:jp+wy,ip-wx:ip+wx],yy[jp-wy:jp+wy,ip-wx:ip+wx],vy[jp-wy:jp+wy,ip-wx:ip+wx],
                         levels=(0,))

        px = cx.collections[0].get_paths()
        py = cy.collections[0].get_paths()

        nx = np.shape(px)[0]
        ny = np.shape(py)[0]

        if(nx>1 or ny>1):
            print "x-contours=",nx,"y-contours=",ny
            s0,s1 = self.get_stagnation(w=w/2.0,tol=None)
            return s0,s1

        temp = []

        for i in range(nx):
            cx = px[i].vertices
#            plt.plot(cx[:,0],cx[:,1],'ko')
            for j in range(ny):
                cy = py[j].vertices
#                plt.plot(cy[:,0],cy[:,1],'ro')
                for k in cx:
                    dif = np.sqrt((k[0]-cy[:,0])**2 + (k[1]-cy[:,1])**2)
                    if np.any(dif < tol):
#                        index = np.where(dif<tol)
#                        print k[0],k[1]
#                        temp.append((cy[index,0], cy[index,1]))
                        temp.append((k[0],k[1]))
#                        return temp

        sx = []
        sy = []
        for l in temp:
#            sx.append(l[0][0][0])
#            sy.append(l[1][0][0])
            sx.append(l[0])
            sy.append(l[1])
            dif = np.sqrt((np.asarray(sx)-xp)**2+(np.asarray(sy)-yp)**2)
        try:
            index = np.where(dif == dif.min())[0]
            return sx[index],sy[index]
        except ValueError:
            print "Error getting the stagnation point! Try to increase 'tol' value."
            return 0,0

    def get_separatrix(self,niter=10,tol=None,noise=10.0,w=0.1,frac=0.05,x0=None,y0=None):
        """
        Noise is given in fraction of cells
        """

        separatrix = []
        dx = (self.xmax-self.xmin)/self.nx
        dy = (self.ymax-self.ymin)/self.ny
        
        if(x0 == None and y0==None):
            s0,s1 = self.get_stagnation(w=w,tol=tol)
            if(s0 == 0 and s1 == 0):
                return
            print "Stagnation point was founded:", s0,s1
        else:
            s0,s1 = x0,y0
        for i in range(niter):
            s = self.__get_stream(s0+noise*dx*(-0.5+np.random.rand()),
                                  s1+noise*dy*(-0.5+np.random.rand()),
                                  bidirectional=True,frac=frac)
            separatrix.append(s)
        return separatrix

    def plot_stagnation(self,s):
        plt.plot(s[0],s[1],'ro',ms=10)

    def plot_streams(self,streams):
        """
        Ploting the streams computed with get_streams method.
        """
        for stream in streams:
            plt.plot(stream[0],stream[1],'k')

    def plot_separatrix(self,separatrix):
        for s in separatrix:
            plt.plot(s[0][0],s[0][1],'r',linewidth=2.0)
            plt.plot(s[1][0],s[1][1],'r',linewidth=2.0)

    def get_map(self,alpha=0.1):

        ymin = self.planet[1]*(1-alpha)
        ymax = self.planet[1]*(1+alpha)
        streams1 = self.get_streams(-np.pi,ymin,ymax,n=30,frac=0.5)
        streams2 = self.get_streams(np.pi-2*(self.xmax-self.xmin)/self.nx,ymin,ymax,n=30,frac=0.5)
        self.plot_streams(streams1)
        self.plot_streams(streams2)
        separ = self.get_separatrix(niter=15,frac=0.1,noise=0.1)
        stag = self.get_stagnation()
        self.plot_stagnation(stag)
        self.plot_separatrix(separ)

        data = {}
        data['separatrix'] = separ
        data['stagnation'] = stag
        data['streams_left'] = streams1
        data['streams_right'] = streams2
        return data

class Reader():
    class __Fargo3d():
        def __init__(self, n, name=None,
                     dims_name="dims.dat",directory=None):
            if directory == None:
                input_file = name + "{0:06d}.dat".format(n)
            else:
                input_file = directory + name + "{0:06d}.dat".format(n)
                dims_name = directory + dims_name
            if name == None:
                return None
            if type(name) != str:
                print "Error! input_file must be a string"
                return
            self.directory = directory
            self.n = n
            self.name = input_file
            self.parameters = self.__shape(dims_name)
            self.axis   = self.__axis(self.parameters)
            self.domain = self.__domain(self.axis,self.parameters)
            self.data   = self.__reader(input_file,
                                      self.parameters,
                                      self.axis)
            self.planet = self.__planet(self.n)

        def __reader(self,input_file, param, axes):
            dim = []
            for i in axes:
                dim.append(int(param['n'+i]))
            return np.fromfile(input_file).reshape(dim)
            
        def __shape(self,dims_name):
            """
            Determines the shape of the simulation.
            """
            dims = open(dims_name,"r")
            lines = dims.readlines()
            name = lines[0].split("\t")
            vals = lines[1].split("\t")
            parameters = {}
            for i,j in zip(name,vals):
                parameters[i.lower()] = float(j)
            return parameters
        
        def __domain(self, axis,parameters):
            d = {}
            for i in axis:
                if i != 'x':
                    a = parameters['ngh'+i]
                    b = parameters['n'+i]
                else:
                    a = 0
                    b = -1
                if self.directory == None:
                    d[i] = np.loadtxt("domain_{0}.dat".format(i))[a:a+b]
                else:
                    d[i] = np.loadtxt(self.directory+"domain_{0}.dat".format(i))[a:a+b]
            return d
        
        def __axis(self,p):
            """
            Determines the dimensions of the simulation,
            ie: x, y, z, xy, yz, xz, xz, xyz
            """
            ax = ['z', 'y', 'x']
            if p['nx'] == 1: ax.remove('x')
            if p['ny'] == 1: ax.remove('y')
            if p['nz'] == 1: ax.remove('z')
            return ax

        def __planet(self,n):
            try:
		if self.directory == None:
                   p = np.loadtxt("planet0.dat")
            	else:
		   p = np.loadtxt(self.directory + "planet0.dat")
            except IOError: 
                return (0,0)

            np.where(p[:,0] == n)
            index = p[0][0]
            xp = p[index,1]; yp = p[index,2]
            r = np.sqrt(xp**2+yp**2)
            phi = np.arctan2(yp,xp)
            return (phi,r)

    def __init__(read_class,n,name=None,
                 dims_name="dims.dat",
                 reader='fargo3d',directory=None):
        """
        Note that here read_class is equivalent to self.
        """
        if(reader == 'fargo3d'):
            return read_class.__Fargo3d(n=n,name=name, dims_name=dims_name,directory=directory)

class Field(Reader):
    def __init__(self, n, name, reader='fargo3d',directory=None):
        reader = Reader.__init__(self, n=n, name=name, reader = 'fargo3d',directory=directory)
        self.name = reader.name
        self.n = n
        self.parameters = reader.parameters
        self.axis = reader.axis
        self.domain = reader.domain
        self.data = reader.data
        self.planet = reader.planet
        self.dx = (self.parameters['xmax'] - 
                   self.parameters['xmin'])/self.parameters['nx']
        self.dy = (self.parameters['ymax'] - 
                   self.parameters['ymin'])/self.parameters['ny']
        self.extent = [self.parameters['xmin'],
                       self.parameters['xmax'],
                       self.parameters['ymin'],
                       self.parameters['ymax']]
        self.nghy = self.parameters['nghy']
        self.nghz = self.parameters['nghz']

        return

    def Get_properties(self):
        print "Input file is:", self.name
        print "Number of snapshot si:", self.n
        print "The matrix of the data is:", self.data
        print "The parameters of the simulations are:", self.parameters
        print "The active axes are:", self.axis
        print "The domain of each axis is:", self.domain
        

    def __projection(self,x_old,y_old,data,x_new,y_new,zero=1e-8):
        """
        rmask = 0 --> Number of bad radial boundary cells (ideal for to do plots of
        zero  = value of the mask
        """
        from matplotlib.mlab import griddata

        data[0,:] = zero
        new_data = griddata(x_old.flat,
                        y_old.flat,
                        data.T.flat, # No olvidar trasponer!
                        x_new,
                        y_new,
                        interp='linear')
        return new_data
    
    def __pol2cart(self,r,t,data,frac=0.7, n=1000, multiperiodic=1):
        """
        frac=0.7 is the fraction of rmax that will be used for the box size.
        Box size is 2*frac*rmax
        n=1000 is the output resolution
        """
        rr, tt = np.meshgrid(r,t)
        xx = rr*np.cos(multiperiodic*tt)
        yy = rr*np.sin(multiperiodic*tt)
        rmax = frac*r.max()
        x = np.arange(-rmax,rmax,2.0*rmax/float(n))
        extent = [-rmax,rmax,-rmax,rmax]
        return  self.__projection(xx,yy,data,x,x),extent
        

    def Plot2d(self, log=False,cmap=plt.cm.hot, norm=None, aspect=None, interpolation=None,
               alpha=None, vmin=None, vmax=None, origin='lower', extent=None,
               shape=None, filternorm=1, filterrad=4.0, imlim=None, resample=None, 
               url=None, hold=None, rmask = 0, 
               projection='pol', n=1000, frac=0.7, multiperiodic=1.0):

        """
        For now, only works with 2d-fields
        Local parameters:
        
        log = False --> It activates logscale. (log10(abs(data)))
        rmask = 0 --> Number of bad radial boundary cells (ideal for to do plots of 
        simulations with strong radial boundary conditions.)
        projection=pol is the coordinate sistem for plotting. Options: {'pol','cart'}
        """

        if log == True:
            data = np.log10(np.abs(self.data[rmask:-(1+rmask),:]))
        else:
            data = self.data[rmask:-(1+rmask),:]

        if projection == 'cart':
            data,extent = self.__pol2cart(self.domain['y'][0+rmask:-(1+rmask)],self.domain['x'],data,
                                          frac=frac,n=n, multiperiodic=multiperiodic)            
            plt.imshow(data,
                       cmap=cmap, norm=norm, aspect=aspect, 
                       interpolation=interpolation, alpha=alpha, 
                       vmin=vmin, vmax=vmax, origin=origin,
                       shape=shape, filternorm=filternorm,
                       filterrad=filterrad, imlim=imlim, 
                       resample=resample, url=url, hold=hold,
                       extent=extent)

        if projection == 'pol':
            if extent == None:
                extent = [self.domain['x'][0],
                          self.domain['x'][-1],
                          self.domain['y'][0+rmask],
                          self.domain['y'][-(1+rmask)]]
                plt.xlabel(self.axis[1])
            plt.ylabel(self.axis[0])
            plt.imshow(data,
                       cmap=cmap, norm=norm, aspect=aspect,
                       interpolation=interpolation, alpha=alpha,
                       vmin=vmin, vmax=vmax, origin=origin,
                       shape=shape, filternorm=filternorm,
                       filterrad=filterrad, imlim=imlim,
                       resample=resample, url=url, hold=hold,
                       extent=extent)
        plt.show()
