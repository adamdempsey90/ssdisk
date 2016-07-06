import numpy as np
import matplotlib.pyplot as plt
import copy
import re
import os

NAVIGATION=1
EVALUATION=2

def common(F1, F2):
    
    Field1 = copy.deepcopy(F1)
    Field2 = copy.deepcopy(F2)

    if Field1.stag.count('x') !=  Field2.stag.count('x'):
        data = Field2.get_data()
        if Field2.stag.count('x')>0:
            Field2.set_stag(Field2.stag.replace('x',''))
            Field2.xn -= 1
        else:
            Field2.set_stag(Field2.stag.replace('c','')+'x')
            Field2.x0 += 1
        Field2.set_data(0.5*(data[:,1:]+data[:,:-1]))

    if Field1.stag.count('y') !=  Field2.stag.count('y'):
        data = Field2.get_data()
        if Field2.stag.count('y')>0:
            Field2.set_stag(Field2.stag.replace('y',''))
            Field2.yn -= 1
        else:
            Field2.set_stag(Field2.stag.replace('c','')+'y')
            Field2.y0 += 1
        Field2.set_data(0.5*(data[1:,:]+data[:-1,:]))

    if (Field1.x0>Field2.x0):
        Field2.x0 = Field1.x0
    else:
        Field1.x0 = Field2.x0
        
    if (Field1.y0>Field2.y0):
        Field2.y0 = Field1.y0
    else:
        Field1.y0 = Field2.y0

    if (Field1.xn>Field2.xn):
        Field1.xn = Field2.xn
    else:
        Field2.xn = Field1.xn

    if(Field1.yn>Field2.yn):
        Field1.yn = Field2.yn
    else:
        Field2.yn = Field1.yn

    return Field1, Field2

def dx(F):
    """
    Diferenciate the field along the x direction (phi). This operation return 
    a wrong value at the last x-cells...
    """

    Dx = copy.deepcopy(F)

    X,Y = np.meshgrid(F.get_x(),F.get_y())
    data = (F.get_data())
    dx = (data[:,1:] - data[:,:-1])/(X[:,1:]-X[:,:-1])
    
    if F.stag.count('x')>0:
        Dx.xn -= 1
        Dx.set_stag(F.stag.replace('x',''))
    else:
        Dx.x0 += 1
        Dx.set_stag(F.stag.replace('c','') + 'x')

    Dx.set_data(dx)
    return Dx

def dy(F):
    """
    Diferenciate the field along the x direction (phi). This operation return 
    a wrong value at the last x-cells...
    """

    Dy = copy.deepcopy(F)

    X,Y = np.meshgrid(F.get_x(),F.get_y())
    data = (F.get_data())
    dy = (data[1:,:] - data[:-1,:])/(Y[1:,:]-Y[:-1,:])
    
    if F.stag.count('y')>0:
        Dy.yn -= 1
        Dy.set_stag(F.stag.replace('y',''))
    else:
        Dy.y0 += 1
        Dy.set_stag(F.stag.replace('c','') + 'y')

    Dy.set_data(dy)

    return Dy

def div(vx, vy):
    T,R = np.meshgrid(vy.get_x(), vy.get_y())
    rvr=vy*R
    t1 = f.dy(rvr)
    t2 = f.dx(vx)
    d = t1 + t2
    T,R = np.meshgrid(d.get_x(), d.get_y())
    return d/R

class Parameters(object):
    """
    Class for reading the simulation parameters.  

    Input: paramfile: the paramters file, normally called
    variables.par 

    Output: The Prameter Object is filled with the parameters readed.
    The type of each one is infered from a cast-error procedure.
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
    name of the domain files. By default is assumed that the names are
    domain_i.dat

    Input: xfile, yfile, zfile -> domain files.  By default is assumed
    that y-z domain files have 3 ghost cells.
    """

    def __init__(self, xfile=None, yfile=None, zfile=None, ngh=3):

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
                if domain_y.shape[0] > ngh:
                    self.ymin = domain_y[ngh:-ngh]
                    self.ymed = 0.5*(self.ymin[:-1]+self.ymin[1:])
                else:
                    self.ymin = domain_y
                    self.ymed = 0.5*(self.ymin[:-1]+self.ymin[1:])
            except IOError:
                print "IOError with ", yfile
                pass

        if zfile != None:
            try:
                domain_z = np.loadtxt(zfile)
                if domain_z.shape[0] > ngh:
                    self.zmin = domain_z[ngh:-ngh]
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
    output number. The first occurrence of n in the planet file is
    taken.
    """
    
    def __init__(self, planetfile, n):

        p = np.loadtxt(planetfile)
        index = np.where(p[:,0] == n)[0]
        index = index[0]
        xp = p[index,1]
        yp  = p[index,2]
        r = np.sqrt(xp**2+yp**2)
        phi = np.arctan2(yp,xp)
        self.x = float(xp)
        self.y = float(yp)
        self.phi = float(phi)
        self.r = float(r)
        self.z  = float(p[index,3])
        self.vx = float(p[index,4])
        self.vy = float(p[index,5])
        self.vz = float(p[index,6])
        self.m  = float(p[index,7])
        self.fileline = int(index)

class SnapShot(object):
    def __init__(self, directory, n):
        
        """
        The main class of pyfargo3d.py. It stores all the data of a
        simulation, including the magnetic fields.
        
        Parameters:

        Input: directory (string): where the simulation was stored.
        Output: A SnapShot class.

        """

        if directory[-1] != '/':
            directory += '/'

        print "Reading " + directory

        self.planet = Planet(directory+"planet0.dat",n)

        domfile = directory + "domain"
        self.Domain = Domain(domfile+"_x.dat",domfile+"_y.dat",domfile+"_z.dat")
        self.Params = Parameters(directory+"variables.par")

        targets = {"gasdens":"c", "gasenergy":"c",
                   "gasvx":"x", "gasvy":"y", "gasvz":"z",
                   "bx":"x", "by":"y", "bz":"z"}
        
        files = os.listdir(directory)

        for target in targets:
            full_target =  target+"{0:d}.dat".format(n)
            if full_target in files:
                exec("self." + target + "=Field('" + 
                     directory + target + "{0:d}.dat'".format(n) + 
                     ",'"+targets[target]+"',self.Params, self.Domain)")
                print full_target+" read..."

    def new_field(self, fieldname, staggered='c'):
        exec("self." + fieldname + "= Field(None,'"+ staggered +"',self.Params, self.Domain)")

    def __explore(self,event,key,stype,color,nticks):

        """
        Allowed stypes are 'v' -> velocity,'b' -> magnetic, 'd' d-field.
        """

        if event.key.lower() == key:
            nav = plt.figure(NAVIGATION)
            p = self.get_stream(event.xdata,event.ydata,stype=stype)
            plt.plot(p[0],p[1],color,linewidth=1)
            ntot = len(p[0])
            if event.key == key.upper():
                for i in range(1,nticks):
                    plt.text(p[0][ntot/nticks*i],p[1][ntot/nticks*i], str(i), color=color)
                ev = plt.figure(EVALUATION)
                d,e = self.evaluate(p[0], p[1], self.navfield, stype=stype)
                ntot = len(d)
                plt.plot(d,e,color)
                for i in range(1,nticks):
                    plt.text(d[ntot/nticks*i],e[ntot/nticks*i], str(i), color=color)

    def __onpress(self, event):

        nticks = 10
        
        self.__explore(event,'v','v','r',nticks)
        self.__explore(event,'b','b','b',nticks)
        if event.key.lower() == 'd':
            try:
                self.dx
            except: 
                self.compute_d()
            self.__explore(event,'d','d','g',nticks)
        
        if event.key.lower() == 'q':
            try:
                self.__delete(NAVIGATION, nticks)
            except:
                pass
            try:
                self.__delete(EVALUATION, nticks)
            except:
                pass

        print event.xdata,"," ,event.ydata

    def __delete(self,fig,nticks):

        plt.figure(fig)
        ax = plt.gca()
        ax.lines.pop()
        for i in range(nticks-1):
            ax.texts.pop()
        plt.draw()
                
    def navigate(self, Field, **karg):
        fig = plt.figure(NAVIGATION)
        self.navfield = Field
        Field.plot(**karg)
        cid1 = fig.canvas.mpl_connect('key_press_event', self.__onpress)

    def get_stream(self, x0, y0, nmax=100000, lmax=None, stype='v'):
        import stream as s

        if lmax != None:
            lmax_f = lmax
        else:
            lmax_f = 2*np.pi*y0

        if stype == 'v':
            v1 = self.gasvx.get_harddata()
            x1 = self.gasvx.get_x()
            y1 = self.gasvx.get_y()
            v2 = self.gasvy.get_harddata()
            x2 = self.gasvy.get_x()
            y2 = self.gasvy.get_y()
        elif stype == 'b':
            v1 = self.bx.get_harddata()
            x1 = self.bx.get_x()
            y1 = self.bx.get_y()
            v2 = self.by.get_harddata()
            x2 = self.by.get_x()
            y2 = self.by.get_y()
        elif stype == 'd':
            v1 = self.dx.get_harddata()
            x1 = self.dx.get_x()
            y1 = self.dx.get_y()
            v2 = self.dy.get_harddata()
            x2 = self.dy.get_x()
            y2 = self.dy.get_y()
        else:
            print "The allowed options are: v, b & d."
            return None

        stream1 = s.stream(v1, x1, y1, v2, x2, y2, x0, y0, nmax, lmax_f, 1)
        stream2 = s.stream(v1, x1, y1, v2, x2, y2, x0, y0, nmax, lmax_f, -1)
        
        s1 = range(0,len(stream1[0]),1)
        s2 = range(0,-len(stream2[0]),-1)

        ind = np.argsort(s1+s2)

        stream_x = np.array(stream1[0]+stream2[0])[ind]
        stream_y = np.array(stream1[1]+stream2[1])[ind]

        if lmax != None:
            x = stream_y*np.cos(stream_x)
            y = stream_y*np.sin(stream_x)
            dx = x[1:]-x[:-1]
            dy = y[1:]-y[:-1]
            l = np.cumsum(np.sqrt(dx**2+dy**2))
            i = np.where(l<lmax)
            return [stream_x[i],stream_y[i]]
        else:
            return [stream_x,stream_y]

    def get_portrait(self, n=100, portype='v'):
        """
        portype is teh type of portrait (v,b or d).
        """
        portrait = []
        for i in range(n):
            x0 = -3.0+6*np.random.rand()
            y0 = 0.9+0.2*i/n
            portrait.append(self.get_stream(x0, y0, sign=1, stype=portype))
            portrait.append(self.get_stream(x0, y0, sign=-1, stype=portype))
        return portrait

    def __get_stagnation(self, field1, field2, x0, y0, w=0.1, plot=False, **karg):

        from get_contour import get_contour
        from intersect_curves import intersect_curves

        """
        Getting the zero velocity point for a staggered mesh
        Don't forget to use a corotating reference frame.
        
        Parameters: 
        ----------- 
        
        v1, v2:  Field class.
        x0, y0:  float
                 The center of the searching box.
        w = 0.1: float.  
                 Is the width of the target zone. Is
                 in units of the width of the main box
        plot = False
                 If it is true, you will see the zero velocity curves. For debuggin only.
                 **karg are the arguments supported by the matplotlib plot function.

        Output: 
        ------- 
        A list with all the stagnation points around
        the planet and contours [(s0x,s0y),...,(snx,sny)].

        """

        nx = self.Params.nx;    ny = self.Params.ny

        x1min = field1.xmin;    x1max = field1.xmax
        x2min = field2.xmin;    x2max = field2.xmax
        
        y1min = field1.ymin;    y1max = field1.ymax
        y2min = field2.ymin;    y2max = field2.ymax
        
        i1 = int((x0-x1min)/(x1max-x1min)*nx);  j1 = int((y0-y1min)/(y1max-y1min)*ny)
        i2 = int((x0-x2min)/(x1max-x1min)*nx);  j2 = int((y0-y2min)/(y2max-y2min)*ny)
        
        w *= 0.5
        
        w1xp  = int(i1+w*nx);    w1xm  = int(i1-w*nx)
        w1yp  = int(j1+w*ny);    w1ym  = int(j1-w*ny)
        
        w2xp  = int(i2+w*nx);    w2xm  = int(i2-w*nx)
        w2yp  = int(j2+w*ny);    w2ym  = int(j2-w*ny)
        
        X1, Y1 = np.meshgrid(field1.x,field1.y)
        contour1 = get_contour(X1[w1ym:w1yp,w1xm:w1xp],
                               Y1[w1ym:w1yp,w1xm:w1xp],
                               field1.data[w1ym:w1yp,w1xm:w1xp],0.0)
        
        X2, Y2 = np.meshgrid(field2.x,field2.y)
        contour2 = get_contour(X2[w2ym:w2yp,w2xm:w2xp],
                               Y2[w2ym:w2yp,w2xm:w2xp],
                               field2.data[w2ym:w2yp,w2xm:w2xp],0.0)
    
        if contour1 == None or contour2 == None:
            print "I can't find zero velocities. Try to play with w parameter..."
            return None

        if np.shape(contour1)[0] > 1 or np.shape(contour2)[0] > 1:
            print "Error: I found more than one contour of zero velocity," \
                " set a small w and try to avoid this problem."
            if plot:
                for c in contour1:
                    plt.plot(c[0],c[1],'k')
                for c in contour2:
                    plt.plot(c[0],c[1],'r')
            return None

        x1 = np.array(contour1[0][0])
        y1 = np.array(contour1[0][1])
        x2 = np.array(contour2[0][0])
        y2 = np.array(contour2[0][1])
        
        stagnations = intersect_curves(x1,y1,x2,y2,1)
    
        if plot:
            for c1 in contour1:
                plt.plot(c1[0],c1[1],**karg)
                for c2 in contour2:
                    plt.plot(c2[0],c2[1],**karg)
            for s in stagnations:
                plt.plot(s[0],s[1],'ko',**karg)
            plt.show()

        return stagnations#,contour1,contour1]

    def get_stagnation(self, x0, y0, w=0.1, plot=False, stype='v', **karg):
        if stype == 'v':
            v1 = self.gasvx
            v2 = self.gasvy
        elif stype == 'b':
            v1 = self.bx
            v2 = self.by
        elif stype == 'd':
            v1 = self.dx
            v2 = self.dy
        else:
            print "The allowed options are: v, b & d."
            return None
        return self.__get_stagnation(v1, v2, x0, y0, w=0.1, plot=False, **karg)

    def get_separatrix(self, x0=None, y0=None, n=10, noise=0.1, w=0.1, stype='v'):
        """
        Noise is given in fraction of cells
        """

        separatrix = []

        dx = (self.Params.xmax-self.Params.xmin)/self.Params.nx
        dy = (self.Params.ymax-self.Params.ymin)/self.Params.ny
        
        if x0 == None: 
            x0 = 0;
        if y0 == None:
            y0 = 1.0;
            
        s0,s1 = self.get_stagnation(x0,y0,w=w,stype=stype)[0]
        if(s0 == 0 and s1 == 0):
            return
#        print "Stagnation point was found:", s0,s1
        for i in range(n):
            x0 = s0+noise*dx*(-0.5+np.random.rand())
            y0 = s1+noise*dy*(-0.5+np.random.rand())
            s = self.get_stream(x0, y0, stype=stype)
            separatrix.append(s)
        return separatrix

    def __curl(self, Field1, Field2):
        X1,Y1 = np.meshgrid(Field1.get_x(),Field1.get_y())
        T1 = dy(Field1*Y1)
        T2 = dx(Field2)
        curl = T1-T2
        X,Y = np.meshgrid(curl.get_x(), curl.get_y())
        curl = curl/Y
        return curl

    def compute_vorticity(self):
        self.vorticity = self.__curl(self.gasvx,self.gasvy) + 2*self.Params.omegaframe

    def compute_current(self):
        self.current = self.__curl(self.bx,self.by)

    def compute_vortencity(self):
        self.check_field("vorticity")
        self.vortencity = self.vorticity/self.gasdens
        
    def compute_v2(self):
        self.v2 = self.gasvx**2+self.gasvy**2

    def compute_bernoulli(self, xp=None, yp=None, pmass=2e-5):

        """
        This routine computes the bernoulli's invariant, defined as:
        B = 0.5*v**2 + eta + phi + , where eta = cs**2*log(rho).
        xp, yp are the cartesian coordinates of of the planet.
        pmass is the mass of the planet.
        """

        if xp == None:
            xp = self.planet.x
        if yp == None:
            yp = self.planet.y

        self.check_field("v2")

        t1 = self.v2*0.5
        logdens = copy.deepcopy(self.gasdens)
        logdens.data = np.log(logdens.data)
        t2 = self.gasenergy**2*logdens

        temp = copy.deepcopy(self.gasdens)

        T,R = np.meshgrid(temp.get_x(), temp.get_y())
        X = R*np.cos(T)
        Y = R*np.sin(T)

        rp = np.sqrt(xp**2+yp**2)
        h  = self.Params.aspectratio
        f  = self.Params.flaringindex
        tm = self.Params.thicknesssmoothing
        omega = self.Params.omegaframe

        temp.set_data(-1.0/np.sqrt(X**2+Y**2) - pmass/np.sqrt((X-xp)**2+(Y-yp)**2+(h*rp**f*rp*tm)**2))
        temp.data -= 0.5*omega**2*R**2
        self.bernoulli = t1+t2+temp

    def compute_d(self):
        
        """
        Computing the D-Field.
        D = v - |\alpha|/alpha*j/(rho*w)*B
        donde alpha = J/omega
        """

        self.check_field("vorticity")
        self.check_field("current")

        alpha = self.current/self.vorticity
        
        factor = self.current/(self.gasdens*self.vorticity)

        self.dx = self.gasvx - factor*self.bx
        self.dy = self.gasvy - factor*self.by

    def compute_angular(self):
        """
        The angular momentum is given in the inertial frame
        """
        
        X,Y = np.meshgrid(self.gasvx.get_x(),self.gasvx.get_y())
        vx = self.gasvx.data + Y*self.Params.omegaframe
        self.angular = self.gasvx*Y
        
    def check_field(self, field):
        try:
            exec("self."+field)
        except AttributeError:
            print "Computing "+field+"..."
            exec("self.compute_"+field+"()")

    def evaluate(self, x, y, Field, stype='v'):
        """
        Make a cut along the curve x. The result is an array d-y,
        where y is the value of the field, and d is the distance to
        the stagnation point (of d-field or v-field).
        """

        ev = [];   d  = []
        dx = 0
        dy = 0
        l = 0
        xold = y[0]*np.cos(x[0])
        yold = y[0]*np.sin(x[0])
        for x0,y0 in zip(x,y):
            ev.append(Field.get_value(x0,y0))
            dx = y0*np.cos(x0) - xold
            dy = y0*np.sin(x0) - yold
            l += np.sqrt(dx**2+dy**2)
            d.append(l)
            xold = y0*np.cos(x0)
            yold = y0*np.sin(x0)
        return np.array(d),np.array(ev)

class Field():

    def __neg__(F1):
        F = copy.deepcopy(F1)
        F.set_data(-F1.get_data())
        return F

    def __abs__(F1):
        F = copy.deepcopy(F1)
        F.set_data(np.abs(F1.get_data()))
        return F

    def __add__(F1,F2):
        if type(F2) == type(F1):
            Field1, Field2 = common(F1,F2)
            F = copy.deepcopy(Field1)
            F.set_data(Field1.get_data() + Field2.get_data())
        else:
            F = copy.deepcopy(F1)
            F.set_data(F1.get_data() + F2)
        F.config_field()
        return F

    def __sub__(F1,F2):
        if type(F2) == type(F1):
            Field1, Field2 = common(F1,F2)
            F = copy.deepcopy(Field1)
            F.set_data(Field1.get_data() - Field2.get_data())
        else:
            F = copy.deepcopy(F1)
            F.set_data(F1.get_data() - F2)
        F.config_field()
        return F

    def __mul__(F1,F2):
        if type(F2) == type(F1):
            Field1, Field2 = common(F1,F2)
            F = copy.deepcopy(Field1)
            F.set_data(Field1.get_data() * Field2.get_data())
        else:
            F = copy.deepcopy(F1)
            F.set_data(F1.get_data() * F2)
        F.config_field()
        return F

    def __div__(F1,F2):
        if type(F2) == type(F1):
            Field1, Field2 = common(F1,F2)
            F = copy.deepcopy(Field1)
            F.set_data(Field1.get_data() /  Field2.get_data())
        else:
            F = copy.deepcopy(F1)
            F.set_data(F1.get_data() / F2)
        F.config_field()
        return F

    def __pow__(F1, n):
        F = copy.deepcopy(F1)
        F.set_data(F1.get_data()**2)
        return F

    def merge_domain(self):

        NX = len(self.__xmed)
        NY = len(self.__ymed)
        self.__x = np.ndarray(2*NX+1)
        for i in range(NX):
            self.__x[2*i] = self.__xmin[i]
            self.__x[2*i+1] = self.__xmed[i]
        self.__x[2*NX] = self.__xmin[NX]

        self.__y = np.ndarray(2*NY+1)
        for i in range(NY):
            self.__y[2*i] = self.__ymin[i]
            self.__y[2*i+1] = self.__ymed[i]
        self.__y[2*NY] = self.__ymin[NY]

    def get_x(self):
        return copy.deepcopy(self.__x[2*(self.x0)+self.__xs:2*self.xn:2])

    def get_y(self):
        return copy.deepcopy(self.__y[2*(self.y0)+self.__ys:2*self.yn:2])

    def get_extent(self):
        return copy.deepcopy([self.x0,self.xn,self.y0,self.yn])

    def set_data(self,data):
        self.data[self.y0:self.yn,self.x0:self.xn] = data

    def get_harddata(self):
        return copy.deepcopy(self.data[self.y0:self.yn,self.x0:self.xn])

    def get_data(self):
        return self.data[self.y0:self.yn,self.x0:self.xn]

    def set_stag(self,stag):
        if stag.count('x')>0:
            self.__xs = 0
        else:
            self.__xs = 1
        if stag.count('y')>0:
            self.__ys = 0
        else:
            self.__ys = 1
        self.stag = stag

    def get_stag(self,stag):
        return self.stag

    def config_field(self):
        self.xmin = self.get_x()[0]
        self.xmax = self.get_x()[-1]
        self.ymin = self.get_y()[0]
        self.ymax = self.get_y()[-1]
        self.x = self.get_x()
        self.y = self.get_y()
        return

    def __init__(self, field, stag, Params, Domain):

        domain = copy.deepcopy(Domain)
        
        self.set_stag(stag)

        self.__xmin = domain.xmin
        self.__xmed = 0.5*(self.__xmin[1:]+self.__xmin[:-1])
        self.__ymin = domain.ymin
        self.__ymed = 0.5*(self.__ymin[1:]+self.__ymin[:-1])

        self.x0 = 0
        self.xn = Params.nx
        self.y0 = 0
        self.yn = Params.ny

        self.merge_domain()

        self._Params = Params
        self._Domain = Domain

        self.stag = stag
        self.dtype = Params.realtype
        self.field = field

        self.data = self.__open_field(Params, field)

        self.config_field()

    def __open_field(self, Params, field):
        """
        Reading and reshaping a field.
        """
        nx = Params.nx
        ny = Params.ny
        nz = Params.nz
        if field != None:
            f = np.fromfile(self.field,dtype=Params.realtype)
        else:
            f = np.ndarray(nx*ny*nz)
        shape = [nz,ny,nx]
        try:
            shape.remove(1)
        except ValueError:
            pass
        return f.reshape(shape)

    def __get_indices(self,x0,x1,y0,y1):

        extent = self.get_extent()

        nx = extent[1]-extent[0]
        ny = extent[3]-extent[2]

        i0 = int((x0-self.xmin)/(self.xmax-self.xmin)*nx)
        j0 = int((y0-self.ymin)/(self.ymax-self.ymin)*ny)
        i1 = int((x1-self.xmin)/(self.xmax-self.xmin)*nx)
        j1 = int((y1-self.ymin)/(self.ymax-self.ymin)*ny)

        if i0 < 0 or i0 > nx:
            i0 = 0
        if j0 < 0 or j0 > ny:
            j0 = 0
        if i1 > nx or i1 < 0:
            i1 = nx-1
        if j1 > ny or j1 < 0:
            j1 = ny-1
        
        return [i0,i1,j0,j1]

    def plot(self, cmap='cubehelix', extent=None, log=False, aspect='auto', **karg):

        if extent == None:
            extent = [self.xmin,self.xmax,self.ymin,self.ymax]
            ext = self.get_extent()
            nx = ext[1]-ext[0]
            ny = ext[3]-ext[2]
            j0 = 0; j1 = ny
            i0 = 0; i1 = nx
        else:            
            i0,i1,j0,j1 = self.__get_indices(*extent)
            extent = [self.get_x()[i0],self.get_x()[i1],self.get_y()[j0],self.get_y()[j1]]
        
        if log:
            data = np.log10(np.abs(self.data))
        else:
            data = self.get_data()
        
        plt.imshow(data[j0:j1,i0:i1],#interpolation='nearest',
                   origin='lower',cmap=cmap,
                   extent=extent,
                   aspect=aspect,
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

        data = self.get_data()
        ny = data.shape[0]
        nx = data.shape[1]


        i = int((x-self.xmin)/(self.xmax-self.xmin)*nx)
        j = int((y-self.ymin)/(self.ymax-self.ymin)*ny)
	
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

    def compute_gradient(self):
        self.gy = dy(self)
        dx_temp = dx(self)
        X,Y = np.meshgrid(dx_temp.get_x(),dx_temp.get_y())
        self.gx = dx_temp/Y
