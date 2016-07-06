from matplotlib import _cntr as cntr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from subprocess import call
from scipy.optimize import curve_fit
from scipy.integrate import cumtrapz
import copy
import glob

class Mesh():
    """
    Mesh class, for keeping all the mesh data.
    Input: directory [string] -> place where the domain files are.
    """
    def __init__(self, directory=""):
        if len(directory) > 1:
            if directory[-1] != '/':
                directory += '/'
        try:
            domain_x = np.loadtxt(directory+"domain_x.dat")
        except IOError:
            print "IOError with domain_x.dat"
        try:
            #We have to avoid ghost cells!
            domain_y = np.loadtxt(directory+"domain_y.dat")[3:-3]
        except IOError:
            print "IOError with domain_y.dat"
        self.xm = domain_x #X-Edge
        self.ym = domain_y #Y-Edge

        self.xmed = 0.5*(self.xm[:-1] + self.xm[1:]) #X-Center
        self.ymed = 0.5*(self.ym[:-1] + self.ym[1:]) #Y-Center

        self.dx = diff(self.xm)
        self.dy = diff(self.ym)
        self.surfy = outer(self.ym[:-1],self.dx)
        self.surfx = outer(self.dy,ones(self.dx.shape))

        self.vol = outer(.5*(self.ym[1:]**2-self.ym[:-1]**2),self.dx)

        self.nx = len(self.xmed)
        self.ny = len(self.ymed)-2

        #(Surfaces taken from the edges)
        #First we make 2D arrays for x & y, that are (theta,r)
#        T,R = meshgrid(self.xm, self.ym)
#        R2  = R*R
#        self.surf = 0.5*(T[:-1,1:]-T[:-1,:-1])*(R2[1:,:-1]-R2[:-1,:-1])


class Parameters():
    """
    Class for reading the simulation parameters.
    input: string -> name of the parfile, normally variables.par
    """
    def __init__(self, directory=''):
        if len(directory) > 1:
            if directory[-1] != '/':
                directory += '/'
        try:
            params = open(directory+"variables.par",'r') #Opening the parfile
        except IOError:                  # Error checker.
            print  paramfile + " not found."
            return
        lines = params.readlines()     # Reading the parfile
        params.close()                 # Closing the parfile
        par = {}                       # Allocating a dictionary
        for line in lines:             #Iterating over the parfile
            name, value = line.split() #Spliting the name and the value (first blank)
            try:
                float(value)           # First trying with float
            except ValueError:         # If it is not float
                try:
                    int(value)         #                   we try with integer
                except ValueError:     # If it is not integer, we know it is string
                    value = '"' + value + '"'
            par[name] = value          # Filling the dictory
        self._params = par             # A control atribute, actually not used, good for debbuging
        for name in par:               # Iterating over the dictionary
            exec("self."+name.lower()+"="+par[name]) #Making the atributes at runtime
        try:
            self.alpha
        except AttributeError:
            self.alpha = self.alphaviscosity


class Field(Mesh, Parameters):
    """
    Field class, it stores all the mesh, parameters and scalar data
    for a scalar field.
    Input: field [string] -> filename of the field
           staggered='c' [string] -> staggered direction of the field.
                                      Possible values: 'x', 'y', 'xy', 'yx'
           directory='' [string] -> where filename is
           dtype='float64' (numpy dtype) -> 'float64', 'float32',
                                             depends if FARGO_OPT+=-DFLOAT is activated
    """
    Q_dict = {'dens':'gasdens{0:d}.dat','vx':'gasvx{0:d}.dat','vy':'gasvy{0:d}.dat'}
    Q_dict['mdot'] = 'temp_files/mdot.dat'
    Q_dict['fd'] = 'temp_files/fd.dat'
    Q_dict['fw'] = 'temp_files/fw.dat'
    Q_dict['lambda_dep'] = 'temp_files/lambda_dep.dat'
    Q_dict['lambda_ex'] = 'temp_files/lambda_ex.dat'
    Q_dict['rhostar'] = 'temp_files/rhostar.dat'
    Q_dict['taurr'] = 'temp_files/tensor.dat'
    Q_dict['taupp'] = 'temp_files/tensor.dat'
    Q_dict['taurp'] = 'temp_files/tensor.dat'
    Q_dict['lstar'] = 'temp_files/lstar.dat'
    Q_dict['Ld'] = 'temp_files/Ld.dat'
    Q_dict['Lw'] = 'temp_files/Lw.dat'
    Q_dict['pot'] = 'temp_files/pot.dat'
    Q_dict['rhoslope'] = 'temp_files/rhoslope.dat'

    name_dict={'dens':'$\\Sigma$', 'vx': '$v_\\phi$', 'vy': '$v_r$'}
    name_dict['mdot'] = '$\\dot{M}$'
    name_dict['fd'] = '$F_{d}$'
    name_dict['fw'] = '$F_{w}$'
    name_dict['lambda_dep'] = '$\\Lambda_{dep}$'
    name_dict['lambda_ex'] = '$\\Lambda_{ex}$'
    name_dict['rhostar'] = '$\\Sigma^*$'
    name_dict['taurr'] = '$\\Pi_{rr}$'
    name_dict['taupp'] = '$\\Pi_{\\phi\\phi}$'
    name_dict['taurp'] = '$\\Pi_{r\\phi}$'
    name_dict['lstar'] = '$l^*$'
    name_dict['Ld'] = '$L_d$'
    name_dict['Lw'] = '$L_w$'
    name_dict['pot'] = '$\\partial_\\phi \\Phi$'

    name_dict['rhoslope'] = '$\\partial_r \\Sigma^*$'

    def __init__(self, Q, num, staggered='c', directory='', dtype='float64'):
        if len(directory) > 1:
            if directory[-1] != '/':
                directory += '/'
        Mesh.__init__(self, directory) #All the Mesh attributes inside Field!
        Parameters.__init__(self, directory) #All the Parameters attributes inside Field!
        self.staggered = staggered
        #Now, the staggering:
        if staggered.count('x')>0:
            self.x = self.xm[:-1] #Do not dump the last element
        else:
            self.x = self.xmed
        if staggered.count('y')>0:
            self.y = self.ym[:-1]
        else:
            self.y = self.ymed

        try:
            field = self.Q_dict[Q].format(num)
        except KeyError:
            print '%s not a valid Field'%Q
            print 'Please choose one of:, ',self.Q_dict.keys()
            raise
        self.xx,self.yy = meshgrid(self.x,self.y)
        self.name = Q
        self.math_name = self.name_dict[Q]
        if 'tau' in self.name.lower():
            tens=self.name
        else:
            tens=None
        print 'Loading ',directory+field
        self.data = self.__open_field(directory+field,dtype,tens=tens) #The scalar data is here.
      #  self.data = self.set_boundary()
        try:
            self.alpha
        except AttributeError:
            self.alpha = self.alphaviscosity
        self.recalculate()
        try:
            self.wkz = self.ymin + (self.ymax-self.ymin)*.0497
            self.wkzr = self.ymax - (self.ymax-self.ymin)*.19
        except AttributeError:
            self.wkzr = self.ymax
            self.wkz = self.ymin
        #self.ft = fft.rfft(self.data,axis=1)/self.nx
        #self.avg = self.data.mean(axis=1)


    def __open_field(self, f, dtype,tens=None):
        """
        Reading the data
        """
        try:
            field = fromfile(f, dtype=dtype)
        except IOError:
            print "Couldn't find %s, trying again with un-merged file"%f
            try:
                field = fromfile(f.replace('.dat','_0.dat'),dtype=dtype)
            except IOError:
                raise
        if tens is None:
            return field.reshape(self.ny, self.nx)
        else:
            if tens == 'taurr':
                return field[:self.nx*self.ny].reshape(self.ny, self.nx)
            elif tens == 'taupp':
                return field[self.nx*self.ny:2*self.nx*self.ny].reshape(self.ny, self.nx)
            elif tens == 'taurp':
                return field[-self.nx*self.ny:].reshape(self.ny, self.nx)
            else:
                print 'Invalid tensor specified'
                return None


    def recalculate(self):
        self.avg = self.data.mean(axis=1)
        self.ft = fft.rfft(self.data,axis=1)/self.nx
        self.grad = self.gradient()
        self.power = (self.ft*conj(self.ft)*self.dy[:,newaxis]).sum(axis=0).real

    def gradient(self):
        q = copy.copy(self.data)
        res = zeros(q.shape)
        one_dim = False
        try:
            q.shape[1]
        except IndexError:
            one_dim = True


        for i in range(1,q.shape[0]-1):
            if one_dim:
                res[i] = (q[i+1]-q[i-1])/(self.y[i+1]-self.y[i-1])
            else:
                res[i,:] = (q[i+1,:]-q[i-1,:])/(self.y[i+1]-self.y[i-1])

        if one_dim:
            res[0] = (q[1]-q[0])/(self.y[1]-self.y[0])
            res[-1] = (q[-1]-q[-2])/(self.y[-1]-self.y[-2])
        else:
            res[0,:] = (q[1,:]-q[0,:])/(self.y[1]-self.y[0])
            res[-1,:] = (q[-1,:]-q[-2,:])/(self.y[-1]-self.y[-2])
        return res


    def set_boundary(self,accretion=True):
        if self.name == 'vx':
            iact = 1; igh = 0;
            inner_ring = (self.data[iact,:] + self.ymed[iact]*self.omegaframe)*sqrt(self.ymed[iact]/self.ymed[igh]) - self.ymed[igh]*self.omegaframe
            iact = -2; igh = -1;
            outer_ring = (self.data[iact,:] + self.ymed[iact]*self.omegaframe)*sqrt(self.ymed[iact]/self.ymed[igh]) - self.ymed[igh]*self.omegaframe
        if self.name == 'vy':
            if accretion:
                iact=1; igh=0;
                inner_ring = -1.5*self.alpha*self.aspectratio*self.aspectratio*self.ymed[igh]**(2*self.flaringindex-.5)
                inner_ring *= ones(self.data[igh,:].shape)
                iact=-2; igh=-1;
                outer_ring = self.data[iact,:]
            else:
                inner_ring = self.data[iact,:]
                outer_ring = self.data[iact,:]

        if self.name == 'dens':
            if accretion:
                iact=1;igh=0;
                inner_ring = self.data[iact,:]* (self.ymed[iact]/self.ymed[igh])**(2*self.flaringindex+.5)
                iact=-2;igh=-1;
                nu_0 = self.alpha*self.aspectratio*self.aspectratio
                fac = 3*pi*nu_0* (self.ymed[iact])**(2*self.flaringindex+1)
                outer_ring = (self.data[iact,:]*fac  + self.mdot*(sqrt(self.ymed[igh])-sqrt(self.ymed[iact])))/(3*pi*nu_0*self.ymed[igh]**(2*self.flaringindex+1))
            else:
                iact=1;igh=0;
                inner_ring = self.data[iact,:]
                iact=-2;igh=-1;
                outer_ring = self.data[iact,:]

        return vstack( (inner_ring, self.data, outer_ring) )

    def center_to_edge(self,direction='y'):
        newfield = copy.deepcopy(self)
        newfield.staggered = ''
        if direction.count('y') > 0:
            newfield.data = vstack( (newfield.data, newfield.data[-1,:]) )
            newfield.data = .5*(newfield.data[1:,:]+  newfield.data[:-1,:])
            newfield.staggered += 'y'
            newfield.y = newfield.ym[:-1]
        if direction.count('x') > 0:
            newfield.data = hstack( ( newfield.data, newfield.data[:,0].reshape(len(newfield.data[:,0]),1)) )
            newfield.data = .5*(newfield.data[:,1:]+  newfield.data[:,:-1])
            newfield.x = newfield.xm[:-1]
            newfield.staggered += 'x'
        newfield.recalculate()
        return newfield
    def edge_to_center(self):
        newfield = copy.deepcopy(self)
        if newfield.staggered.count('y') > 0:
            newfield.data = vstack( (newfield.data, newfield.data[-1,:]) )
            newfield.data = .5*(newfield.data[1:,:] + newfield.data[:-1,:])
            newfield.staggered.strip('y')
        if newfield.staggered.count('x') > 0:
            newfield.data = hstack( (newfield.data[:,-1].reshape(len(newfield.data[:,-1]),1), newfield.data, newfield.data[:,0].reshape(len(newfield.data[:,0]),1)) )
            newfield.data = .5*(newfield.data[:,1:] + newfield.data[:,:-1])
            newfield.staggered.strip('x')
        newfield.avg = newfield.data.mean(axis=1)
        return newfield


    def plotmode_summary(self,mag=False,planet=None,xlims=None,ylims=None,log=False,**karg):


        fig,axes=subplots(3,3,sharex=True)
        subplots_adjust(hspace=0)

        for m,ax in zip(range(9),axes.flatten()):
            self.plotmode(m,ax=ax,planet=planet,xlims=xlims,ylims=ylims,mag=mag,logy=log,**karg)

        for ax in axes.flatten():
            ax.set_title('')
            ax.set_ylabel('')
            ax.set_xlabel('')

        for ax in axes.flatten()[-3:]:
            if planet:
                ax.set_xlabel('$(r-a)/H(a)$',fontsize=20)
            else:
                ax.set_xlabel('$r$',fontsize=20)

        if mag:
            axes[1,0].set_ylabel('$|$'+self.math_name+'$|^2$',fontsize=20,rotation=0,labelpad=20)
        else:
            axes[1,0].set_ylabel(self.math_name,fontsize=20,rotation=0,labelpad=20)
        return fig,axes

    def plotmode(self,m,ax=None,norm=1,shift=0,mag=False,planet=None,xlims=None,ylims=None,logx=False,logy=False,**karg):

        try:
            m[0]
        except TypeError:
            m= array([m])


        if planet is None:
            y = copy.copy(self.y)
            xstr = '$r$'
        else:
            y = (self.y-planet)/(self.aspectratio*planet)
            xstr = '$(r-a)/H$'

        fontsize=karg.pop('fontsize',20)
        figsize = karg.pop('figsize',(10,8))
        if ax is None:
            fig=figure(figsize=figsize)
            ax=fig.add_subplot(111)

        color = karg.pop('color','k')

        if mag:
            for i in m:
                ax.plot(y,(self.ft[:,i]*conj(self.ft[:,i])).real,'-',label='$m=%d$'%i,**karg)
            ax.set_ylabel('$|\\hat{' + self.math_name.strip('$') + '}|^2$' ,fontsize=fontsize)

        else:
            for i in m:
                ax.plot(y,self.ft[:,i].real,'-',label='$m=%d$'%i,**karg)
            ax.set_prop_cycle(None)
            for i in m:
                ax.plot(y,self.ft[:,i].imag,'--',**karg)
            ax.set_ylabel('$\\hat{' + self.math_name.strip('$') + '}$' ,fontsize=fontsize)



        ax.set_xlabel(xstr,fontsize=fontsize)

        if len(m) < 6:
            ax.set_title('$m=$' + ','.join(['%d'%i for i in m]),fontsize=fontsize)
            ax.legend(loc='best')
        else:
            ax.set_title('$m = %d...%d$'%(min(m),max(m)),fontsize=fontsize)

        if logx and planet is None:
            ax.set_xscale('log')

        if logy and mag:
            ax.set_yscale('log')
        ax.minorticks_on()
        if xlims is not None:
            ax.set_xlim(xlims)
        if ylims is not None:
            ax.set_ylim(ylims)

    def plotpower(self,ax=None,log=False,logx=False,logy=False,xlims=None,ylims=None,norm=False,**karg):
        if ax is None:
            fig=figure()
            ax=fig.add_subplot(111)

        dat  = copy.copy(self.power)
        if norm and dat[0] != 0:
            dat /= dat[0]

        ax.plot(dat,**karg)
        ax.set_xlabel('$m$',fontsize=20)
        ax.set_ylabel('$|$' + self.math_name + '$|_m^2$',fontsize=20)
        if logx or log:
            ax.set_xscale('log')
        if logy or log:
            ax.set_yscale('log')

        if xlims is not None:
            ax.set_xlim(xlims)
        if ylims is not None:
            ax.set_ylim(ylims)



    def plotavg(self,ax=None,log=False,logx=False,logy=False,xlims=None,ylims=None,planet=None,norm=1,shift=0,**karg):
        fontsize=karg.pop('fontsize',20)
        figsize = karg.pop('figsize',(10,8))


        if ax is None:
            fig=figure(figsize=figsize)
            ax=fig.add_subplot(111)
        y = copy.copy(self.y)
#        if self.staggered.count('y')>0:
#            y = copy.copy(self.ym[:-1])
#        else:
#            y = copy.copy(self.ymed)

        xstr = '$r$'
        if planet is not None:
            y = (y - planet)/(self.aspectratio*planet)
            xstr = '$ (r-a)/H $'
            logx=False
            if log:
                logy=True
                log=False

        ax.plot(y,(self.avg-shift)/norm,**karg)
        if logx or log:
            ax.set_xscale('log'),
        if logy or log:
            ax.set_yscale('log')
        ax.set_xlabel(xstr,fontsize=fontsize)
        ax.set_ylabel(self.math_name,fontsize=fontsize)
        ax.minorticks_on()
        if xlims is not None:
            ax.set_xlim(xlims)
        if ylims is not None:
            ax.set_ylim(ylims)



    def plot(self, norm=1,shift=0,ax=None,clrbar=True,log=False,logx=False, abslog=False,cartesian=False, title=None,cmap='viridis',sample=1, rasterized=False,**karg):
        """
        A layer to plt.imshow or pcolormesh function.
        if cartesian = True, pcolormesh is launched.
        """

        if self.x[0] == pi or self.x[-1] != pi:
            dp = diff(self.x)[0]
            dp *= 2
            yn =self.data[:,-1] + (pi - self.x[-1])*(self.data[:,0]-self.data[:,-1])/dp
            data = (copy.copy(self.data)-shift)/norm
            data[:,-1] = yn
            data[:,0] = yn
            x = copy.copy(self.x)
            x[-1] = pi
            x[0] = -pi

        else:
            data = (copy.copy(self.data)-shift)/norm
            x =copy.copy(self.x)

        ny,nx = data.shape
        indy = range(0,ny,sample)
        indx = range(10) + range(10,nx-10,sample)+range(nx-10,nx)
        data = data[indy,:][:,indx]
        x = x[indx]
        y = self.y[indy]
        shading = karg.pop('shading','gouraud')
        interpolation = karg.pop('interpolation','bilinear')
        fontsize = karg.pop('fontsize',20)
        xlims = karg.pop('xlims',None)
        ylims = karg.pop('ylims',None)
        T,R = meshgrid(x,y)
        if ax == None:
            fig=figure()
            ax=fig.add_subplot(111)
        if log:
            data = np.log10(data)
        if abslog:
            data = np.log10(np.abs(data))
        if cartesian:
            X = R*cos(T)
            Y = R*sin(T)
            line2d=ax.pcolormesh(X,Y,data,cmap=cmap,shading=shading,**karg)
            if xlims is not None:
                ax.set_xlim(xlims)
            if ylims is not None:
                ax.set_ylim(ylims)

        else:
            if logx:
                line2d=ax.pcolormesh(log10(R),T,data, cmap = cmap,shading=shading,**karg)
                ax.set_xlabel('$\\log_{10}(r)$',fontsize=fontsize)
                if xlims is not None:
                    ax.set_xlim(log10(xlims[0]),log10(xlims[1]))
                else:
                    ax.set_xlim(log10(y[0]),log10(y[-1]))
            else:
                line2d=ax.pcolormesh(R,T,data, cmap = cmap,shading=shading,**karg)
                ax.set_xlabel('$r$',fontsize=fontsize)
                if xlims is not None:
                    ax.set_xlim(xlims)
                else:
                    ax.set_xlim(y[0],y[-1])

            #origin='lower',aspect='auto',interpolation=interpolation,extent=[x[0],x[-1],y[0],y[-1]],**karg)
            ax.set_ylabel('$\\phi$',fontsize=fontsize)
            if ylims is not None:
                ax.set_ylim(ylims)
            else:
                ax.set_ylim(self.x[0],self.x[-1])
        if clrbar:
            cbar = colorbar(line2d,ax=ax)

        if title is None:
            if log:
                ax.set_title('$\\log_{10}$'+self.math_name,fontsize=fontsize)
            elif abslog:
                ax.set_title('$\\log_{10}|$'+self.math_name+'$|$',fontsize=fontsize)
            else:
                ax.set_title(self.math_name,fontsize=fontsize)
        else:
            ax.set_title(title,fontsize=fontsize)

    def contour(self, log=False, cartesian=False, **karg):
        if log:
            data = np.log(self.data)
        else:
            data = self.data
        ax = gca()
        T,R = meshgrid(self.x,self.y)
        if cartesian:
            X = R*cos(T)
            Y = R*sin(T)
            contour(X,Y,data,**karg)
        else:
            contour(T,R,data,**karg)


    def shift_field(Field,direction):
        import copy
        """
        Half cell shifting along the direction provided by direction
        direction can be ('x','y', 'xy', 'yx').

        After a call of this function, Field.xm/xmed has not
        sense anymore (it is not hard to improve).
        """
        F = copy.deepcopy(Field)
        if direction.count('x')>0:
            F.data = 0.5*(Field.data[:,1:]+Field.data[:,:-1])
            F.x = 0.5*(Field.x[1:]+Field.x[:-1])
        if direction.count('y')>0:
            F.data = 0.5*(F.data[1:,:]+F.data[:-1,:])
            F.y = 0.5*(F.y[1:]+F.y[:-1])

        F.nx = len(F.x)
        F.ny = len(F.y)

        return F


    def cut_field(Field, direction, side):
        """
        Cutting a field:
        Input: field --> a Field class
               axis  --> 'x', 'y' or 'xy'
               side  --> 'p' (plus), 'm' (minnus), 'pm' (plus/minnus)
        """

        cutted_field = copy.deepcopy(Field)
        ny,nx = Field.ny, Field.nx
        mx = my = px = py = 0

        if direction.count('x')>0:
            if side.count('m')>0:
                mx = 1
            if side.count('p')>0:
                px = 1
        if direction.count('y')>0:
            if side.count('m')>0:
                my = 1
            if side.count('p')>0:
                py = 1

        cutted_field.data = Field.data[my:ny-py,mx:nx-px]
        cutted_field.x = cutted_field.x[mx:nx-px]
        cutted_field.y = cutted_field.y[my:ny-py]

        return cutted_field

    def vector_field(vx,vy, **karg):
        nsx = nsy = 3
        T,R = meshgrid(vx.x[::nsx],vx.y[::nsy])
        X = R*cos(T)
        Y = R*sin(T)
        vx = vx.data[::nsy,::nsx]
        vy = vy.data[::nsy,::nsx]
        U = vy*cos(T) - vx*sin(T)
        V = vy*sin(T) + vx*cos(T)
        ax = gca()
        ax.quiver(X,Y,U,V,scale=5,pivot='midle', **karg)

    def euler(vx, vy, x, y, reverse):
        """
        Euler integrator for computing the streamlines.
        Parameters:
        ----------

        x,y: Float.
             starter coordinates.
        reverse: Boolean.
                 If reverse is true, the integratin step is negative.
        Reverse inverts the sign of the velocity

        Output
        ------

        (dx,dy): (float,float).
                 Are the azimutal and radial increment.
                 Only works for cylindrical coordinates.
        """
        sign = 1.0
        if reverse:
            sign = -1
        vphi = get_v(vx,x,y)
        vrad = get_v(vy,x,y)
        if vphi == None or vrad == None: #Avoiding problems...
            return None,None
        l = np.min((((vx.xmax-vx.xmin)/vx.nx),((vx.ymax-vx.ymin)/vx.ny)))
        h = 0.5*l/np.sqrt((vphi**2+vrad**2))

#        return sign*h*np.array([vphi/y,vrad])
        return sign*h*np.array([vphi,vrad])


    def bilinear(x,y,f,p):
        """
        Computing bilinear interpolation.
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

    def get_v(v, x, y):
        """
        For a real set of coordinates (x,y), returns the bilinear
        interpolated value of a Field class.
        """

        i = int((x-v.xmin)/(v.xmax-v.xmin)*v.nx)
        j = int((y-log(v.ymin))/(log(v.ymax/v.ymin))*v.ny)

        if i<0 or j<0 or i>v.data.shape[1]-2 or j>v.data.shape[0]-2:
            return None

        f11 = v.data[j,i]
        f12 = v.data[j,i+1]
        f21 = v.data[j+1,i]
        f22 = v.data[j+1,i+1]
        try:
            x1  = v.x[i]
            x2  = v.x[i+1]
            y1  = log(v.y[j])
            y2  = log(v.y[j+1])
            return bilinear((x1,x2),(y1,y2),(f11,f12,f21,f22),(x,y))
        except IndexError:
            return None



    def get_stream(vx, vy, x0, y0, nmax=1000000, maxlength=4*np.pi, bidirectional=True, reverse=False):
        """
        Function for computing a streamline.
        Parameters:
        -----------

        x0,y0: Floats.
              Initial position for the stream
        nmax: Integer.
              Maxium number of iterations for the stream.
        maxlength: Float
                   Maxium allowed length for a stream
        bidirectional=True
                      If it's True, the stream will be forward and backward computed.
        reverse=False
                The sign of the stream. You can change it mannualy for a single stream,
                but in practice, it's recommeneded to use this function without set reverse
                and setting bidirectional = True.

        Output:
        -------

        If bidirectional is False, the function returns a single array, containing the streamline:
        The format is:

                                          np.array([[x],[y]])

        If bidirectional is True, the function returns a tuple of two arrays, each one with the same
        format as bidirectional=False.
        The format in this case is:

                                (np.array([[x],[y]]),np.array([[x],[y]]))

        This format is a little bit more complicated, and the best way to manipulate it is with iterators.
        For example, if you want to plot the streams computed with bidirectional=True, you should write
        some similar to:

        stream = get_stream(x0,y0)
        ax.plot(stream[0][0],stream[0][1]) #Forward
        ax.plot(stream[1][0],stream[1][1]) #Backward

        """

        if bidirectional:
            s0 = get_stream(vx, vy, x0, y0, reverse=False, bidirectional=False, nmax=nmax,maxlength=maxlength)
            s1 = get_stream(vx, vy, x0, y0, reverse=True,  bidirectional=False, nmax=nmax,maxlength=maxlength)
            return (s0,s1)

        l = 0
        x = [x0]
        y = [y0]

        for i in xrange(nmax):
            ds = euler(vx, vy, x0, y0, reverse=reverse)
            if ds[0] == None:
                if(len(x)==1):
                    print "There was an error getting the stream, ds was NULL (see get_stream)."
                break
            l += np.sqrt(ds[0]**2+ds[1]**2)
            dx = ds[0]
            dy = ds[1]
            if(np.sqrt(dx**2+dy**2)<1e-13):
                print "Warning (get_stream): ds is very small, maybe you're in a stagnation point."
                print "Try selecting another initial point."
                break
            if l > maxlength:
                print "maxlength reached: ", l
                break
            x0 += dx
            y0 += dy
            x.append(x0)
            y.append(y0)
        return np.array([x,y])

    def get_random_streams(vx, vy, xmin=None, xmax=None, ymin=None, ymax=None, n=30, nmax=100000):
        if xmin == None:
            xmin = vx.xmin
        if ymin == None:
            ymin = vx.ymin
        if xmax == None:
            xmax = vx.xmax
        if ymax == None:
            ymax = vx.ymax
        X = xmin + np.random.rand(n)*(xmax-xmin)
        Y = ymin + np.random.rand(n)*(ymax-ymin)
        streams = []
        counter = 0
        for x,y in zip(X,Y):
            stream = get_stream(vx, vy, x, y, nmax=nmax, bidirectional=True)
            streams.append(stream)
            counter += 1
            print "stream ",counter, "done"
        return streams

    def plot_random_streams(streams, cartesian=False, **kargs):
        ax = plt.gca()
        print np.shape(streams)
        for stream in streams:
            for sub_stream in stream:
                if cartesian:
                    ax.plot(sub_stream[1]*cos(sub_stream[0]),sub_stream[1]*sin(sub_stream[0]),**kargs)
                else:
                    ax.plot(sub_stream[0],sub_stream[1],**kargs)

from scipy.interpolate import interp1d,interp2d


class Tensor(Mesh,Parameters):
    def __init__(self,rho,vx,vy,directory=''):
        if directory != '':
            if directory[-1] != '/':
                directory += '/'
        self.directory = directory

        Mesh.__init__(self,directory)
        self.params=Parameters(directory)

        self.rr = copy.deepcopy(rho)
        self.pp = copy.deepcopy(rho)
        self.rp = copy.deepcopy(vy)
        self.divp = copy.deepcopy(vx)
        self.divr = copy.deepcopy(vy)
        self.drFnu = copy.deepcopy(vx)
        nx = rho.nx
        ny = rho.ny
        vol = rho.vol
        ym = rho.ym
        xm  =rho.xm
        ymed = rho.ymed
        xmed =rho.xmed
        dx = diff(xm)[0]
        dy = diff(ym)
        dly = (dy/ymed)[0]


        nu = lambda x: self.params.alpha*self.params.aspectratio*self.params.aspectratio*(x)**(2*self.params.flaringindex+.5)


        rho = copy.copy(rho.data)
        vx = copy.copy(vx.data)
        vy = copy.copy(vy.data)

        for j in range(1,ny+2-1):   # 2 Ghost zones
            for i in range(nx):
                im = i-1 if i-1>=0 else nx-1
                ip = i+1 if i+1<nx else 0

                visc = nu(ymed[j])
                viscm = nu(ym[j])
                surfx = ym[j+1]-ym[j]
                surfy = ym[j]*dx
                surfyp = ym[j+1]*dx
                inv = 2/((ym[j+1]**2 - ym[j]**2)*dx)
                cs2 = self.params.aspectratio * ymed[j]**(self.params.flaringindex-.5)
                cs2 *= cs2
                div_v = 0

                div_v += (vx[j,ip]-vx[j,i])*surfx;
                div_v += (vy[j+1,i]*surfyp-vy[j,i]*surfy);
                div_v *= 2.0/3.0*inv;
                self.pp.data[j,i] = -rho[j,i]*cs2
                self.pp.data[j,i] += visc*rho[j,i]*(2.0*(vx[j,ip]-vx[j,i])/dx- div_v);
                self.pp.data[j,i] += visc*rho[j,i]*(vy[j+1,i]+vy[j,i])/ymed[j];
                self.rr.data[j,i] = -rho[j,i]*cs2
                self.rr.data[j,i] = visc*rho[j,i]*(2.0*(vy[j+1,i]-vy[j,i])/(ym[j+1]-ym[j]) - div_v);
                self.rp.data[j,i] = viscm*.25*(rho[j,i]+rho[j,im]+rho[j-1,i]+rho[j-1,im])*((vy[j,i]-vy[j,im])/(dx*ym[j]) + (vx[j,i]-vx[j-1,i])/(ymed[j]-ymed[j-1])-.5*(vx[j,i]+vx[j-1,i])/ym[j]);

        self.pp.data[0,:] = self.pp.data[2,:]
        self.pp.data[1,:] = self.pp.data[2,:]
        self.pp.data[-1,:] = self.pp.data[-2,:]
        self.rr.data[0,:] = self.rr.data[2,:]
        self.rr.data[1,:] = self.rr.data[2,:]
        self.rr.data[-1,:] = self.rr.data[-2,:]
        self.rp.data[0,:] = self.rp.data[2,:]
        self.rp.data[1,:] = self.rp.data[2,:]
        self.rp.data[-1,:] = self.rp.data[-2,:]

        for j in range(1,ny+2-1):   # 2 Ghost zones
            for i in range(nx):
                im = i-1 if i-1>=0 else nx-1
                ip = i+1 if i+1<nx else 0


                self.drFnu.data[j,i]  = (ym[j+1]*ym[j+1]*self.rp.data[j+1,i]-ym[j]*ym[j]*self.rp.data[j,i])/(ym[j+1]-ym[j]);
                self.divp.data[j,i] = 2.0*(self.pp.data[j,i]-self.pp.data[j,im])/(dx*(rho[j,i]+rho[j,im]));
                self.divp.data[j,i]  += 2.0*self.drFnu.data[j,i]/(ymed[j]*ymed[j]*(rho[j,im]+rho[j,i]));
                self.divp.data[j,i]  += 2.0*(ymed[j]*self.rr.data[j,i]-ymed[j-1]*self.rr.data[j-1,i])/((ymed[j]-ymed[j-1])*(rho[j,i]+rho[j-1,i])*ym[j]);
                self.divr.data[j,i] = 2.0*(ymed[j]*self.rr.data[j,i]-ymed[j-1]*self.rr.data[j-1,i])/((ymed[j]-ymed[j-1])*(rho[j,i]+rho[j-1,i])*ym[j]);
                self.divr.data[j,i]  += 2.0*(self.rp.data[j,ip]-self.rp.data[j,i])/(dx*ym[j]*(rho[j,i]+rho[j-1,i]));
                self.divr.data[j,i]  -= (self.pp.data[j,i]+self.pp.data[j-1,i])/(ym[j]*(rho[j,i]+rho[j-1,i]));

        self.divp.data[0,:] = self.divp.data[2,:]
        self.divp.data[1,:] = self.divp.data[2,:]
        self.divp.data[-1,:] = self.divp.data[-2,:]
        self.divr.data[0,:] = self.divr.data[2,:]
        self.divr.data[1,:] = self.divr.data[2,:]
        self.divr.data[-1,:] = self.divr.data[-2,:]
        self.drFnu.data[0,:] = self.drFnu.data[2,:]
        self.drFnu.data[1,:] = self.drFnu.data[2,:]
        self.drFnu.data[-1,:] = self.drFnu.data[-2,:]

        self.rr.recalculate()
        self.rr.name = 'Pirr'
        self.rr.math_name = '$\\Pi_{rr}$'
        self.pp.recalculate()
        self.pp.name = 'Pipp'
        self.pp.math_name = '$\\Pi_{\\phi\\phi}$'
        self.rp.recalculate()
        self.rp.name = 'Pirp'
        self.rp.math_name = '$\\Pi_{r\\phi}$'
        self.divr.recalculate()
        self.divp.name = 'divPip'
        self.divp.math_name = '$\\frac{(\\nabla \\cdot \\Pi)_\\phi}{\\Sigma}$'
        self.divp.recalculate()
        self.divr.name = 'divPir'
        self.divr.math_name = '$\\frac{(\\nabla \\cdot \\Pi)_r}{\\Sigma}$'
        self.drFnu.recalculate()
        self.drFnu.name = 'drFnu'
        self.drFnu.math_name = '$\\partial_r(r^2 \\Pi_{r\phi})$'
class Sim(Mesh,Parameters):
    def __init__(self,i,directory='',p=0,fargodir='/projects/p20783/amd616/fargo3d/'):
        if directory != '':
            if directory[-1] != '/':
                directory += '/'

        self.directory = directory
        self.fargodir = fargodir

        Mesh.__init__(self,directory)
        self.params=Parameters(directory)
        self.n = i
        self.t = i * self.params.dt * self.params.ninterm
        self.dens = Field('dens',i,directory=directory)
        self.vx = Field('vx',i,directory=directory,staggered='x')
        self.vr = Field('vy',i,directory=directory,staggered='y')


#        self.Pi = Tensor(self.dens,self.vp,self.vr)
#        self.Fw,self.Fd,self.Lambda,self.dTr = self.calc_torques()
#        self.vrf = self.mdot.avg/(-2*pi*self.mdot.y*self.rhos.avg)
#        self.vp = self.vp.shift_field('x')
#        self.vr = self.vr.shift_field('y')
#        self.vp = self.vp.cut_field(direction='y',side='p')
#        self.vr = self.vr.cut_field(direction='x',side='p')
#        self.dens = self.dens.cut_field(direction='xy', side='p')
        self.r  = self.dens.y
        self.phi = self.dens.x
        self.dphi = self.dx
#        self.dlr = diff(log(self.r))[0]
#        self.dr = self.dlr*self.r
        self.dr = self.dy
        self.dlr = self.dy/self.ymed

        self.dlr = self.dlr[:self.dens.ny-1]
        self.dr = self.dr[:self.dens.ny-1]



#        self.pp, self.rr = meshgrid(self.phi,self.r)
#    	self.dbar = self.dens.data.mean(axis=1)
#        self.vrbar0 = self.vr.data.mean(axis=1)
#    	self.vrbar = (self.dens.data*self.vr.data).mean(axis=1)
#    	self.vpbar = (self.dens.data*self.vp.data).mean(axis=1)
#    	self.mdot_full = -2*pi*self.rr*(self.vr.data*self.dens.data)
#        self.mdot = self.mdot_full.mean(axis=1)
        try:
    	    _,self.px,self.py,self.pz,self.pvx,self.pvy,self.pvz,self.mp,_,self.omf  = loadtxt(directory+'planet{0:d}.dat'.format(p))[i,:]
        except IndexError:
            try:
    	        _,self.px,self.py,self.pz,self.pvx,self.pvy,self.pvz,self.mp,_,self.omf  = loadtxt(directory+'planet{0:d}.dat'.format(p))[-1,:]
            except IndexError:
    	        _,self.px,self.py,self.pz,self.pvx,self.pvy,self.pvz,self.mp,_,self.omf  = loadtxt(directory+'planet{0:d}.dat'.format(p))


    	self.a = sqrt(self.px**2  + self.py**2)
        self.K = self.mp**2 / (self.dens.alpha * self.dens.aspectratio**5)
        self.sig0  = self.dens.mdot/(3*pi*self.nu(self.dens.y))
        self.torque_norm = .266 * self.params.aspectratio**3 * self.params.mdot * self.K
        self.wkb_torque = (self.a/abs(self.dens.y-self.a))**4
        self.wkb_torque[self.dens.y <  self.a] *= -1
        self.wkb_torque[abs(self.dens.y-self.a)<=1.3*self.scaleH(self.a)] = 0
        self.wkb_torque *= 1e-2*self.sig0 *pi*self.a *self.mp**2

        self.wkb_torque /= (self.dens.y*2*pi)

      #  self.load_fluxes(i)
        self.vp = copy.deepcopy(self.vx)
        self.vp.data += self.vp.y[:,newaxis]*self.omf
        self.vp.recalculate()

        #self.mdot,self.rhos,self.rhosp = self.calc_flux()
        try:
            self.nu0 = self.dens.alpha*self.dens.aspectratio*self.dens.aspectratio
        except AttributeError:
            self.nu0 = self.dens.alphaviscosity*self.dens.aspectratio*self.dens.aspectratio


    	self.tvisc = self.dens.ymax**2/self.nu(self.dens.ymax)
        self.tviscp = self.a**2/(self.nu(self.a))
        self.torb = 2*np.pi * self.a**(1.5)
        self.rh = (self.mp/3)**(1./3) * self.a
        self.safety_fac  = .5

        self.rayleigh = np.zeros(self.vp.data.shape)
        dat = self.vp.data*self.vp.y[:,newaxis]
        self.rayleigh[1:-1,:] = (dat[2:,:] - dat[:-2,:])/((self.vp.y[2:]-self.vp.y[:-2])[:,newaxis])
        self.rayleigh[0,:] = (dat[1,:]-dat[0,:])/(self.vp.y[1]-self.vp.y[0])
        self.rayleigh[-1,:] = (dat[-1,:]-dat[-2,:])/(self.vp.y[-1]-self.vp.y[-2])
        self.rayleigh = self.rayleigh.mean(axis=1)


        self.vortensity = np.zeros(self.vp.data.shape)
        dx = self.dens.dx[0]
        dat = self.vp.data*self.vp.y[:,newaxis]
        self.vortensity[1:-1,:] = (dat[2:,:] - dat[:-2,:])/((self.vp.y[2:]-self.vp.y[:-2])[:,newaxis])
        self.vortensity[0,:] = (dat[1,:]-dat[0,:])/(self.vp.y[1]-self.vp.y[0])
        self.vortensity[-1,:] = (dat[-1,:]-dat[-2,:])/(self.vp.y[-1]-self.vp.y[-2])
        dat = self.vr.data
        self.vortensity[:,1:-1] -= (dat[:,2:] - dat[:,:-2])/(2*dx)
        self.vortensity[:,0] -= (dat[:,1]-dat[:,0])/dx
        self.vortensity[:,-1] -= (dat[:,-1]-dat[:,-2])/dx

        self.vortensity /= (self.dens.data*self.dens.y[:,newaxis])

        self.rwi = (.5/self.vortensity).mean(axis=1)





        try:
            trqname = 'torque%d.dat'%self.n
            dat = np.fromfile(self.directory+trqname)
            print 'Loaded torque file from %s'%trqname
            attrs = ['y','dbar','Lt','Ltn','Ld','Ldn','Lw','Lwn',
                    'LdS','LdSn',
                    'drFt','drFd','mdotl','drFnu','drFdB','drFw', 'drFwB',
                    'Lamex','Lamdep','LamdepB',
                    'dtLt','dtLd','dtLw',
                    'dtLds','dtdbar','mdot',
                    'Lamdep1','Lamdep2','Lamdep3','Lamdep4','Lamdep5','Lamdep6','Lamdep7',
                    'dtLd_rhs']

            for i,a in enumerate(attrs):
                setattr(self,'trq_'+a,dat[i*self.params.ny:(i+1)*self.params.ny])
            self.trq_Lamtot = np.zeros(self.trq_Lamdep.shape)
            for i in range(1,8):
                self.trq_Lamtot += getattr(self,'trq_Lamdep%d'%i)



            ind = (self.trq_y>self.dens.wkz)&(self.trq_y<self.dens.wkzr)
            indL = (self.trq_y <= self.dens.wkz + self.dlr[0]*self.dens.wkz)
            indR = (self.trq_y >= self.dens.wkzr - self.dlr[0]*self.dens.wkzr)
            self.trq_intLam = np.zeros(self.trq_y.shape)
            self.trq_intLam[ind] = (2*np.pi*self.dens.y*self.dens.dy*self.trq_Lamdep)[ind].cumsum()
            self.trq_intLam[ind] -= self.trq_intLam[ind][0]
            self.trq_intLam[indR] = self.trq_intLam[ind][-1]
            self.sig_ss = 1 + self.trq_intLam/(self.params.mdot * np.sqrt(self.trq_y))

            self.trq_intmdotl = np.zeros(self.trq_y.shape)
            self.trq_intmdotl[ind] = (2*np.pi*self.dens.y*self.dens.dy*self.trq_mdotl)[ind].cumsum()
            self.trq_intmdotl[ind] -= self.trq_intmdotl[ind][0]
            self.trq_intmdotl[self.trq_y>=self.dens.wkzr] = self.trq_intmdotl[ind][-1]

            self.trq_intdrFnu = np.zeros(self.trq_y.shape)
            self.trq_intdrFnu[ind] = (2*np.pi*self.dens.y*self.dens.dy*self.trq_drFnu)[ind].cumsum()
            self.trq_intdrFnu[ind] -= self.trq_intdrFnu[ind][0]
            self.trq_intdrFnu[self.trq_y>=self.dens.wkzr] = self.trq_intdrFnu[ind][-1]

            self.trq_intdtLd = np.zeros(self.trq_y.shape)
            self.trq_intdtLd[ind] = (2*np.pi*self.dens.y*self.dens.dy*self.trq_dtLd)[ind].cumsum()
            self.trq_intdtLd[ind] -= self.trq_intdtLd[ind][0]
            self.trq_intdtLd[self.trq_y>=self.dens.wkzr] = self.trq_intdtLd[ind][-1]

            dat = np.fromfile(self.directory+'torque_m%d.dat'%self.n)
            self.mmax = 30
            self.trq_drFdm = dat[:self.params.ny*(self.mmax+2)].reshape(self.mmax+2,self.params.ny)
            self.trq_drFdm[1,:] = self.trq_drFd - self.trq_drFdm[2:,:].sum(axis=0)

            self.trq_Lamexm = dat[self.params.ny*(self.mmax+2):2*self.params.ny*(self.mmax+2)].reshape(self.mmax+2,self.params.ny)

            self.trq_dtLtm = dat[2*self.params.ny*(self.mmax+2):3*self.params.ny*(self.mmax+2)].reshape(self.mmax+2,self.params.ny)
            self.trq_dtLtm[1,:] = self.trq_dtLt - self.trq_dtLtm[2:,:].sum(axis=0)

            self.trq_drFtm = self.trq_Lamexm - self.trq_dtLtm
            self.trq_drFwm =self.trq_drFtm - self.trq_drFdm
            self.trq_dtLdm = np.zeros(self.trq_drFdm.shape)
            self.trq_dtLdm[0,:] = self.trq_dtLd
            self.trq_dtLdm[1,:] = self.trq_dtLd
            self.trq_dtLwm = self.trq_dtLt - self.trq_dtLdm
            self.trq_Lamdepm = self.trq_dtLdm + self.trq_drFdm
        except IOError:
            pass

        print 'The time is %.2e, %.2e outer viscous times, %.2e planet viscous times, %.2e planet orbits' %(self.t,self.t/self.tvisc, self.t/self.tviscp,self.t/self.torb)
#        self.dTr,self.Lambda,self.Fh=self.calc_torques()

##    	self.omega = zeros(self.vp.data.shape)
##    	self.vpc = zeros(self.vp.data.shape)
##        self.l = zeros(self.vp.data.shape)
##        for i in range(len(self.r)):
##            self.omega[i,:] = self.vp.data[i,:]/self.r[i] + self.a**(-1.5)
##            self.vpc[i,:] = self.omega[i,:]*self.r[i]
##            self.l[i,:] = self.r[i] * self.vpc[i,:]
##
##        self.Fj_full = -2*pi*self.dens.data*self.nu(self.rr)*self.rr**3 * self.grad(self.omega)
##        self.Fj = self.Fj_full.mean(axis=1)
##        ind_excl = (self.r>=self.a*(1-self.dens.thicknesssmoothing*self.dens.aspectratio))&(self.r<=self.a*(1+self.dens.thicknesssmoothing*self.dens.aspectratio))
##        self.dTr = (-self.rr*self.dens.data * self.dp_potential(self.rr,self.pp))
##        self.dTr_excl = self.dTr.mean(axis=1)
##        self.dTr_excl[ind_excl] = 0
##        self.dTr_tot_excl = 2*pi*trapz(self.dTr_excl* self.dr)
##        self.dTr_mean = self.dTr.mean(axis=1)
##        self.dTr_total = 2*pi*trapz(self.dTr_mean*self.dr)
##        self.dTr_out = self.dphi*trapz((self.dTr.sum(axis=1)*self.dr)[self.r>=self.a],x=self.r[self.r>=self.a])
##        self.dTr_in = self.dphi*trapz((self.dTr.sum(axis=1)*self.dr)[self.r<=self.a],x=self.r[self.r<=self.a])
##
##        self.mdotdlr = (self.mdot_full * self.grad(self.l)).mean(axis=1)
##
##        ##        try:
##            self.dbar0 = self.dens.mdot/(3*pi*self.nu(self.r))
##        except AttributeError:
##            pass
##        self.vr_visc = -1.5*self.nu(self.r)/self.r
##        self.steady_state_calc()
##        self.mdoti,self.dbar0=self.inner_disk_sol((self.dens.ymin*1.2,self.a*.5))
##
##        self.A = self.dTr_total/self.mdoti
##        self.A_excl = self.dTr_tot_excl/self.mdoti

    def steadystate_plot(self,ylims=None,savefig=None,logx=False,logy=False):
        fig,axes = plt.subplots(3,1,sharex=True,figsize=(10,10))

        ind = (self.trq_y > self.dens.wkz)&(self.trq_y < self.dens.wkzr)
        indL = (self.trq_y <= self.dens.wkz + self.dlr[0]*self.dens.wkz)
        indR = (self.trq_y >= self.dens.wkzr - self.dlr[0]*self.dens.wkzr)

        intLam = (self.trq_y*self.dens.dy*2*pi*self.trq_Lamdep).cumsum()
        intLam -= intLam[0]
        temp2 = (self.trq_y*self.dens.dy*2*pi*self.trq_Lamdep)[ind].cumsum()
        temp2 -= temp2[0]
        intLam2 = np.zeros(self.trq_y.shape)
        intLam2[ind] = temp2
        intLam2[indR] = temp2[-1]

        sig = 1 + intLam/(self.dens.mdot *np.sqrt(self.trq_y))
        sig2 = 1 + intLam2/(self.dens.mdot *np.sqrt(self.trq_y))

        mdot = self.trq_mdot/self.params.mdot
        mdot2 = self.trq_mdot[ind]/self.params.mdot
        mdotL = self.trq_mdot[indL]/self.params.mdot
        mdotR = self.trq_mdot[indR]/self.params.mdot

        ldep = self.trq_Lamdep.copy()
        ldep2 = self.trq_Lamdep[ind]
        ldepR = self.trq_Lamdep[indR]
        ldepL = self.trq_Lamdep[indL]
        lex = self.trq_Lamex.copy()
        lex2 = self.trq_Lamex[ind]
        lexL = self.trq_Lamex[indL]
        lexR = self.trq_Lamex[indR]

        drFw = self.trq_drFw.copy()
        drFw2 = self.trq_drFw[ind]
        drFwL = self.trq_drFw[indL]
        drFwR = self.trq_drFw[indR]

        y = self.trq_y.copy()
        y2 = self.trq_y[ind]
        yL = self.trq_y[indL]
        yR = self.trq_y[indR]

        dbar = self.dens.avg/self.sig0
        wkzL = self.dens.wkz
        wkzR = self.dens.wkzr


        axes[0].plot(y,dbar,'.k',label='Fargo, $t=%.3f t_{visc}^p$'%(s.t/s.tviscp))
        axes[0].plot(y,sig2,'-b',linewidth=2,label='Steady State')

        axes[1].plot(y2,lex2,'-k',label='$\\Lambda_{ex}$')
        axes[1].plot(y2,ldep2,'-b',label='$\\Lambda_{dep}$')
        axes[1].plot(y2,drFw2,'-r',label='$\\nabla_rF_w$')

        axes[1].plot(yL,lexL,'--k')
        axes[1].plot(yL,ldepL,'--b')
        axes[1].plot(yL,drFwL,'--r')
        axes[1].plot(yR,lexR,'--k')
        axes[1].plot(yR,ldepR,'--b')
        axes[1].plot(yR,drFwR,'--r')

        axes[2].plot(y2,mdot2,'-k')
        axes[2].plot(yL,mdotL,'-k')
        axes[2].plot(yR,mdotR,'-k')
        axes[2].axhline(1,color='k',linestyle='--')

        axes[2].set_xlabel('$r$',fontsize=20)
        axes[2].set_ylim(-4,4)
        axes[0].axvline(wkzL,linestyle='--',color='k')
        axes[1].axvline(wkzL,linestyle='--',color='k')
        axes[2].axvline(wkzL,linestyle='--',color='k')
        axes[0].axvline(wkzR,linestyle='--',color='k')
        axes[1].axvline(wkzR,linestyle='--',color='k')
        axes[2].axvline(wkzR,linestyle='--',color='k')

        axes[0].set_xlim(0,self.dens.ymax)


        axes[0].set_ylabel('$\\Sigma/\\Sigma_0$',fontsize=20)
        axes[2].set_ylabel('$\\dot{M}/\\dot{M}_o$',fontsize=20)

        axes[0].legend(loc='lower right')
        axes[1].legend(loc='lower right')

        if ylims is not None:
            axes[1].set_ylim(ylims)
        if logy:
            axes[0].set_yscale('log')

        if logx:
            for ax in axes.flatten():
                ax.set_xscale('log')


        for ax in axes.flatten():
            ax.minorticks_on()
        plt.subplots_adjust(wspace=0)
        plt.setp(axes[2].get_xticklabels()[-1],visible=False)

        if savefig is not None:
            fig.savefig(savefig)

    def torque_summary(self,axes=None,log=False,logx=False,logy=True,savefig=None,tstr=None,exclude=True,ylims=None):
        if axes is None:
            fig,axes=plt.subplots(2,2,sharex='col',figsize=(14,9))

        #if tstr is None:
        #    tstr = '$q=%.e, \\alpha=%.e, K = %.2e$'%(self.mp,self.params.alpha,self.K)

        if exclude:
            ind = (self.dens.y > self.dens.wkz)&(self.dens.y<self.dens.wkzr)
        else:
            ind = np.ones(self.dens.y.shape).astype(bool)

        axes[0,0].plot(self.dens.y[ind],self.dens.avg[ind]/self.sig0[ind])
        axes[0,1].plot(self.trq_y[ind],self.trq_mdot[ind]/self.params.mdot)

        axes[1,0].plot(self.trq_y[ind],self.trq_Lamdep[ind],'-b')
        axes[1,0].plot(self.trq_y[ind],self.trq_Lamex[ind],'-k')
        axes[1,0].plot(self.trq_y[ind],self.trq_drFw[ind],'-r')
        axes[1,0].plot(self.trq_y[ind],self.trq_dtLw[ind],'-m')

        axes[1,1].plot(self.trq_y[ind],self.trq_Lamdep[ind],'-b')
        #axes[1,1].plot(self.trq_y[ind],self.trq_mdotl[ind],'-g')
        axes[1,1].plot(self.trq_y[ind],self.trq_drFnu[ind],'-r')
        #axes[1,1].plot(self.trq_y[ind],self.trq_dtLd[ind],'-m')
        axes[1,1].plot(self.trq_y[ind],self.trq_mdotl[ind]+self.trq_dtLd[ind],'-k')

        axes[1,0].set_xlabel('$r$',fontsize=20)
        axes[1,1].set_xlabel('$r$',fontsize=20)

        axes[0,1].set_ylim(-5,5)

        #axes[0,0].set_title(tstr,fontsize=20)

        axes[0,1].axhline(1,color='k')

        axes[0,0].set_ylabel('$\\Sigma/\\Sigma_0$',fontsize=20)
        axes[0,1].set_ylabel('$\\dot{M}/\\dot{M}_0$',fontsize=20)

        axes[1,0].set_ylabel('Torques',fontsize=15)

        if ylims is not None:
            axes[1,0].set_ylim(ylims)

        if log or logy:
            axes[0,0].set_yscale('log')

        if log or logx:
            for ax in axes.flatten():
                ax.set_xscale('log')

        for ax in axes.flatten():
            ax.minorticks_on()

        fig.canvas.draw()

        if savefig is not None:
            print 'Saving figure to ', savefig
            fig.savefig(savefig)

        return

    def calc_mdot(self):
        dqm = zeros(self.dens.data.shape)
        dqp = zeros(self.dens.data.shape)
        slope = zeros(self.dens.data.shape)
        mdot = zeros(self.dens.data.shape)
        denstar = zeros(self.dens.data.shape)

        dy = (self.dens.ym[1:]-self.dens.ym[:-1])
        r = self.dens.ymed

        dqm[1:-1,:] =(self.dens.data[1:-1,:]  - self.dens.data[:-2,:])/dy[1:-1,newaxis]
        dqp[1:-1,:] =(self.dens.data[2:,:]  - self.dens.data[1:-1,:])/dy[2:,newaxis]
        ind = sign(dqm*dqp)
        slope = (ind>0).astype(int) * 2*dqm*dqp/(dqm+dqp)
        denstar[1:-1,:] = (self.vr.data[1:-1,:]>0).astype(int)*(self.dens.data[:-2,:]+.5*slope[:-2,:]*dy[:-2,newaxis]) + (self.vr.data[1:-1,:]<=0).astype(int)*(self.dens.data[1:-1,:]-.5*slope[1:-1,:]*dy[1:-1,newaxis])

        mdot = denstar * self.vr.data * -2*pi*self.dens.ym[:-1,newaxis]
        return mdot

    def dens_evolution(self):

        dat = np.fromfile(self.directory + 'mass_1d_Y_raw.dat')
        nt = len(dat)/self.params.ny
        dat = dat.reshape(nt,self.params.ny)

        for i in range(nt):
            dat[i,:] /= (2*np.pi*self.vol[:,0] * self.sig0)


        return self.params.dt*arange(nt), dat.min(axis=1)

    def load_fluxes(self,i):
        execdir = self.fargodir+'utils/python/'
        execname = 'load_fluxes'
        call(['mkdir',self.directory + 'temp_files'])

        lines='%d\n' % self.params.nx
        lines +='%d\n' % self.params.ny
        lines +='%d\n' %  self.params.nz
        lines +='%f\n' %  self.params.alpha
        lines +='%f\n' %  self.mp
        lines +='%f\n' %  self.a
        lines +='%f\n' %  self.omf
        lines +='%f\n' %  self.params.aspectratio
        lines +='%f\n' %  self.params.flaringindex
        lines +='%f\n' %  self.params.mdot
        lines +='%f\n' %  self.params.thicknesssmoothing

        with open(self.directory + 'param_file.txt','w') as f:
            f.write(lines)

        callstr = [execdir+execname,'%d'%i,self.directory]
        print ' '.join(callstr)
        call(callstr)


        self.fw = Field('fw',0,directory=self.directory,staggered='y')
        self.fd = Field('fd',0,directory=self.directory,staggered='y')
        self.lambda_dep = Field('lambda_dep',0,directory=self.directory,staggered='y')
        self.lambda_ex = Field('lambda_ex',0,directory=self.directory,staggered='y')
        self.mdot = Field('mdot',0,directory=self.directory,staggered='y')
        self.rhostar = Field('rhostar',0,directory=self.directory,staggered='y')
        self.taurr = Field('taurr',0,directory=self.directory)
        self.taupp = Field('taupp',0,directory=self.directory)
        self.taurp = Field('taurp',0,directory=self.directory,staggered='xy')
        self.lstar = Field('lstar',0,directory=self.directory,staggered='y')
        self.pot = Field('pot',0,directory=self.directory,staggered='y')
        self.rhoslope = Field('rhoslope',0,directory=self.directory,staggered='y')

        self.Ld = Field('Ld',0,directory=self.directory)
        self.Lw = Field('Lw',0,directory=self.directory)

        callstr = ['rm',self.directory+'param_file.txt']
        call(callstr)
        callstr = ['rm','-rf',self.directory+'temp_files/']
        call(callstr)



        indL = self.lambda_ex.y<=self.lambda_ex.wkz
        indR = self.lambda_ex.y>self.lambda_ex.wkz
        self.lambda_ex.integral = cumtrapz(self.lambda_ex.avg[indR]*self.lambda_ex.dy[indR])
        self.lambda_ex.integral = hstack( ( zeros(self.lambda_ex.ny-len(self.lambda_ex.integral)),self.lambda_ex.integral))
        indL = self.lambda_dep.y<=self.lambda_dep.wkz
        indR = self.lambda_dep.y>self.lambda_dep.wkz
        self.lambda_dep.integral = cumtrapz(self.lambda_dep.avg[indR]*self.lambda_dep.dy[indR])
        self.lambda_dep.integral = hstack( ( zeros(self.lambda_dep.ny-len(self.lambda_dep.integral)),self.lambda_dep.integral))

        return

    def dump_torque(self,exclude_wkz=False):
        fname = self.directory + 'trq_%d_out.dat'%self.n

        if exclude_wkz:
            ind = (self.trq_y > self.dens.wkz)&(self.trq_y < self.dens.wkzr)
        else:
            #ind = range(len(self.trq_y))
            ind = range(10,len(self.trq_y)-10)


        header = np.array( [ float(len(self.trq_y[ind])),
            self.params.alpha,
            self.params.mdot,
            self.params.aspectratio,
            self.params.ymin,
            self.params.ymax])


        np.hstack((header,
            self.trq_y[ind],
            self.trq_Lamdep[ind],
            self.trq_Lamex[ind],
            self.trq_mdot[ind],
            self.trq_drFd[ind],
            self.trq_drFw[ind],
            self.trq_Ld[ind],
            self.trq_dbar[ind]) ).tofile(fname)

    def plot_resonance(self,ax=None):
        if ax is None:
            fig=figure()
            ax=fig.add_subplot(111)
        omd = self.vx.data/self.ymed[:,newaxis] + self.omf
        kapd = (omd[2:,:]*self.ymed[2:,newaxis]**2 - omd[:-2,:]*self.ymed[:-2,newaxis]**2)/(2*self.dlr[0])
        omd = omd[1:-1,:]
        y = self.ymed[1:-1]
        kapd *= 2*omd/((y*y)[:,newaxis])
        omd = omd.mean(axis=1)
        kapd = kapd.mean(axis=1)
        res = (kapd/(omd - self.omf)**2)
        x = (y-self.a)/(self.params.aspectratio*self.a)
        ax.plot(x[res>0],sqrt(res[res>0]))
        ax.set_xlabel('$(r-a)/H(a)$',fontsize=15)
        ax.set_ylabel('$m$',fontsize=15)
        ax.set_yscale('log')


    def fpred(self,x):
        scalar = False
        try:
            x[1]
        except IndexError:
            scalar = True
        ind = self.r<=x
        res = self.dTr_mean*2*pi + self.mdotdlr
        if scalar:
            return trapz( (res*self.dr)[ind])
        else:
            return array([trapz( (res*self.dr)[self.r<=i]) for i in x])


    def calc_torques(self):
        res_dTr = copy.deepcopy(self.rhos)

        res_dTr.data *=  (-res_dTr.yy * self.dp_potential(yy,xx))

        res_dTr.name = 'Lambda_ex'
        res_dTr.math_name = '$\\Lambda_{ex}$'
        res_dTr.recalculate()

        ym = self.dens.ym
        xm = self.dens.xm
        ymed = self.dens.ymed
        xmed = self.dens.xmed

        vx = copy.copy(self.vp.data)
        vy = copy.copy(self.vr.data)
        rho = copy.copy(self.dens.data)

        for j in range(1,self.dens.ny+2-1):
            for i in range(nx):
                im = i-1 if i-1>=0 else nx-1
                ip = i+1 if i+1<nx else 0
                l_ed = (ym[j]*vx[j,im] + ym[j]*vx[j,ip] + ym[j-1]*vx[j-1,ip] + ym[j-1]*vx[j-1,im])*.25
                dlr_ed = .5*((ym[j]*vx[j,im] - ym[j-1]*vx[j-1,im])/(ym[j]-ym[j-1]) + (ym[j]*vx[j,ip] - ym[j-1]*vx[j-1,ip])/(ym[j]-ym[j-1]))
                prp = .5*(self.Pi.rp[j,im] + self.Pi.rp[j,ip])
                dpp = .25*( self.Pi.divp[j,im] + self.Pi.divp[j,ip] + self.Pi.divp[j-1,im] + self.Pi.divp[j-1,ip] )

                mdl = self.mdot.avg[j]*dlr_ed
                res_Fw[j,i] = -self.mdot[j,i] + mdl
                res_Fd[j,i] =  - 2*pi*ym[j]*ym[j]*prp - mdl


        res_Trp = copy.deepcopy(self.Pi.rp)
        res_Fd = copy.deepcopy(self.rhosp)
        res_Fw = copy.deepcopy(self.rhosp)

        fac = -self.mdot.avg*self.vp.avg*self.vp.y
        res_Fd.data = -fac[:,newaxis] - 2*pi*self.Pi.rp.data *self.Pi.rp.yy
        res_Fw.data  =  - (self.mdot*self.vp*self.vp.yy +fac[:,newaxis])
        res_Trp.data = res_Fd.data + res_Fw.data

        res_dep = copy.deepcopy(self.dens)
        # Assume a keplerian rotation for now
        res_Fd.data = self.mdot.avg*self.rhos.avg (3*pi*self.nu(yy) * yy**(.5))
        res_Fd.name = 'Fh'
        res_Fd.math_name = '$F_{H,\\nu}$'
        res_Fd.recalculate()
        res_dep.data  = -self.mdot.data /(2*sqrt(yy)) + res_fH.grad()
        res_dep.name = 'Lambda_dep'
        res_dep.math_name = '$\\Lambda_{dep}$'
        res_dep.recalculate()

        return res_dTr,res_dep,res_fH
    def calc_dt(self):
        cfl = self.params.cfl
        c2 = 4*sqrt(2)
        vy = abs(self.vr.data)
        vx = self.vp.data - self.omf*self.vp.yy
        vx = (vx - vx.mean(axis=1)[:,newaxis])
        visc = 4*self.nu(self.dens.yy)

        dt2_y = (self.vr.surfy/vy).min()
        dt2_x = (self.vp.surfx/vx).min()
        dt4_y = (self.dens.surfy/visc).min()
        dt4_x = (self.dens.surfx/visc).min()

        dt3_x = (c2*abs(self.vp.surfx[:,:-1]/(vx[:,1:]-vx[:,:-1]))).min()
        dt3_y = (c2*abs(self.vr.surfy[:-1,:]/(vy[1:,:]-vx[:-1,:]))).min()

        dt_x = sqrt( dt2_x**(-2) + dt3_x**(-2) + dt4_x**(-2) )
        dt_y = sqrt( dt2_y**(-2) + dt3_y**(-2) + dt4_y**(-2) )

        return cfl*min(1/dt_x,1/dt_y)

    def momentum_flux(self):

        rho = copy.copy(self.dens.data)
        vy = copy.copy(self.vr.data)
        vx = copy.copy(self.vp.data)

        surfx = self.dens.surfx
        surfy = self.dens.surfy
        vol = self.dens.vol
        ymed = self.dens.ymed
        ym = self.vr.ym[:-1]
        xmed =self.dens.xmed
        xm = self.vp.xm[:-1]

        vx -= vx.mean(axis=1)[:,newaxis]

        Piy_m = rho[:-1,:]*vy[:-1,:]
        Piy_p = rho[:-1,:]*vy[1:,:]
        Pix_m = rho[:,:-1]*vx[:,:-1]
        Pix_p = rho[:,:-1]*vx[:,1:]





        # First y
        dxx,dyy = meshgrid(self.dx,self.dy)
        fluxp = self.transport(Piy_p,vy[1:,:],dyy[:-1,:],surfy[1:,:])
        fluxm = self.transport(Piy_m,vy[1:,:],dyy[:-1,:],surfy[1:,:])

#        Q =  copy.copy(Piy_p)
#        Qs = copy.copy(Piy_p)
#        v = vy[:-1,:]
#        dqm = zeros(Q.shape)
#        dqp = zeros(Q.shape)
#        slope = zeros(Q.shape)
#        flux = zeros(Q.shape)
#        _,dyy = meshgrid(self.dx,self.dy[:-1])
#        surfy = copy.copy(surfy[:-1,:])
#
#
#        dqm[1:-1,:] = (Q[1:-1,:]-Q[:-2,:])/dyy[1:-1,:]
#        dqp[1:-1,:] = (Q[2:,:]-Q[1:-1,:])/dyy[2:,:]
#        slope[1:-1,:] = ( (dqm*dqp)[1:-1,:] > 0).astype(int) * (2*(dqp*dqm)/(dqp+dqm))[1:-1,:]
#        Qs[1:-1,:] = (v[1:-1,:]>0).astype(int) * (Q[:-2,:]+.5*slope[:-2,:]*dyy[:-2,:])
#        Qs[1:-1,:] += (v[1:-1,:]<=0).astype(int) * (Q[1:-1,:]-.5*slope[1:-1,:]*dyy[1:-1,:])
#        flux[1:-1,:] = -self.nx*v[1:-1,:] * Qs[1:-1,:] * surfy[1:-1,:]

        return fluxp,fluxm

    def transport(self,Q,v,dx,surfx):
        Qs = copy.copy(Q)
        dqm = zeros(Q.shape)
        dqp = zeros(Q.shape)
        slope = zeros(Q.shape)
        flux = zeros(Q.shape)


        dqm[1:-1,:] = (Q[1:-1,:]-Q[:-2,:])/dx[1:-1,:]
        dqp[1:-1,:] = (Q[2:,:]-Q[1:-1,:])/dx[2:,:]
        slope[1:-1,:] = ( (dqm*dqp)[1:-1,:] > 0).astype(int) * (2*(dqp*dqm)/(dqp+dqm))[1:-1,:]
        Qs[1:-1,:] = (v[1:-1,:]>0).astype(int) * (Q[:-2,:]+.5*slope[:-2,:]*dx[:-2,:])
        Qs[1:-1,:] += (v[1:-1,:]<=0).astype(int) * (Q[1:-1,:]-.5*slope[1:-1,:]*dx[1:-1,:])
        flux[1:-1,:] = -v[1:-1,:] * Qs[1:-1,:] * surfx[1:-1,:]
        return flux

    def calc_visc_tens(self):
        vr = copy.copy(self.vp.data)/self.vp.yy
        vp = copy.copy(self.vr.data)
        rho = copy.copy(self.dens.data)
        taurp = self.dens.center_to_edge('xy')

        #ddx,ddy=meshgrid(self.dx,self.dy)
        #dlr = log(self.dy)[0]
        #rho = self.dens.center_to_edge('xy')

        #dvrdp = (vy[:,1:]-vy[:,:-1])/self.dx[0]
        #dvrdp /= self.vr.yy[:,:-1]
        #dvpdr = (vx[1:,:]-vx[:-1,:])/dlr
        #taurp.data[1:-1,1:-1] = self.nu(taurp.yy[1:-1,:1:-1])*rho.data[1:-1,1:-1]*(dvrdp[1:-1,:-1]+dvpdr[:-1,1:-1])


        dx = self.dx[0]
        dy = self.dy
        ymed = self.dens.ymed
        ym = self.vr.ym[:-2]
        surfx = self.dens.surfx
        surfy = self.dens.surfy
        vol = self.dens.vol
        invym = 1./ym
        rhoc = .25*(rho[1:,1:] + rho[:-1,1:] + rho[:-1,:-1] + rho[1:,:-1])
        dvp = vp[1:,:-1]-vp[:-1,:-1]
        dvp /= (ymed[1:]-ymed[:-1])[:,newaxis]
        dvr = vr[:-1,1:]-vr[:-1,:-1]
        dvr *= invym[:,newaxis]
        dvr /= dx
        visc = self.nu(ym)
        taurp.data[1:,1:] =rhoc*( dvr + dvp -.5*(vp[1:,:-1]+vp[:-1,:-1])*invym[:,newaxis])*visc[:,newaxis]
        taurp.data[0,:] = taurp.data[1,:]
        taurp.data[:,0] = taurp.data[:,1]
        taurp.staggered='xy'
        taurp.y = taurp.ym[:-1]
        taurp.x = taurp.xm[:-1]
        taurp.recalculate()

        taupp = copy.deepcopy(self.dens)
        div_v = (vp[1:,1:]-vp[1:,:-1])*surfx[1:,:]
        div_v += (vr[2:,:-1]*surfy[2:,:-1] - vr[1:,:-1]*surfy[1:,:-1])
        div_v *= (2./3)/vol[1:,:-1]
        visc = self.nu(self.dens.yy[1:,:-1])
        taupp.data[1:,:-1] = 2*(vp[1:,1:]-vp[1:,:-1])/dx - div_v
        taupp.data[1:,:-1] += (vr[2:,:-1]-vr[1:,:-1])/self.dens.yy[1:,:-1]
        taupp.data[1:,:-1] *= visc*self.dens.data[1:,:-1]
        taupp.recalculate()


        return taurp,taupp

    def calc_flux(self):
        res = copy.deepcopy(self.vr)
        res1 = copy.deepcopy(self.vr)
        rho = copy.copy(self.dens.data)
        vy = copy.copy(self.vr.data)
        surfy = copy.copy(self.surfy)

        dqm = zeros(rho.shape)
        dqp = zeros(rho.shape)
        slope = zeros(rho.shape)
        mdot = zeros(rho.shape)
        rhos = copy.copy(self.dens.data)
        _,dyy = meshgrid(self.dx,self.dy)

        dqm[1:-1,:] = (rho[1:-1,:]-rho[:-2,:])/dyy[1:-1,:]
        dqp[1:-1,:] = (rho[2:,:]-rho[1:-1,:])/dyy[2:,:]
        slope[1:-1,:] = ( (dqm*dqp)[1:-1,:] > 0).astype(int) * (2*(dqp*dqm)/(dqp+dqm))[1:-1,:]
        rhos[1:-1,:] = (vy[1:-1,:]>0).astype(int) * (rho[:-2,:]+.5*slope[:-2,:]*dyy[:-2,:])
        rhos[1:-1,:] += (vy[1:-1,:]<=0).astype(int) * (rho[1:-1,:]-.5*slope[1:-1,:]*dyy[1:-1,:])
        mdot[1:-1,:] = -self.nx*vy[1:-1,:] * rhos[1:-1,:] * surfy[1:-1,:]
        mdot[0,:] = mdot[1,:]
        mdot[-1,:] = mdot[-2,:]

        res.data = mdot
        res.name = 'Mdot'
        res.math_name = '$\\dot{M}$'
        res1.data = rhos
        res1.name = 'Rhostar'
        res1.math_name = '$\\Sigma^*$'
        res.recalculate()
        res1.recalculate()

        # phi direction
        res2 = copy.deepcopy(self.vp)
        rho = copy.copy(self.dens.data)
        ns = (len(rho[:,0]),1)
        dx = self.dx[0]


        nrho = hstack( (rho[:,-1].reshape(ns),rho,rho[:,0].reshape(ns)) )

        dqm = zeros(nrho.shape)
        dqp = zeros(nrho.shape)
        slopep = zeros(nrho.shape)

        dqm[:,1:] = (nrho[:,1:]-nrho[:,:-1])/dx
        dqp[:,:-1] = (nrho[:,1:]-nrho[:,:-1])/dx
        dqm[:,0] = (nrho[:,0]-rho[:,-2])/dx
        dqp[:,-1] = (rho[:,1]-nrho[:,-1])/dx

        slopep = ( (dqm*dqp) > 0).astype(int) * (2*(dqp*dqm)/(dqp+dqm))

        res2.data = nrho[:,:-2] + .5*slopep[:,:-2]*dx
#        res2.data[:,:-1] = .5*(self.dens.data[:,1:] + self.dens.data[:,:-1])
#        res2.data[:,-1] = .5*(self.dens.data[:,-1] + self.dens.data[:,0])
        res2.name = 'Rhostarp'
        res2.math_name = '$\\Sigma_\\phi^*$'
        res2.recalculate()


        return res,res1,res2
    def grad(self,q):
        res = zeros(q.shape)
        one_dim = False
        try:
            q.shape[1]
        except IndexError:
            one_dim = True


        for i in range(1,q.shape[0]-1):
            if one_dim:
                res[i] = (q[i+1]-q[i-1])/(self.r[i+1]-self.r[i-1])
            else:
                res[i,:] = (q[i+1,:]-q[i-1,:])/(self.r[i+1]-self.r[i-1])

        if one_dim:
            res[0] = (q[1]-q[0])/(self.r[1]-self.r[0])
            res[-1] = (q[-1]-q[-2])/(self.r[-1]-self.r[-2])
        else:
            res[0,:] = (q[1,:]-q[0,:])/(self.r[1]-self.r[0])
            res[-1,:] = (q[-1,:]-q[-2,:])/(self.r[-1]-self.r[-2])
        return res
    def nu(self,x):
        return self.dens.aspectratio**2 * self.dens.alpha  * x**(2*self.dens.flaringindex+.5)
    def vr_nu(self,x):
        return -1.5*self.nu(x)/x
    def scaleH(self,x):
        return self.dens.aspectratio * x**(self.dens.flaringindex + 1)
    def cs(self,x):
        return self.scaleH(x) * x**(-1.5)
    def inner_disk_sol(self,rlims):
        ind = (self.r>=rlims[0])&(self.r<=rlims[1])
        popt,pcov = curve_fit(lambda x,a: a/(3*pi*self.nu(x)),self.r[ind],self.dbar[ind])
        return popt[0],popt[0]/(3*pi*self.nu(self.r))

    def plot2d(self,q,axex=None,fig=None,logscale=False,logr=True,rlims=None,plims=None,norm=False,**kargs):
        if fig == None:
            #fig,axes=subplots(1,2)
            figsize = kargs.pop('figsize',(20,15))
            fig = figure(figsize=figsize);
            gs = gridspec.GridSpec(3,4)
            axcart = fig.add_subplot(gs[:2,:-2])
            axcyl = fig.add_subplot(gs[:2,-2:])
            axdbar = fig.add_subplot(gs[-1,:])


        if q=='dens':
            fld=self.dens
            dat = copy.copy(self.dens.data.transpose())
            dstr0 = '$\\Sigma(r)$'
            dbar = copy.copy(self.dbar)
            if norm:
                dat0 = copy.copy(self.dbar0)
        elif q=='vr':
            fld =self.vr
            dat = copy.copy(self.vr.data.transpose())
            dstr0 = '$v_r$'

            dbar = copy.copy(self.vrbar)
            if norm:
                dat0 = copy.copy(self.vr_visc)
        elif q=='vp':
            fld=self.vp
            dat = copy.copy(self.vp.data.transpose())
            dbar = copy.copy(self.vpbar)
            dstr0 = '$v_\\phi$'
            if norm:
                dat0 = pow(self.r,-1.5)

        else:
            print '%s is not a valid option' % q
            return

        if logscale:
            dstr = '$\\log_{10}$'+ dstr0
        else:
            dstr = dstr0
        if rlims == None:
            rlims = (self.dens.ymin,self.dens.ymax)
        if plims == None:
            plims = (-pi,pi)

        print rlims,plims

        rinds = (self.r<=rlims[1])&(self.r>=rlims[0])
        pinds = (self.phi<=plims[1])&(self.phi>=plims[0])

        print self.phi[pinds]
        dat = dat[:,rinds][pinds,:]

        if norm:
            if q=='vp':
                dbar = (dbar- dat0)[rinds]
                for i in range(dat.shape[0]):
                    dat[i,:] -= dat0[rinds]
            else:
                dbar =  (dbar/dat0)[rinds]
                for i in range(dat.shape[0]):
                    dat[i,:] /= dat0[rinds]


        if logscale:
            dat = log10(dat)
        r = self.r[rinds]
        lr = log(r)
        phi = self.phi[pinds]

        print r,phi
        rlims0 = copy.copy(rlims)
        if logr:
            rstr = '$\\ln r$'
            rlims = (log(rlims[0]),log(rlims[1]))
            rr,pp = meshgrid(lr,phi)
        else:
            rstr = '$r$'
            rr,pp = meshgrid(r,phi)

        print rlims

        cmap = kargs.pop('cmap','viridis')
        line2d = axcyl.pcolormesh(rr,pp,dat,cmap=cmap,shading='gouraud')
        cbar = colorbar(line2d,ax=axcyl)
        cbar.set_label(dstr,fontsize=20)
        axcyl.set_xlim(rlims)
        axcyl.set_ylim(plims)
        axcyl.set_xlabel(rstr,fontsize=20)
        axcyl.set_ylabel('$\\phi$',fontsize=20)

        axdbar.plot(r,dbar)
        axdbar.set_xlabel('$r$',fontsize=20)
        axdbar.set_ylabel(dstr0,fontsize=20)
        if norm:
            axdbar.axhline(1,color='k')
        if logr:
            axdbar.set_xscale('log')
        if logscale:
            axdbar.set_yscale('log')
        axdbar.set_xlim(rlims0)
        fld.plot(cartesian=True,ax=axcart,log=logscale)


        return fig,[axcart,axcyl,axdbar]

    def potential(self,r,phi):
        if not self.dens.rochesmoothing:
            smoothing = self.dens.thicknesssmoothing*self.scaleH(self.a)
        else:
            smoothing = self.dens.thicknesssmoothing*self.rh
        smoothing *= smoothing
        rad = r**2 + smoothing + self.a**2 - 2*r*self.a*cos(phi)
        rad = sqrt(rad)
        return -self.mp/rad
    def dp_potential(self,r,phi):
        if not self.dens.rochesmoothing:
            smoothing = self.dens.thicknesssmoothing*self.scaleH(self.a)
        else:
            smoothing = self.dens.thicknesssmoothing*self.rh
        smoothing *= smoothing
        rad = r**2 + self.a**2 - 2*r*self.a*cos(phi) + smoothing
        rad = rad**(1.5)
        return self.mp * r*self.a*sin(phi)/rad
    def dr_potential(self,r,phi):
        if not self.dens.rochesmoothing:
            smoothing = self.dens.thicknesssmoothing*self.scaleH(self.a)
        else:
            smoothing = self.dens.thicknesssmoothing*self.rh
        smoothing *= smoothing
        rad = r**2 + self.a**2 - 2*r*self.a*cos(phi) + smoothing
        rad = rad**(1.5)
        return self.mp * (r-self.a*cos(phi))/rad


    def conv_fargo_tq_1d(self,r,tq):
        tq *= -self.mp/(r**2 * self.dlr*2*pi)
        return tq
    def conv_fargo_tq_tot(self,tq):
        tq *= -self.mp/(2*pi)
        return tq

    def animate_mdot(self,irange,logspacing=True,cmap='viridis',**kargs):

        self.__init__(0,directory=self.directory)
        temp0,temp1,temp2,temp3,temp4,temp5,temp6  = self.calculate_total_mdot(logspacing)
        dens_list = [log10(temp0)]
        divm_list = [temp1]
        divm_mean_list = [temp2]
        drmx_mean_list = [temp3]
        dpmy_mean_list = [temp4]
        mx_mean_list = [temp5]
        my_mean_list = [temp6]
        t = [self.t]

        fig,axes = subplots(1,4,figsize=(15,10))

        axes[3].set_title('t = %.1f' % t[0],color='w')
        line0=axes[0].imshow(dens_list[0].transpose(),aspect='auto',origin='lower',cmap=cmap,**kargs)
        line1=axes[1].imshow(divm_list[0].transpose(),aspect='auto',origin='lower',cmap=cmap,**kargs);

        line2,=axes[2].plot(self.r,divm_mean_list[0],label='< $ \\nabla \\cdot \\dot{M} = \\dot{\\Sigma}$ >')
        line3,=axes[2].plot(self.r,drmx_mean_list[0],label='<$\\nabla_r (r \\Sigma v_r)$>')
        line4,=axes[2].plot(self.r,dpmy_mean_list[0],label='-<$\\nabla_\\phi(\\Sigma v_\\phi)$>')
        axes[2].legend(loc='best')
        axes[2].set_ylim(-.0001,.0001)
        line5,=axes[3].plot(self.r,mx_mean_list[0],label='<-$r \\Sigma v_r$>')
        line6,=axes[3].plot(self.r,my_mean_list[0],label='-<$\\Sigma v_\\phi$>')
        axes[3].legend(loc='best')

        for i in irange:
            if i > 0:
                self.__init__(i,directory=self.directory)
                t.append(self.t)
                temp0,temp1,temp2,temp3,temp4,temp5,temp6  = self.calculate_total_mdot(logspacing)
                dens_list.append(log10(temp0).transpose())
                divm_list.append(temp1.transpose())
                divm_mean_list.append(temp2)
                drmx_mean_list.append(temp3)
                dpmy_mean_list.append(temp4)
                mx_mean_list.append(temp5)
                my_mean_list.append(temp6)


        for i in range(len(t)):
            axes[2].set_title('t = %.1f' % t[i])
            line0.set_data(dens_list[i])
            line1.set_data(divm_list[i])
            line2.set_ydata(divm_mean_list[i])
            line3.set_ydata(drmx_mean_list[i])
            line4.set_ydata(dpmy_mean_list[i])
            line5.set_ydata(mx_mean_list[i])
            line6.set_ydata(my_mean_list[i])
            fig.canvas.draw()


    def plot_mdot(self,logspacing=True):
        dens,divm, divm_mean, drmx_mean, dpmy_mean, mx_mean,my_mean = self.calculate_total_mdot(logspacing)
        fig,axes = subplots(1,4,figsize=(15,10))

        axes[0].imshow(log10(dens),aspect='auto',origin='lower');
        axes[1].imshow(divm,aspect='auto',origin='lower');
        axes[3].plot(self.r,mx_mean,self.r,my_mean)
        axes[3].legend(['<-$r \\Sigma v_r$>', '-<$\\Sigma v_\\phi$>'],loc='best')
        axes[2].plot(self.r,divm_mean,self.r,drmx_mean,self.r,dpmy_mean)
        axes[2].legend(['< $ \\nabla \\cdot \\dot{M} = \\dot{\\Sigma}$ >', '-<$\\nabla_r (r \\Sigma v_r)$>', '-<$\\nabla_\\phi(\\Sigma v_\\phi)$>'],loc='best')


    def calculate_total_mdot(self,logspacing=True):
        pp,rr = meshgrid(self.phi,self.r)
        if logspacing:
            dr = diff(log(self.r))[0]
            norm = rr
        else:
            dr = diff(self.r)[0]
            norm = 1

        dp = diff(self.phi)[0]

        dens = self.dens.data
        vp = self.vpc
        vr = self.vr.data

        mx = -rr * dens*vr
        my = -dens*vp

        drmx,_ = gradient(mx,dlr,dp)
        _,dpmy = gradient(my,dlr,dp)

        drmx /= (rr*norm)
        dpmy /= rr
        divm = drmx + dpmy
        return dens,divm, divm.mean(axis=1),drmx.mean(axis=1), dpmy.mean(axis=1), mx.mean(axis=1), my.mean(axis=1)



    def summary(self):
    	fig,(axd,axm,axv) = subplots(3,1,sharex='col')
    	lined, = axd.plot(self.r,self.dbar,linewidth=3)
    	linev, = axv.plot(self.r,self.vrbar,linewidth=3)
    	linem, = axm.plot(self.r,self.mdot,linewidth=3)
    	axv.set_xlabel('$r$',fontsize=20)
    	axv.set_ylabel('$<v_r$>',fontsize=20)
    	axd.set_ylabel('$<\\Sigma>$',fontsize=20)
    	axm.set_ylabel('<$\\dot{M}$>',fontsize=20)
    #		axm.set_yscale('symlog',linthreshy=1e-7)
    #		axv.set_ylim(-.001,.001)
    	axd.set_title('t = %.1f = %.1f P = %.1f t_visc' % (self.t,self.t/(2*pi*self.a**(1.5)),self.t/self.tvisc))
    def streams(self,rlims=None,plims=None,ax=None,planet=None,noise=.1,clrbar=True,softening=False,sample=1,rasterized=False,**kargs):
        npoints= kargs.pop('npoints',10)
        draw_flag = False
        if ax == None:
            fig=figure()
            ax=fig.add_subplot(111)
            draw_flag = True

        if rlims == None:
            rlims = (self.dens.ymin,self.dens.ymax)
        if plims == None:
            plims = (-pi,pi)


        vr = .5*(self.vr.data[1:,:-1] + self.vr.data[:-1,:-1])
        vp = .5*(self.vx.data[:-1,1:] + self.vx.data[:-1,:-1])
        dens = copy.copy(self.dens.data[:-1,:-1])
        y = copy.copy(self.dens.y[:-1])
        x = copy.copy(self.dens.x[:-1])
        rinds = (y<=rlims[1])&(y>=rlims[0])
        pinds = (x<=plims[1])&(x>=plims[0])
        y = y[rinds]
        x = x[pinds]
        vr = vr[rinds,:][:,pinds]
        vp = vp[rinds,:][:,pinds]
        dens = dens[rinds,:][:,pinds]
        #vr = self.vr.data[rinds,:][:,pinds]
        #vp = self.vp.data[rinds,:][:,pinds]
        #dens = self.dens.data[rinds,:][:,pinds]
        rr,pp = meshgrid(y,x)
        lr = log(y)
        phi = x
        #rr,pp = meshgrid(self.r[rinds],self.phi[pinds])
        #lr = log(self.r[rinds])
        #phi = self.phi[pinds]

        cmap = kargs.pop('cmap','viridis')
        color = kargs.pop('color','w')

        line2d= ax.pcolormesh(log(rr[::sample,::sample]),pp[::sample,::sample],log10(dens[::sample,::sample].transpose()),cmap=cmap,shading='gouraud')
        if clrbar:
            cbar = colorbar(line2d,ax=ax)
            cbar.set_label('$\\log_{10}{\\Sigma}$',fontsize=20)
        ax.streamplot(lr,phi,vr.transpose(),vp.transpose(),color=color,**kargs)
        ax.set_xlabel('$\ln(r)$',fontsize=20)
        ax.set_ylabel('$\phi$',fontsize=20)
        ax.set_ylim(plims)
        ax.set_xlim((log(rlims[0]),log(rlims[1])))
        ax.axvline(log(self.a+2*self.rh),color='k',linewidth=3)
        ax.axvline(log(self.a-2*self.rh),color='k',linewidth=3)


        if softening:
            rh,ph = self.calculate_circle(self.rh,rlims)
            rs,ps = self.calculate_circle(self.dens.aspectratio*self.a*self.dens.thicknesssmoothing,rlims)

            ax.plot(log(rh), ph,'-r',linewidth=3,rasterized=rasterized)
            ax.plot(log(rs),ps,'--r',linewidth=3,rasterized=rasterized)

        lb,rb,sep_lines = self.horseshoe_width(noise=noise,npoints=npoints,rlims=rlims,plims=plims)

        for line in sep_lines:
            ax.plot(line[:,0],line[:,1],'-w',linewidth=2,rasterized=rasterized)

        if planet is None:
            xlbl = ['%.1f'%(exp(v)) for v in ax.get_xticks()]
            ax.set_xticklabels(xlbl)
            ax.set_xlabel('$r$',fontsize=20)
        else:
            xlbl = ['%.1f'%((exp(v)-planet)/(self.params.aspectratio*planet)) for v in ax.get_xticks()]
            ax.set_xticklabels(xlbl)
            ax.set_xlabel('$(r-a)/H$',fontsize=20)
        if draw_flag:
            fig.canvas.draw()
#        return log(self.r[rinds]),self.phi[pinds],rr,pp,vr.transpose(),vp.transpose()
        return lb,rb

    def horseshoe_width(self,noise=.1,npoints=10,rlims=None,plims=(-2,2)):
        if rlims == None:
            rlims = (self.dens.ymin,self.dens.ymax)
        if plims == None:
            plims = (-pi,pi)


        vr = .5*(self.vr.data[1:,:-1] + self.vr.data[:-1,:-1])
        vp = .5*(self.vx.data[:-1,1:] + self.vx.data[:-1,:-1])
        dens = copy.copy(self.dens.data[:-1,:-1])
        y = copy.copy(self.dens.y[:-1])
        x = copy.copy(self.dens.x[:-1])
        rinds = (y<=rlims[1])&(y>=rlims[0])
        pinds = (x<=plims[1])&(x>=plims[0])
        y = y[rinds]
        x = x[pinds]
        vr = vr[rinds,:][:,pinds]
        vp = vp[rinds,:][:,pinds]
        dens = dens[rinds,:][:,pinds]
        #rr,pp = meshgrid(y,x)
        lr = log(y)
        phi = x
        sep_lines = self.separatrix(lr,phi,vr.transpose(),vp.transpose(),noise=noise,npoints=10)

        right_b = max([max(line[:,0]) for line in sep_lines])
        left_b = min([min(line[:,0]) for line in sep_lines])

        return left_b,right_b,sep_lines
    def summary_plot(self,savefile=None):
        tstr1 = '$q = %.1e, \\alpha = %.1e, K = %.1e$'%(self.mp,self.params.alpha,self.K)
        tstr2 = '$t = %.2e = %.2e P = %.2e t_{visc,p}$'%(self.t,self.t/self.torb,self.t/self.tviscp)
        mmr2 = self.a * 2**(2./3)
        mmr3 = self.a  * 3**(2./3)

        fig,axes = subplots(1,2,figsize=(20,10))
        lb,rb=self.streams(rlims=(.6,2),plims=(-2,2),ax=axes[1])
        axes[0].plot(self.dens.y,self.dens.avg,'-k')
        axes[0].plot(self.dens.y,self.sig0,'--k')
        axes[0].axvline(mmr2,color='k',linestyle='--')
        axes[0].axvline(mmr3,color='k',linestyle='--')
        axes[0].set_xlabel('$r$',fontsize=15)
        axes[0].set_ylabel('$\\langle \\Sigma \\rangle$',fontsize=15)
        axes[0].axvline(lb+self.a,color='k',linewidth=2)
        axes[0].axvline(rb+self.a,color='k',linewidth=2)

        ax = plt.axes([.26,.65,.2,.2])
        ax.plot(self.dens.y,self.dens.avg,'-k')
        ax.plot(self.dens.y,self.sig0,'--k')
        ax.axvline(mmr2,color='k',linestyle='--')
        ax.axvline(mmr3,color='k',linestyle='--')
        ax.axvline(lb+self.a,color='k',linewidth=2)
        ax.axvline(rb+self.a,color='k',linewidth=2)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.minorticks_on()

        try:
            xt = (self.trq_y -self.a)/self.scaleH(self.a)
            left_bound = (self.dens.ymin + (self.dens.ymax-self.dens.ymin)*.0497 - self.a)/self.scaleH(self.a)
            right_bound = 10.
            ind = (xt>left_bound)&(xt<right_bound)
            ax2 = plt.axes([.26,.4,.2,.2])
            ax2.plot(xt[ind],self.trq_Lamex[ind],'-k',linewidth=2)
            ax2.plot(xt[ind],self.trq_Lamdep[ind],'-r',linewidth=2)
            ax2.plot(xt[ind],self.trq_drFw[ind],'-b',linewidth=2)
            mmr2x = (mmr2-self.a)/self.scaleH(self.a)
            mmr3x = (mmr3-self.a)/self.scaleH(self.a)
            hsL = (lb+self.a-self.a)/self.scaleH(self.a)
            hsR = (rb+self.a-self.a)/self.scaleH(self.a)
            if mmr2x > left_bound and mmr2x < right_bound:
                ax2.axvline(mmr2x,color='k',linestyle='--')
            if mmr3x > left_bound and mmr3x < right_bound:
                ax2.axvline(mmr3x,color='k',linestyle='--')
            if hsL > left_bound and hsL < right_bound:
                ax2.axvline(hsL,color='k',linewidth=2)
            if hsR > left_bound and hsR < right_bound:
                ax2.axvline(hsR,color='k',linewidth=2)
            ax2.minorticks_on()
            ax2.set_xlim(left_bound,right_bound)
            ax2.set_xlabel('$(r-r_p)/H(r_p)$',fontsize=15)
            print 'Added torque to plot'
        except AttributeError:
            print "Could't find a torque file for this output!"

        axes[0].set_title(tstr1,fontsize=15)
        axes[1].set_title(tstr2,fontsize=15)
        if mmr2 < 2:
            axes[1].axvline(log(mmr2),color='r',linestyle='-',linewidth=2)
        if mmr3 < 2:
            axes[1].axvline(log(mmr3),color='r',linestyle='-',linewidth=2)

        for axp in axes:
            axp.minorticks_on()
        fig.canvas.draw()

        if savefile is not None:
            fig.savefig(savefile)
        return fig,axes

    def calculate_circle(self,d,rlims=None):
        if rlims == None:
            rlims = (self.dens.ymin,self.dens.ymax)

        rinds = (self.r>=rlims[0])&(self.r<=rlims[1])
        newinds = (self.r[rinds]>=(self.a-d))&(self.r[rinds]<=(self.a+d))
        circle_upper = arccos((self.r[rinds][newinds]**2 + self.a**2 - d**2)/(2*self.a*self.r[rinds][newinds]))
        circle_lower=  -circle_upper

        circle_r = hstack((array([self.a-d]),self.r[rinds][newinds]))
        circle_r = hstack((circle_r,array([self.a+d])))
        circle_r = hstack((circle_r,self.r[rinds][newinds][::-1]))

        circle_phi = hstack((array([0]),circle_upper))
        circle_phi = hstack((circle_phi,array([0])))
        circle_phi = hstack((circle_phi,circle_lower[::-1]))
        return circle_r, circle_phi

    def stagnation(self,lr,phi,vr,vp):
        cr = contour(lr,phi,vr,levels=(0,),colors='w')
        cp = contour(lr,phi,vp,levels=(0,),colors='w')

        stag_points = []
        for pathr in cr.collections[0].get_paths():
            vertr = pathr.vertices
            segsr =[ vstack( (vertr[i,:],vertr[i+1,:])) for i in range(vertr.shape[0]-1)]
            for pathp in cp.collections[0].get_paths():
                vertp = pathp.vertices
                segsp =[ vstack( (vertp[i,:],vertp[i+1,:])) for i in range(vertp.shape[0]-1)]
                for sr in segsr:
                    for sp in segsp:
                        if intersect(sr,sp):
                            stag_points.append(get_intersection(sr,sp))

        return stag_points

    def separatrix(self,lr,phi,vr,vp,noise=0.1,npoints=10):
        ivp = interp2d(lr,phi,vp)
        ivr = interp2d(lr,phi,vr)
        dlr = lr[1]-lr[0]
        dp = phi[1]-phi[0]

        plt.figure()
        stag_points = self.stagnation(lr,phi,vr,vp)
        plt.close()
        uvals=np.array([x[0] for x in stag_points])
        inds = np.argsort(uvals)
#        for x in stag_points:
#            plot(x[0],x[1],'ow',markersize=20)
        sL = stag_points[inds[0]]
        sR = stag_points[inds[-1]]
        lines=[]

        tempx = noise*np.sqrt(dlr**2 + dp**2)

        for i in range(npoints*2):
            theta = 2*np.pi *  i/(float(2*npoints))
            u0 = sL[0] + tempx * np.cos(theta)
            p0 = sL[1] + tempx*np.sin(theta)

#            u0 = stag_points[0][0] + dlr*noise*2*(-.5 + rand())
#            p0 = stag_points[0][1] + dp*noise*2*(-.5+rand())
            lines.append(self.get_stream(u0,p0,ivp,ivr,lr,phi,False))
            lines.append(self.get_stream(u0,p0,ivp,ivr,lr,phi,True))
            u0 = sR[0] + tempx * np.cos(theta)
            p0 = sR[1] + tempx*np.sin(theta)

#            u0 = stag_points[0][0] + dlr*noise*2*(-.5 + rand())
#            p0 = stag_points[0][1] + dp*noise*2*(-.5+rand())
            lines.append(self.get_stream(u0,p0,ivp,ivr,lr,phi,False))
            lines.append(self.get_stream(u0,p0,ivp,ivr,lr,phi,True))


        return lines
    def get_stream(self,u0,p0,ivp,ivr,lr,phi,reverse=False):
        dlr = lr[1]-lr[0]
        dp = phi[1]-phi[0]

        lr_min,lr_max = (lr.min(),lr.max())
        phi_min,phi_max = (phi.min(),phi.max())

        l = np.min((dlr,dp))
        nmax=1000

        res = np.array([u0,p0])
        uj = u0
        pj=  p0
        breakflag = False
        for j in range(nmax):
            dx, dy = self.euler_step(uj,pj,l,ivp,ivr,reverse)
            uj += dx
            pj += dy
            if uj >= lr_max:
                uj = lr_max
#                print 'Stream lnr > lnr_max'
                breakflag = True
            if uj <= lr_min:
 #               print 'Stream lnr < lnr_min'
                uj = lr_min
                breakflag = True
            if pj >= phi_max:
  #              print 'Stream phi > phi_max'
                pj = phi_max
                breakflag = True
            if pj <= phi_min:
   #             print 'Stream phi < phi_min'
                pj = phi_min
                breakflag = True
            if breakflag:
                break
            else:
               res=np.vstack( (res,np.array([uj,pj])))
    #    print 'Iterated past max iterations of %d' % nmax
        return res

    def euler_step(self,u0,p0,l,ivp,ivr,reverse=False):
        safety_fac = self.safety_fac
        h = 1.0
        if reverse:
            h = -1.0

        vp_val = ivp(u0,p0)[0]
        vr_val = ivr(u0,p0)[0]

        h *= safety_fac*l/np.sqrt( vp_val**2 + vr_val**2)
        return h*vr_val, h*vp_val



    def animate(self,irange):
    	fig,(axd,axm,axv) = subplots(3,1,sharex='col')
    	lined, = axd.plot(self.r,self.dbar,linewidth=3)
    	linev, = axv.plot(self.r,self.vrbar,linewidth=3)
    	linem, = axm.plot(self.r,self.mdot,linewidth=3)

        s0 = Sim(0,directory=self.directory)
        axd.plot(s0.r,s0.dbar,'--k',linewidth=3)
    	axv.plot(s0.r,s0.vrbar,'--k',linewidth=3)
    	axm.plot(s0.r,s0.mdot,'--k',linewidth=3)
        axv.set_ylim(-.05,.05)

    	axv.set_xlabel('$r$',fontsize=20)
    	axv.set_ylabel('<$v_r$>',fontsize=20)
    	axd.set_ylabel('<$\\Sigma$>',fontsize=20)
    	axm.set_ylabel('<$\\dot{M}$>',fontsize=20)

    	dbar = zeros((len(self.dbar),len(irange)))
    	vrbar = zeros((len(self.dbar),len(irange)))
    	mdot = zeros(dbar.shape)
    	t = zeros(len(irange))
    	a = zeros(len(irange))

    	t[0] = self.t
    	a[0] = self.a
    	dbar[:,0] = self.dbar
    	vrbar[:,0] = self.vrbar
    	mdot[:,0] = self.mdot
    	for i,j in enumerate(irange[1:],start=1):
            self.__init__(j,directory=self.directory)
            t[i] = self.t
            dbar[:,i] = self.dbar
            vrbar[:,i] = self.vrbar
            mdot[:,i] = self.mdot
            a[i] = self.a

#    	axv.set_ylim(-.001,.001)
#    	axm.set_ylim((mdot.min(),mdot.max()))
#    	axm.set_yscale('symlog',linthreshy=1e-8)
    	for i in range(len(t)):
    		lined.set_ydata(dbar[:,i])
    		linev.set_ydata(vrbar[:,i])
    		linem.set_ydata(mdot[:,i])
    		axd.set_title('t = %.1f = %.1f P' % (t[i],t[i]/(2*pi*a[i]**(1.5))),color='w')
    		fig.canvas.draw()


    def steady_state_calc(self):
        self.mdot0_ss = 3*pi*self.nu0*(self.dens.ymax*self.dens.sigmaout-self.dens.ymin*self.dens.sigmain)/(np.sqrt(self.dens.ymax)-np.sqrt(self.dens.ymin))
        self.lam0_ss = 2*pi*self.dens.sigmain*self.dens.ymax + 2*self.mdot0_ss*(np.sqrt(self.r)-np.sqrt(self.dens.ymin))/(3*self.nu0)
        self.vr0_ss = -self.mdot0_ss/self.lam0_ss
        self.dbar0_ss = self.lam0_ss/(2*pi*self.r)


    def refine_grid(self,levels=3,logspacing=True,fname='domain_y.dat',save=False,savedir='./'):
        if levels > 3:
            print 'Too many levels, maximum is 3'
            return
        rad = loadtxt(fname)
        hs_zone = lambda r: abs(r-self.a) <= self.rh*2
        rh_zone = lambda r: abs(r-self.a) <= self.rh
        eps = self.dens.thicknesssmoothing * self.dens.aspectratio * self.a
        soft_zone = lambda r: abs(r-self.a) <= eps

        print 'Planet at %f\nHs zone at %f\nHill zone at %f\nSoft zone at %f' %(self.a,2*self.rh,self.rh,eps)

        if eps > self.rh:
            if eps > 2*self.rh:
                ind_func = [soft_zone,hs_zone,rh_zone]
            else:
                ind_func = [hs_zone,soft_zone,rh_zone]
            print 'smoothing zone larger than hill zone'
        else:
            ind_func = [hs_zone,rh_zone,soft_zone]

        for j in range(levels):
            ind = ind_func[j](rad)
            ri = rad[ind][0]
            ro = rad[ind][-1]

            if logspacing:
                spacing = self.dlr * 2.**(-(j+1))
                rad1 = np.exp(np.linspace(np.log(ri),np.log(ro),np.log(ro/ri)/spacing+1))
            else:
                spacing = self.dr[0] * 2.**(-(j+1))
                rad1 = np.linspace(ri,ro,(ro-ri)/spacing+1)

            leftr = rad[ rad < ri]
            rightr = rad[ rad > ro]

            rad = np.hstack( (leftr,rad1,rightr) )

        print 'Went from %d points to %d points' % (self.dens.ny,len(rad)-7)

        if save:
            call(['cp',fname,fname+'.backup'])
            with open(savedir+fname,'w') as f:
                f.write('\n'.join(rad.astype(str)))

        return rad
    def refine_phi(self,j=2,fname='domain_x.dat',save=False,savedir='./'):
        x = np.loadtxt(fname)
        dx = diff(x)[0]
        x = np.linspace(x[0],x[-1],(x[-1]-x[0])/(dx/j)+1)


        print 'Went from %d points to %d points' % (self.dens.ny,len(x)-1)

        if save:
            call(['cp',fname,fname+'.backup'])
            with open(savedir+fname,'w') as f:
                f.write('\n'.join(x.astype(str)))

        return x
    def linear_torque(self,r,eps=.2,c=2./3):
        norm = eps*self.mp**2/(2*r)
        ha = self.a*self.dens.aspectratio
        x = (r-self.a)/ha

        sgn = (x>=0).astype(int)
        sgn[ sgn == 0 ] = -1

        indR = x >=c
        indL = x <= -c
        res = zeros(r.shape)
        res[indR] = (self.a/abs(r[indR]-self.a))**4
        res[indL] = (r[indL]/(r[indL]-self.a))**4

        return res * norm * sgn

    def load_torque(self,dt=6.28):
        execdir = self.fargodir+'utils/python/'
        execname = 'run_step'

        lines='%d\n' % self.params.nx
        lines +='%d\n' % self.params.ny
        lines +='%d\n' %  self.params.nz
        lines +='%f\n' %  self.params.alpha
        lines +='%f\n' %  self.mp
        lines +='%f\n' %  self.a
        lines +='%f\n' %  self.params.aspectratio
        lines +='%f\n' %  self.params.flaringindex
        lines +='%f\n' %  self.params.mdot
        lines +='%f\n' %  self.params.thicknesssmoothing
        lines +='%d\n' % (1-self.params.nz)

        with open(self.directory + 'param_file.txt','w') as f:
            f.write(lines)

        dirval = self.directory
        if dirval == '':
            dirval = './'

        callstr = [execdir + execname,'%d'%self.n,'%f'%dt,dirval]
        print callstr
        #call(callstr)

        self.torque = Torque()
        return


    def compute_fft(self):
        import numpy.fft as ft
        dhat = ft.rfft(self.dens.data,axis=1)/(self.dens.nx-1)
        vrhat = ft.rfft(self.vr.data,axis=1)/(self.dens.nx-1)
        vphat = ft.rfft(self.vp.data,axis=1)/(self.dens.nx-1)
        p_dhat = trapz(dhat*conj(dhat)*self.dr[:,np.newaxis],axis=0)
        p_vrhat = trapz(vrhat*conj(vrhat)*self.dr[:,np.newaxis],axis=0)
        p_vphat = trapz(vphat*conj(vphat)*self.dr[:,np.newaxis],axis=0)
        return dhat,vrhat,vphat,p_dhat,p_vrhat,p_vphat

def ccw(a,b,c):
    return (c[1]-a[1])*(b[0]-a[0]) > (b[1]-a[1])*(c[0]-a[0])

def intersect(segment_1, segment_2):
    a1 = segment_1[0,:]
    a2 = segment_1[1,:]
    b1 = segment_2[0,:]
    b2 = segment_2[1,:]
    return ccw(a1,b1,b2) != ccw(a2,b1,b2) and ccw(a1,a2,b1) != ccw(a1,a2,b2)
def get_intersection(segment_1, segment_2):
    if intersect(segment_1,segment_2):
#        print 'Calculating Intersection'
        m1 = segment_1[1,:]-segment_1[0,:]
        m1 = m1[1]/m1[0]
        m2 = segment_2[1,:] - segment_2[0,:]
        m2 = m2[1]/m2[0]
        lhs = np.array([[-m1,1.],[-m2,1.]])
        rhs = np.array([[ segment_1[0,1] - m1*segment_1[0,0]],[segment_2[0,1]-m2*segment_2[0,0]]])
        return dot( np.linalg.inv(lhs),rhs).reshape(2)
    else:
#        print 'Lines do not intersect'
        return False
def boundary_progress(i,fig=None,axes=None,j=0,alpha=0.01,h=0.05,flaringindex=0):
    config_flag = True
    if fig is None:
        fig,axes=subplots(3,1,sharex=True)
        config_flag = False
    try:
        rho = fromfile('gasdens{0:d}_0.dat'.format(i))
        vr = fromfile('gasvy{0:d}_0.dat'.format(i))
        ym = loadtxt('domain_y.dat')
        xm = loadtxt('domain_x.dat')
        ymed = (ym[1:] + ym[:-1])/2
        xmed = (xm[1:] + xm[:-1])/2
        ny = len(ymed)-6
        nx = len(xmed)
        vr = vr.reshape(ny+6,nx)
        rho = rho.reshape(ny+6,nx)
    except IOError:
        rho = fromfile('gasdens{0:d}.dat'.format(i))
        vr = fromfile('gasvy{0:d}.dat'.format(i))
        ym = loadtxt('domain_y.dat')
        xm = loadtxt('domain_x.dat')
        ymed = (ym[1:] + ym[:-1])/2
        xmed = (xm[1:] + xm[:-1])/2
        ny = len(ymed)-6
        nx = len(xmed)
        rho = rho.reshape(ny,nx)
        vr = vr.reshape(ny,nx)
#    mdot = 2*pi*ymed * (rho*vr).mean(axis=1)
    #mdot = 2*ymed[:-1]*pi*(rho[:-1,:]*((vr[:1,:]+vr[:-1,:])/2)).mean(axis=1)
    rho = rho.mean(axis=1)
    vr = vr.mean(axis=1)
    nu = alpha*h*h*ymed**(2*flaringindex+0.5)
    fj = 3*pi*nu*rho*sqrt(ymed)

    axes[0].plot(ymed,rho,'.-')
    axes[1].plot(ymed,fj,'.-')
    axes[2].plot(ym[:-1],vr,'.-')

    axes[0].set_ylabel('$\\Sigma$',fontsize=20)
    axes[1].set_ylabel('$F_\\nu$',fontsize=20)
    axes[2].set_ylabel('$v_r$',fontsize=20)


    axes[2].set_xlabel('$r$',fontsize=20)
#    axes[2].plot(ymed,mdot)

    if not config_flag:
        for ax in axes:
            ax.minorticks_on()
#            ax.set_xscale('log')
#        axes[0].set_yscale('log')
        subplots_adjust(hspace=0)
    return fig,axes,ymed,ym,rho,fj,vr



def time_avg_mdot(irange):

    ym = loadtxt('domain_y.dat')[3:-3]
    ymed = .5*(ym[1:] + ym[:-1])
    ymed = ymed[:-1]
    fac = 2*pi*ymed
    facm = 2*pi*ym[1:-1]
    res = zeros((len(ymed),len(irange)))
    res_rho = zeros((len(ymed),len(irange)))
    for i,t in enumerate(irange):
        rho,momy,_,_,_,momys = load_single_time_data(t)
        res[:,i] = (momys).mean(axis=1)
        res_rho[:,i] = rho.mean(axis=1)

    return -facm*res.mean(axis=1), -fac*res.std(axis=1),res_rho.mean(axis=1)*fac,res_rho.std(axis=1)*fac,ymed



def load_single_time_data(i):
    vx,vy,rho = load_single_time(i)

    vxc = 0.5*(vx.data[:,:-1] + vx.data[:,1:])
    vyc = 0.5*(vy.data[:-1,:] + vy.data[1:,:])

    momy = vyc*rho.data[:-1,:]
    momx = vxc*rho.data[:,:-1]

    momys = .5*(rho.data[1:,:] + rho.data[:-1,:])*vy.data[1:,:]



    return rho.data[:-1,:],momy,vxc,vyc,rho.ymed[:-1],momys

def vortencity(rho,vx,vy):
    Tx,Rx = meshgrid(vx.x,vx.y)
    Ty,Ry = meshgrid(vy.x,vy.y)
    rvx = Rx*(vx.data)

    curldata = (( rvx[1:,1:]-rvx[:-1,1:])/(Rx[1:,1:]-Rx[:-1,1:]) -
            (vy.data[1:,1:] - vy.data[1:,:-1])/(Ty[1:,1:]-Ty[1:,:-1]))

    curl = copy.deepcopy(vx)
    curl.nx = curl.nx-1
    curl.ny = curl.ny-1
    curl.x = vx.x[:-1]
    curl.y = vy.y[:-1]

    rho_corner = .25*(rho.data[1:,1:] + rho.data[:-1,:-1] + rho.data[1:,:-1] + rho.data[:-1,1:])

    T,R = meshgrid(curl.x,curl.y)
    curl.data = (curldata/R + 2*rho.omegaframe)/rho_corner
    return curl

def load_single_time(i):

    rho = Field('gasdens{0:d}.dat'.format(i))
    vx = Field('gasvx{0:d}.dat'.format(i),staggered='x')
    vy = Field('gasvy{0:d}.dat'.format(i),staggered='y')


    return vx,vy,rho




def load_flux(i):

    ym = loadtxt('domain_y.dat')[3:-3]
    ymed = .5*(ym[1:] + ym[:-1])

    xm = loadtxt('domain_x.dat')
    xmed = .5*(xm[1:] + xm[:-1])
    dx = diff(xm)
    dy = diff(ym)

    ny = len(ymed)
    nx  = len(xmed)


    rho = fromfile('gasdens{0:d}.dat'.format(i)).reshape(ny,nx)
    vy = fromfile('gasvy{0:d}.dat'.format(i)).reshape(ny,nx)

    dqm = zeros(rho.shape)
    dqp = zeros(rho.shape)
    slope = zeros(rho.shape)
    mdot = zeros(rho.shape)

    surfy = outer(ym,dx)
    dxx,dyy = meshgrid(dx,dy)
    dqm = zeros(rho.shape)
    dqp = zeros(rho.shape)
    slope1 = zeros(rho.shape)
    mdot1=zeros(rho.shape)
    surfy1 = outer(ym,dx)
    dqm[1:-1] = (rho[1:-1,:]-rho[:-2,:])/dyy[1:-1,:]
    dqp[1:-1] = (rho[2:,:]-rho[1:-1,:])/dyy[2:,:]
    slope1[1:-1,:] = ( (dqm*dqp)[1:-1,:] > 0).astype(int) * (2*(dqp*dqm)/(dqp+dqm))[1:-1,:]
    mdot[1:-1,:] = (vy[1:-1,:]>0).astype(int) * vy[1:-1,:]*(rho[:-2,:]+.5*slope[:-2,:]*dyy[:-2,:])*surfy1[1:-2,:]
    mdot[1:-1,:] += (vy[1:-1,:]<=0).astype(int) * vy[1:-1,:]*(rho[1:-1,:]-.5*slope[1:-1,:]*dyy[1:-1,:])*surfy1[1:-2,:]

    return ym[1:-2],-mdot.sum(axis=1)[1:-1]

def calc_mass_flux(i,cfl=0.5,h=0.05,dt=None,fac=1,plot_flag=False):


    ym = loadtxt('domain_y.dat')[3:-3]
    ymed = .5*(ym[1:] + ym[:-1])

    xm = loadtxt('domain_x.dat')
    xmed = .5*(xm[1:] + xm[:-1])
    dx = diff(xm)
    dy = diff(ym)

    ny = len(ymed)
    nx  = len(xmed)

    rho = fromfile('gasdens{0:d}.dat'.format(i)).reshape(ny,nx)
    vx = fromfile('gasvx{0:d}.dat'.format(i)).reshape(ny,nx)
    vy = fromfile('gasvy{0:d}.dat'.format(i)).reshape(ny,nx)
    rho1 = fromfile('gasdens{0:d}.dat'.format(i+1)).reshape(ny,nx)
    vy1 = fromfile('gasvy{0:d}.dat'.format(i+1)).reshape(ny,nx)



    flux_simple = zeros(rho.shape)
    rhonew_simple = copy.copy(rho)



    if dt is None:
        dt = cfl * min( dy/(h/sqrt(ymed)))

    print 'Using dt = %f' % dt

    rhostar = zeros(rho.shape)
    slope = zeros(rho.shape)

    zone_size_y = lambda j: ym[j+1]-ym[j]

    for i,yc in enumerate(ymed[1:-1],start=1):
        for j,xc in enumerate(xmed):
            dqm = (rho[i,j]-rho[i-1,j]) / zone_size_y(i)
            dqp = (rho[i+1,j]-rho[i,j]) / zone_size_y(i+1)
            if dqm*dqp <= 0:
                slope[i,j] = 0
            else:
                slope[i,j] = 2*dqp*dqm/(dqp+dqm)


 #   dqm = zeros(rho.shape)
 #   dqp = zeros(rho.shape)
 #   slope = zeros(rho.shape)
 #   mdot = zeros(rho.shape)

 #   surfy = outer(dy,dx)
 #   dxx,dyy = meshgrid(dx,dy)
 #   vol = outer(.5*(ym[1:]**2-ym[:-1]**2),dx)

 #   dqm[1:-1,:] =(rho[1:-1,:]  - rho[:-2,:])/dy[1:-1,:]
 #   dqp[1:-1,:] =(rho[2:,:]  - rho[1:-1,:])/dy[2:,:]
 #   ind = sign(dqm*dqp)
 #   slope = (ind>0).astype(int) * 2*dqm*dqp/(dqm+dqp)
 #   mdot[1:-1,:] = (vy>0).astype(int)*(rho[:-2,:]+.5*slope[:-2,:]*dy[:-2,:]) + (vy<=0).astype(int)*(rho[1:,:]-.5*slope[1:,:]*dy[1:,:])
 #   mdot[1:-1,:] *= surfy[1:-1,:]*vy[1:-1,:]



    for i,yc in enumerate(ymed[1:-1],start=1):
        for j,xc in enumerate(xmed):
            if vy[i,j] > 0:
                rhostar[i,j] = rho[i-1,j] + .5*slope[i-1,j]*( zone_size_y(i-1) - vy[i,j]*dt)
            else:
                rhostar[i,j] = rho[i,j] - .5*slope[i,j]*( zone_size_y(i) + vy[i,j]*dt)


    rhonew = copy.copy(rho)
    flux = zeros(rho.shape)
    vol = (2*pi/ny)*.5*(ym[1:]**2 - ym[:-1]**2)
    surfy = (2*pi/ny)*ym
    inv = 1./vol

    for i,yc in enumerate(ymed[1:-1],start=1):
        for j,xc in enumerate(xmed):
            flux[i,j] = -(vy[i,j]*rhostar[i,j]*surfy[i]-vy[i+1,j]*rhostar[i+1,j]*surfy[i+1])
            rhonew[i,j] -= dt*inv[i]*flux[i,j]

    for i,yc in enumerate(ymed[1:-1],start=1):
        for j,xc in enumerate(xmed):
            flux_simple[i,j] = surfy[i+1]*.5*(rho[i,j] + rho[i+1,j])*vy[i+1,j] - surfy[i]*.5*(rho[i-1,j]+rho[i,j])*vy[i,j]
            rhonew_simple[i,j] += - dt*inv[i] * flux_simple[i,j]


    fluxSp = zeros(rho.shape)
    fluxSm = zeros(rho.shape)
    fluxp = zeros(rho.shape)
    fluxm = zeros(rho.shape)
    for i,yc in enumerate(ymed[1:-1],start=1):
        for j,xc in enumerate(xmed):
            fluxSp[i,j] = surfy[i+1]*.5*(rho[i,j] + rho[i+1,j])*vy[i+1,j]
            fluxSm[i,j] = surfy[i]*.5*(rho[i-1,j]+rho[i,j])*vy[i,j]
            fluxp[i,j] = vy[i+1,j]*rhostar[i+1,j]*surfy[i+1]
            fluxm[i,j] = vy[i,j]*rhostar[i,j]*surfy[i]


    indt = (vy>0).astype(int)
    indf = (~(vy>0)).astype(int)

    pi_p = zeros(rho.shape)
    pi_m = zeros(rho.shape)
    fluxup = zeros(rho.shape)
    for i,yc in enumerate(ymed[1:-1],start=1):
        for j,xc in enumerate(xmed):
            pi_p[i,j] = vy[i,j] * rho[i,j]*surfy[i]
            pi_m[i,j] = vy[i,j] * rho[i-1,j]*surfy[i]
            ap = slope[i-1,j]
            am = slope[i,j]
            if vy[i,j] > 0:
                fluxup[i,j] = vy[i,j]*(rho[i-1,j] + .5*ap*zone_size_y(i-1) ) * surfy[i]
            else:
                fluxup[i,j] = vy[i,j]* (rho[i,j] - .5*am*zone_size_y(i) )*surfy[i]

    if plot_flag:
        fig,axes=subplots(5,1,sharex=True)
        axes[0].set_title('$\\Delta t = %f$'%dt,fontsize=20)


        axes[4].plot(ym[:-1],vy.mean(axis=1),'-b')
        axes[4].set_ylabel('$v_r$',fontsize=20)

        axes[3].plot(ym[1:][1:-1],fluxp.mean(axis=1)[1:-1],'-b')
        axes[3].plot(ym[:-1][1:-1],fluxm.mean(axis=1)[1:-1],'--b')
        axes[3].plot(ym[1:][1:-1],fluxSp.mean(axis=1)[1:-1],'-r')
        axes[3].plot(ym[:-1][1:-1],fluxSm.mean(axis=1)[1:-1],'--r')
    #    axes[3].plot(ymed[:-1][1:-1],(ymed[:-1]*(2*pi/nx)* (.5*(vy[1:,:]+vy[:-1,:])*rho[:-1,:]).mean(axis=1))[1:-1],'--r')
    #    axes[3].plot(ym[1:-2],pi_p.mean(axis=1)[1:-1],'-k')
    #    axes[3].plot(ym[1:-2],pi_m.mean(axis=1)[1:-1],'--k')
        axes[3].plot(ym[1:-2],fluxup.mean(axis=1)[1:-1],'-m')
        axes[3].set_ylabel('$F$',fontsize=20)

        axes[2].plot(ymed[1:-1], 1.e-16 +abs( ( (rhonew-rho1)/rho1 ).mean(axis=1))[1:-1],'-b')
        axes[2].plot(ymed[1:-1], 1.e-16 +abs( ( (rhonew_simple-rho1)/rho1 ).mean(axis=1))[1:-1],'-r')
        axes[2].set_yscale('log')
        axes[2].set_ylabel('Relative Error',fontsize=15)

        axes[0].plot(ymed[1:-1],(rhonew-rho).mean(axis=1)[1:-1],'-b')
        axes[0].plot(ymed[1:-1],(rhonew_simple-rho).mean(axis=1)[1:-1],'-r')
        axes[0].plot(ymed[1:-1],(rho1-rho).mean(axis=1)[1:-1]*fac,'-g')
        axes[0].set_ylabel('$\\Delta\\Sigma$',fontsize=20)
        axes[0].legend(['Mock step', '$\\Sigma^*=\\Sigma_{avg}$','FARGO'],loc='upper right')


        axes[1].plot(ymed[1:-1],flux.mean(axis=1)[1:-1],'-b')
        axes[1].plot(ymed[1:-1],flux_simple.mean(axis=1)[1:-1],'-r')
        axes[1].set_ylabel('$\\Delta F$',fontsize=20)

        axes[-1].set_xlabel('$r$',fontsize=20)
        for ax in axes:
            ax.minorticks_on()

    return ym[:-1][1:-1],-fluxup.sum(axis=1)[1:-1]

def check_mdot():
    t = loadtxt('mass.dat')
    rho = fromfile('mass_1d_Y_raw.dat')
    momy = fromfile('momy_1d_Y_raw.dat')

    ym = loadtxt('domain_y.dat')[3:-3]
    xm = loadtxt('domain_x.dat')

    ymed = .5*(ym[:-1]+ym[1:])
    xmed = .5*(xm[:-1] + xm[1:])

    surf = .5*(ym[1:]**2 - ym[:-1]**2)
    vol = surf * diff(xm)[0] # Uniform gird in phi

    dy = diff(ym)
    dly = dy/ymed

    ny = len(ymed)
    nx = len(xmed)
    nt = len(rho)/ny

    if len(rho)/ny != len(momy)/ny or len(rho)/ny != len(t[:,0]) or len(momy)/ny != len(t[:,0]):
        print "File's loaded at different times! Try again."
        return

    rho = rho.reshape(nt,ny)
    momy = momy.reshape(nt,ny)
    t = t[:,0]


    for i in range(rho.shape[0]):
        rho[i,:] /= vol
        momy[i,:] /= vol
        momy[i,:] *= -2*pi*ymed
    rho /= nx

    momy /= nx



    return ymed,rho,momy,t


def check_mdot_single(i):
    rho = fromfile('gasdens{0:d}.dat'.format(i))
    vy = fromfile('gasvy{0:d}.dat'.format(i))

    rho0 = fromfile('gasdens{0:d}.dat'.format(i-1))

    ym = loadtxt('domain_y.dat')[3:-3]
    xm = loadtxt('domain_x.dat')

    ymed = .5*(ym[:-1]+ym[1:])
    xmed = .5*(xm[:-1] + xm[1:])

    surf = .5*(ym[1:]**2 - ym[:-1]**2)
    vol = surf * diff(xm)[0] # Uniform gird in phi

    dy = diff(ym)
    dly = dy/ymed

    ny = len(ymed)
    nx = len(xmed)

    rho = rho.reshape(ny,nx)
    rho0 = rho0.reshape(ny,nx)
    vy = vy.reshape(ny,nx)


 #   dbar = rho.mean(axis=1)

#    momy = -2*pi*ymed[:-1]*(rho[:-1,:]*.5*(vy[1:,:] + vy[:-1,:])).mean(axis=1)


#    momys = -2*pi*ym[1:-1]*(  .5*(rho[1:,:] + rho[:-1,:])*vy[1:,:] ).mean(axis=1)


 #   rho_1 = interp1d(ymed,rho,axis=0,kind='slinear')
 #   rho_2 = interp1d(ymed,rho,axis=0,kind='quadratic')
    rho_3 = interp1d(ymed,rho,axis=0,kind='cubic')


#    momy1 =-2*pi*ym[1:-1]* (rho_1(ym[1:-1])*vy[1:,:]).mean(axis=1)
#    momy2 =-2*pi*ym[1:-1]* (rho_2(ym[1:-1])*vy[1:,:]).mean(axis=1)
    momy3 =-2*pi*ym[1:-1]* (rho_3(ym[1:-1])*vy[1:,:]).mean(axis=1)


    params = Parameters()

    dt = params.dt * params.ninterm

    dlamdt = 2*pi*ymed*(rho-rho0).mean(axis=1)/dt


    figure()
    plot(ymed,dlamdt*dy)
    figure();
    plot(ym[1:-2],diff(momy3))



    return ym[1:-1],momy3



def grab_axis(i):
    return figure(i).get_axes()[0]

def problem_size(nx,ny,nz=1):
    return 168. * nx*ny*nz/(1e9)
def get_num_cells(Nh,h,ri,ro):
    ny = int( Nh * log(ro/ri)/h)
    nx = int(2*pi/h * Nh)
    return nx,ny
def get_optimal_num_points(h,ri,ro,totsize=12,nperh=None, threed= False):
    psize = 2*pi*log(ro/ri)/(h*h)
    if nperh is not None:
        Nh = nperh
    else:
        Nh = int(sqrt( totsize * 1e9 /(168*psize) ))
        totsize *= .9
    nx,ny = get_num_cells(Nh,h,ri,ro)
    psize = problem_size(nx,ny)
    print 'For a total memory of %.1f GB, total snapshot size of %.3f GB, you can have nx=%d, ny=%d, Nh=%d' % (psize,4*8*nx*ny/(1e9),nx,ny,Nh)
    return nx,ny

def torque_fit(sims=None,dirs=None,nt=None):
    if sims is None:
        sims=[]
        if dirs is not None and nt is not None:
            for n,d in dirs:
                sims.append(Sim(n,directory=d))
        else:
            print 'Need directory list and snapshot time.'
            return

    torq = lambda x,mu,p,eps: eps*exp( -( log(abs(x)) - mu)**2/(2*p**2))/(sqrt(2*pi)*p)


    lamf_list=[]
    lamex_list=[]
    x_list=[]
    params=zeros((6,len(sims)))
    eparams=zeros((6,len(sims)))
    a_list=[]
    for i,s in enumerate(sims):
        lamex = s.lambda_ex.avg/(2*pi*s.lambda_ex.y*s.rhostar.avg)

        x = (s.lambda_ex.y-s.a)/(s.params.aspectratio*s.a)

        indR = (x>=2./3)&(x<=6)
        indL = (x<=-2./3)&(x>=-6)
        indm = (x>-2./3)&(x<2./3)
        indfL = (x<-6)
        indfR = (x>6)
        pL,eL = curve_fit(torq,x[indL],lamex[indL],p0=(log(2),.3,.003))
        eL = sqrt(diag(eL))
        pR,eR = curve_fit(torq,x[indR],lamex[indR],p0=(log(2),.3,.003))
        eR = sqrt(diag(eR))
        lamf = hstack( (zeros(x[indfL].shape),torq(x[indL],*pL),zeros(x[indm].shape),torq(x[indR],*pR),zeros(x[indfR].shape)))
        print 'Fit for sim %d'%i
        print pL,pR
        x_list.append(x)
        lamf_list.append(lamf)
        lamex_list.append(lamex)
        params[:3,i] = pL
        params[-3:,i] = pR
        eparams[:3,i] = eL
        eparams[-3:,i] = eR
        a_list.append(s.params.alpha)


    return x_list,lamf_list,lamex_list,params,eparams,a_list

def split_conv(f1,f2):
    res = real( f1.ft*conj(f2.ft) )
    res[:,1:] *= 2
    return res
def split_power(f1,f2):
    res = split_conv(f1,f2)
    res *= conj(res)
    return res.sum(axis=0)

def convolve_field(f1, f2,m):

    if m == 0:
        res  = zeros(f1.ft[:,0].shape)
        for n in range(1,f1.ft.shape[1]):
            res += 2*real(f1.ft[:,n]*conj(f2.ft[:,n]))
        return res
    else:
        res  = zeros(f1.ft[:,0].shape)
        for n in range(m):
            res += 2*real(f1.ft[:,n]*f2.ft[:,m-n])
        for n in range(m,n):
            res += 2*real(f1.ft[:,n]*f2.ft[:,n-m])
        return res
def avg(q,axis=0):
    try:
        res = trapz(q,axis=axis)/(q.shape[axis])
        return res
    except:
        return q.mean(axis=axis)



def read_dat_file(fname,ny):
    dat = fromfile(fname)
    ymed = dat[:ny]
    dat = dat[ny:]

    Lt = dat[:ny]
    dat = dat[ny:]

    Ld = dat[:ny]
    dat = dat[ny:]

    Lw = dat[:ny]
    dat = dat[ny:]

    drFt = dat[:ny]
    dat = dat[ny:]
    drFd = dat[:ny]
    dat = dat[ny:]
    drFw = dat[:ny]
    dat = dat[ny:]
    lex = dat[:ny]
    dat = dat[ny:]
    lind = dat[:ny]
    dat = dat[ny:]
    ldep = dat[:ny]
    dat = dat[ny:]
    mdt = dat[:ny]
    dat = dat[ny:]
    mdd = dat[:ny]
    dat = dat[ny:]
    return ymed,Lt,Ld,Lw,drFt,drFd,drFw,lex,lind,ldep,mdt,mdd




class Torque():

    def __init__(self,directory='temp_files/'):

#        dat = np.loadtxt('param_file.txt')
#        nx = dat[0]
#        ny = dat[1]
#        self.nx = nx
#        self.ny = ny
#        self.nz = dat[2]
#        self.alpha = dat[3]
#        self.mp = dat[4]
#        self.a = dat[5]
#        self.h = dat[6]
#      #  self.omf = dat[6]
#        self.flaringindex = dat[7]
#        self.soft = dat[9]
#        self.indirectterm = bool(dat[10])
#      #  self.dt = dat[11]
        dat = np.fromfile('torque%d.dat'%self.n)
        attrs = ['y','dbar','Lt','Ltn','Ld','Ldn','Lw','Lwn',
                'drFt','drFd','drFw',
                'Lamex','indLam','Lamdep',
                'dtLt','dtLd','dtLw','mdot']

        for i,a in enumerate(attrs):
            setattr(self,'trq_'+a,dat[i*ny:(i+1)*ny])
        dat = np.fromfile('torque_m%d.dat'%self.n)
        attrs = ['y','dbar','Lt','Ltn','Ld','Ldn','Lw','Lwn',
                'drFt','drFd','drFw',
                'Lamex','indLam','Lamdep',
                'dtLt','dtLd','dtLw','mdot']
        self.mmax = 30
        self.trq_drFdm = dat[:ny*(self.mmax+2)].reshape(self.mmax+2,ny)
        self.trq_Lamexm = dat[ny*(self.mmax+2):2*ny*(self.mmax+2)].reshape(self.mmax+2,ny)

        self.trq_dtLtm = dat[2*ny*(self.mmax+2):3*ny*(self.mmax+2)].reshape(self.mmax+2,ny)

        self.trq_drFtm = self.trq_Lamexm - self.trq_dtLtm
        self.trq_drFwm =self.trq_drFtm - self.trq_drFdm
        self.trq_dtLdm = np.zeros(self.trq_drFdm.shape)
        self.trq_dtLdm[0,:] = self.trq_dtLd
        self.trq_dtLwm = self.trq_dtLt - self.trq_dtLdm
        self.Lamdepm = self.trq_dtLdm + self.trq_drFdm

        #self.ymin = dat[-2*(ny+1):-(ny+1)]
        #self.mdot= dat[-(ny+1):]

        if not self.indirectterm:
            self.indLam = np.zeros(self.indLam.shape)
        return
    def plot(self,**kargs):
        figsize = kargs.pop('figsize',(20,10))
        fig,axes = plt.subplots(3,2,sharex='col',figsize=figsize)

        fontsize = kargs.pop('fontsize',14)
        logx = kargs.pop('logx',False)

        axes[0,0].plot(self.y,self.dtLt + self.drFt,label='$\\dot{L}_T + div(F_T)$')
        axes[0,0].plot(self.y,self.Lamex+self.indLam,label='$\\Lambda_{ex}$')
        axes[1,0].plot(self.y,self.dtLd + self.drFd,label='$\\dot{L}_d + div(F_d)$')
        axes[1,0].plot(self.y,self.Lamdep,label='$\\Lambda_{dep}$')
        axes[2,0].plot(self.y,self.dtLw + self.drFw,label='$\\dot{L}_w + div(F_w)$')
        axes[2,0].plot(self.y,self.Lamex+self.indLam-self.Lamdep,label='$\\Lambda_{ex}-\\Lambda_{dep}$')

        axes[0,1].plot(self.y,self.dtLt,label='$\\dot{L}_T$')
        axes[0,1].plot(self.y,self.drFt,label='$div(F_T)$')
        axes[0,1].plot(self.y,self.Lamex+self.indLam,label='$\\Lambda_{ex}$')
        axes[1,1].plot(self.y,self.dtLd,label='$\\dot{L}_d$')
        axes[1,1].plot(self.y,self.drFd,label='$div(F_d)$')
        axes[1,1].plot(self.y,self.Lamdep,label='$\\Lambda_{dep}$')
        axes[2,1].plot(self.y,self.dtLw,label='$\\dot{L}_w$')
        axes[2,1].plot(self.y,self.drFw,label='$div(F_w)$')
        axes[2,1].plot(self.y,self.Lamex+self.indLam,label='$\\Lambda_{ex}$')
        axes[2,1].plot(self.y,self.Lamdep,label='$\\Lambda_{dep}$')

        axes[-1,0].set_xlabel('$r$',fontsize=fontsize)
        axes[-1,1].set_xlabel('$r$',fontsize=fontsize)
        axes[0,0].set_ylabel('Total',fontsize=fontsize)
        axes[1,0].set_ylabel('Disk',fontsize=fontsize)
        axes[2,0].set_ylabel('Wave',fontsize=fontsize)

        for ax in axes.flatten():
            ax.minorticks_on()
            if logx:
                ax.set_xscale('log')
            ax.legend(loc='best')
        return
    def plot_lam(self,planet=None,plotdens=True,indirect=False,logx=None,fig=None,ax=None,xlims=None,ylims=None,**kargs):
        figsize = kargs.pop('figsize',(10,8))
        fontsize = kargs.pop('fontsize',15)
        if fig is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)

        if planet is None:
            x = self.y
            xstr = '$r$'
        else:
            x = (self.y - planet)/(planet*self.h)
            xstr = '$(r-a)/H_p$'

        if xlims is not None:
            ind = (x >= xlims[0])&(x<= xlims[1])
        else:
            ind = np.array([True]*len(self.y))

        ax.plot(x[ind],self.Lamdep[ind],'-g',label='$\\Lambda_{dep}$',**kargs)
        ax.plot(x[ind],self.drFw[ind],'--r',label='div($F_w$)',**kargs)
        if indirect:
            ax.plot(x[ind],self.Lamex[ind]+self.indLam[ind],'-b',label='$\\Lambda_{ex}$',**kargs)
            ax.plot(x[ind],self.dtLw[ind],'--m',label='$\\dot{L}_w$',**kargs)

        else:
            ax.plot(x[ind],self.Lamex[ind],'-b',label='$\\Lambda_{ex}$',**kargs)
            ax.plot(x[ind],self.dtLw[ind]-self.indLam[ind],'--m',label='$\\dot{L}_w$',**kargs)
        if logx and planet is not None:
            ax.set_xscale('log')
        if xlims is not None:
            ax.set_xlim(xlims)
        if ylims is not None:
            ax.set_ylim(ylims)
        ax.minorticks_on()

        ax.set_xlabel(xstr,fontsize=fontsize)
        ax.legend(loc='best')

        if plotdens:
            ax1 = ax.twinx()
            ax1.plot(x[ind],self.dbar[ind],'--k')
            ax1.minorticks_on()
            if logx and planet is not None:
                ax1.set_xscale('log')

        return

def make_summary_plots_f(dirlist,numlist):
    tot = len(dirlist)

    print 'Working on:'
    for i,(d,n) in enumerate(zip(dirlist,numlist)):
        fname = d  + d.split('/')[-2] + '_f_%d.png'%n
        print n,d,fname,
        sim = Sim(n,directory=d)
        sim.summary_plot(savefile=fname)
        print 'Finished %d of %d'%(i+1,tot)

    return

def make_summary_plots(dirlist):
    tot = len(dirlist)
    for i,d in enumerate(dirlist):
        dirname= d.split('/')[-2]
        times=glob.glob(d + 'gasdens*.dat')
        nt = len(times) - 1
        if nt >0:
            fname = d + dirname + '_%d.png'%nt
            sim = Sim(nt,directory=d)
            sim.summary_plot(savefile=fname)
            print 'Finished %d of %d'%(i+1,tot)
        else:
            print 'No outputs in %s'%dirname
    return


def mdot_comparison(sims, ylims=(-5,5), savefig=None):
    fig,axes=plt.subplots(1,3,figsize=(20,6))

    alphas = np.array([ s.params.alpha for s in sims])
    q = np.array([ s.mp for s in sims])
    K = np.array([ s.K for s in sims])
    h = np.array([ s.params.aspectratio for s in sims])
    times = np.array([ s.t/s.tviscp for s in sims])
    gdepthFL = .14* (alphas/1e-2)**1.41 * (q/1e-3)**(-2.16) * (h/.05)**(6.61)
    gdepthFH  = 4.7e-3 * (q/5e-3)**(-1) * (alphas/1e-2)**(1.26) *(h/.05)**(6.12)
    gdepthF = (q>=5e-3).astype(int) * gdepthFH + (q < 5e-3).astype(int)*gdepthFL

    gdepthD = 1./(1 + (.451 + h*.561)*K/(3*np.pi))


    for s in sims:
        axes[0].plot(s.dens.y,s.dens.avg/s.sig0)
        axes[1].plot(s.dens.y,s.trq_mdot/s.params.mdot)
        t,d = s.dens_evolution()
        axes[2].plot(t/(2*np.pi),d)

    cvals = [ l.get_color() for l in axes[0].lines]
    for g,gd,c in zip(gdepthF,gdepthD,cvals):
        axes[0].axhline(g,color=c,linestyle='--')

    axes[1].axhline(1,color='k')

    axes[0].set_xlabel('$r$',fontsize=20)
    axes[1].set_xlabel('$r$',fontsize=20)
    axes[2].set_xlabel('t (yrs)',fontsize=15)

    axes[0].set_ylabel('$\\Sigma/\\Sigma_0$',fontsize=20)
    axes[1].set_ylabel('$\\dot{M}/\\dot{M}_0$',fontsize=20)
    axes[2].set_ylabel('$\\left( \\Sigma / \\Sigma_0 \\right)_{min}$',fontsize=20)

    axes[1].set_ylim(ylims)
    axes[0].set_yscale('log')
    axes[2].set_xscale('log')
    axes[2].set_yscale('log')

    if savefig:
        fig.savefig(savefig)


def make_comparison_plots(sim1,sim2,ylims=None,savefig=None):

    sims = [sim1,sim2]

    fig,axes=plt.subplots(3,2,sharex='col',figsize=(20,10))

    for i,s in enumerate(sims):

        axes[0,i].plot(y,dbar,'.k',label='Fargo')
        axes[0,i].plot(y2,sig2,'-b',linewidth=2,label='Excluding WKZ')
        axes[0,i].plot(y,sig,'--r',linewidth=2,label='Including WKZ')

        axes[1,i].plot(y2,lex2,'-k',label='$\\Lambda_{ex}$')
        axes[1,i].plot(y2,ldep2,'-b',label='$\\Lambda_{dep}$')
        axes[1,i].plot(y2,drFw2,'-r',label='$\\frac{1}{r}\\partial_r(r F_w)$')

        axes[1,i].plot(yL,lexL,'--k')
        axes[1,i].plot(yL,ldepL,'--b')
        axes[1,i].plot(yL,drFwL,'--r')
        axes[1,i].plot(yR,lexR,'--k')
        axes[1,i].plot(yR,ldepR,'--b')
        axes[1,i].plot(yR,drFwR,'--r')

        axes[2,i].plot(y2,mdot2,'-k')
        axes[2,i].plot(yL,mdotL,'-k')
        axes[2,i].plot(yR,mdotR,'-k')
        axes[2,i].axhline(1,color='k',linestyle='--')

        axes[2,i].set_xlabel('$r$',fontsize=20)
        axes[2,i].set_ylim(0,3)
        axes[0,i].axvline(wkzL,linestyle='--',color='k')
        axes[1,i].axvline(wkzL,linestyle='--',color='k')
        axes[2,i].axvline(wkzL,linestyle='--',color='k')
        axes[0,i].axvline(wkzR,linestyle='--',color='k')
        axes[1,i].axvline(wkzR,linestyle='--',color='k')
        axes[2,i].axvline(wkzR,linestyle='--',color='k')

        axes[0,i].set_xlim(0,s.dens.ymax)






    axes[0,0].set_ylabel('$\\Sigma/\\Sigma_0$',fontsize=20)
    axes[1,0].set_ylabel('Torque density',fontsize=15)
    axes[2,0].set_ylabel('$\\dot{M}/\\dot{M}_o$',fontsize=20)

    axes[0,0].set_title('Damping $v_R$ and $\\Sigma$',fontsize=20)
    axes[0,1].set_title('Damping $v_R$',fontsize=20)

    axes[0,1].legend(loc='lower right')
    axes[1,1].legend(loc='upper right')
    if ylims is not None:
        axes[1,0].set_ylim(ylims)
        axes[1,1].set_ylim(ylims)


    for ax in axes.flatten():
        ax.minorticks_on()
    plt.subplots_adjust(wspace=0)
    plt.setp(axes[2,0].get_xticklabels()[-1],visible=False)

    if savefig is not None:
        fig.savefig(savefig)

def make_ss_plots(sim1,sim2,ylims=None,savefig=None):

    sims = [sim1,sim2]

    fig,axes=plt.subplots(3,2,sharex='col',sharey='row',figsize=(20,10))

    for i,s in enumerate(sims):
        ind = (s.trq_y > s.dens.wkz)&(s.trq_y < s.dens.wkzr)
        indL = (s.trq_y <= s.dens.wkz + s.dlr[0]*s.dens.wkz)
        indR = (s.trq_y >= s.dens.wkzr - s.dlr[0]*s.dens.wkzr)

        intLam = (s.trq_y*s.dens.dy*2*pi*s.trq_Lamdep).cumsum()
        intLam -= intLam[0]
        intLam2 = (s.trq_y*s.dens.dy*2*pi*s.trq_Lamdep)[ind].cumsum()
        intLam2 -= intLam2[0]

        sig = 1 + intLam/(s.dens.mdot *np.sqrt(s.trq_y))
        sig2 = 1 + intLam2/(s.dens.mdot *np.sqrt(s.trq_y[ind]))

        mdot = s.trq_mdot/s.params.mdot
        mdot2 = s.trq_mdot[ind]/s.params.mdot
        mdotL = s.trq_mdot[indL]/s.params.mdot
        mdotR = s.trq_mdot[indR]/s.params.mdot

        ldep = s.trq_Lamdep.copy()
        ldep2 = s.trq_Lamdep[ind]
        ldepR = s.trq_Lamdep[indR]
        ldepL = s.trq_Lamdep[indL]
        lex = s.trq_Lamex.copy()
        lex2 = s.trq_Lamex[ind]
        lexL = s.trq_Lamex[indL]
        lexR = s.trq_Lamex[indR]

        drFw = s.trq_drFw.copy()
        drFw2 = s.trq_drFw[ind]
        drFwL = s.trq_drFw[indL]
        drFwR = s.trq_drFw[indR]

        y = s.trq_y.copy()
        y2 = s.trq_y[ind]
        yL = s.trq_y[indL]
        yR = s.trq_y[indR]

        dbar = s.dens.avg/s.sig0

        wkzL = s.dens.wkz
        wkzR = s.dens.wkzr


        axes[0,i].plot(y,dbar,'.k',label='Fargo')
        axes[0,i].plot(y2,sig2,'-b',linewidth=2,label='Excluding WKZ')
        axes[0,i].plot(y,sig,'--r',linewidth=2,label='Including WKZ')

        axes[1,i].plot(y2,lex2,'-k',label='$\\Lambda_{ex}$')
        axes[1,i].plot(y2,ldep2,'-b',label='$\\Lambda_{dep}$')
        axes[1,i].plot(y2,drFw2,'-r',label='$\\frac{1}{r}\\partial_r(r F_w)$')

        axes[1,i].plot(yL,lexL,'--k')
        axes[1,i].plot(yL,ldepL,'--b')
        axes[1,i].plot(yL,drFwL,'--r')
        axes[1,i].plot(yR,lexR,'--k')
        axes[1,i].plot(yR,ldepR,'--b')
        axes[1,i].plot(yR,drFwR,'--r')

        axes[2,i].plot(y2,mdot2,'-k')
        axes[2,i].plot(yL,mdotL,'-k')
        axes[2,i].plot(yR,mdotR,'-k')
        axes[2,i].axhline(1,color='k',linestyle='--')

        axes[2,i].set_xlabel('$r$',fontsize=20)
        axes[2,i].set_ylim(0,3)
        axes[0,i].axvline(wkzL,linestyle='--',color='k')
        axes[1,i].axvline(wkzL,linestyle='--',color='k')
        axes[2,i].axvline(wkzL,linestyle='--',color='k')
        axes[0,i].axvline(wkzR,linestyle='--',color='k')
        axes[1,i].axvline(wkzR,linestyle='--',color='k')
        axes[2,i].axvline(wkzR,linestyle='--',color='k')

        axes[0,i].set_xlim(0,s.dens.ymax)






    axes[0,0].set_ylabel('$\\Sigma/\\Sigma_0$',fontsize=20)
    axes[1,0].set_ylabel('Torque density',fontsize=15)
    axes[2,0].set_ylabel('$\\dot{M}/\\dot{M}_o$',fontsize=20)

    axes[0,0].set_title('Damping $v_R$ and $\\Sigma$',fontsize=20)
    axes[0,1].set_title('Damping $v_R$',fontsize=20)

    axes[0,1].legend(loc='lower right')
    axes[1,1].legend(loc='upper right')
    if ylims is not None:
        axes[1,0].set_ylim(ylims)
        axes[1,1].set_ylim(ylims)


    for ax in axes.flatten():
        ax.minorticks_on()
    plt.subplots_adjust(wspace=0)
    plt.setp(axes[2,0].get_xticklabels()[-1],visible=False)

    if savefig is not None:
        fig.savefig(savefig)

def torque_balance(sim):
    fig,axes=plt.subplots(5,1,sharex=True)

    axes[0].plot(sim.trq_y,sim.trq_Lamex)
    axes[0].plot(sim.trq_y,sim.trq_drFt + sim.trq_dtLt)
    axes[1].plot(sim.trq_y,sim.trq_Lamex-sim.trq_Lamdep)
    axes[1].plot(sim.trq_y,sim.trq_drFw + sim.trq_dtLw)
    axes[2].plot(sim.trq_y,sim.trq_Lamdep)
    axes[2].plot(sim.trq_y,sim.trq_drFd + sim.trq_dtLd)

    for i in range(1,7):
        axes[3].plot(sim.trq_y,getattr(sim,'trq_Lamdep%d'%i),label=str(i))

    axes[3].legend(loc='best')

    axes[4].plot(sim.trq_y,sim.trq_Lamtot)
    axes[4].plot(sim.trq_y,sim.trq_Lamdep)


    axes[-1].set_xlabel('$r$',fontsize=20)

    axes[0].set_ylabel('Total',fontsize=15)
    axes[1].set_ylabel('Waves',fontsize=15)
    axes[2].set_ylabel('Disk',fontsize=15)

def compare_torque_ab(sim,savefig=None):
    y = sim.trq_y[ind].copy()
    lamdepA = sim.trq_Lamdep[ind].copy()
    lamdepB = sim.trq_LamdepB[ind].copy()

    drFwA = sim.trq_drFw[ind].copy()
    drFwB = sim.trq_drFwB[ind].copy()

    lamex = sim.trq_Lamex[ind].copy()

    dR = (sim.trq_y*sim.dens.dy)[ind]
    intLA = (dR * lamdepA).cumsum()
    intLB = (dR * lamdepB).cumsum()

    intLex = (dR * lamex).cumsum()
    intLex -= intLex[0];
    intFwA = (dR * drFwA).cumsum()
    intFwA -= intFwA[0]

    intRA = (dR * lamex - dR * drFwA).cumsum()
    intRB = (dR * lamex - dR * drFwB).cumsum()

    intFwB = (dR * drFwB).cumsum()
    intFwB -= intFwB[0]


    intLA -= intLA[0]
    intLB -= intLB[0]
    intRA -= intRA[0]
    intRB -= intRB[0]


    fig,axes = plt.subplots(2,2,sharex=True,figsize=(25,15))
    plt.subplots_adjust(hspace=0)

    axes[0,0].plot(y,lamex,'-k',label='$\\Lambda_{ex}$')
    axes[0,0].plot(y,drFwA,'-r',label='$\\nabla_r (F_w)$')
    axes[0,0].plot(y,lamdepA,'-b',label='$\\Lambda_{dep}$')

    axes[0,1].plot(y,lamex,'-k',label='$\\Lambda_{ex}$')
    axes[0,1].plot(y,drFwB,'-r',label='$\\nabla_r (F_w)$')
    axes[0,1].plot(y,lamdepB,'-b',label='$\\Lambda_{dep}$')

    axes[1,0].plot(y,intRA,'-m',label='$\\int \, dr \\, r \\Lambda_{ex} - r \\nabla_r(F_w)$')
    axes[1,0].plot(y,intLA,'-b',label='$\\int \, dr \\, r \\Lambda_{dep}$')
    axes[1,0].plot(y,intLex,'-k',label='$\\int \, dr \\, r \\Lambda_{ex}$')
    axes[1,0].plot(y,intFwA,'-r',label='$\\int \, dr \\, r \\nabla_r(F_w)$')

    axes[1,1].plot(y,intRB,'-m',label='$\\int \, dr \\, r \\Lambda_{ex} - r \\nabla_r(F_w)$')
    axes[1,1].plot(y,intLB,'-b',label='$\\int \, dr \\, r \\Lambda_{dep}$')
    axes[1,1].plot(y,intLex,'-k',label='$\\int \, dr \\, r \\Lambda_{ex}$')
    axes[1,1].plot(y,intFwB,'-r',label='$\\int \, dr \\, r \\nabla_r(F_w)$')


    axes[1,0].set_xlabel('$r$',fontsize=20)
    axes[1,1].set_xlabel('$r$',fontsize=20)

    axes[0,0].set_title('Option A, $ \\langle \\Sigma v_r \\rangle$',fontsize=20)
    axes[0,1].set_title('Option B, $ \\langle \\Sigma \\rangle \\langle v_r \\rangle$',fontsize=20)

    axes[0,0].legend(loc='lower right')
    axes[1,0].legend(loc='lower right')

    plt.setp(axes[1,0].get_yticklabels()[-1],visible=False)
    plt.setp(axes[1,1].get_yticklabels()[-1],visible=False)
    plt.setp(axes[0,0].get_yticklabels()[0],visible=False)
    plt.setp(axes[0,1].get_yticklabels()[0],visible=False)
    fig.canvas.draw()

    if savefig is not None:
        fig.savefig(savefig)
    return



