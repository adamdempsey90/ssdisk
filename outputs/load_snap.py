import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
class Field():

    def __init__(self,t,nsteps=8,dt=np.pi*1e-5,nx=384,ny=128+6):

        self.ymin = np.loadtxt('domain_y.dat')
        self.xmin = np.loadtxt('domain_x.dat')
        self.xmed = .5*(self.xmin[1:]+self.xmin[:-1])
        self.ymed = .5*(self.ymin[1:]+self.ymin[:-1])
        self.nx = nx
        self.ny = ny
        self.dt = dt
        self.t = t
        self.dens = [None]*nsteps
        self.vx = [None]*nsteps
        self.vy = [None]*nsteps
        self.pot = [None]*nsteps

        s = nx*ny
        for i in range(nsteps):
            dat = np.fromfile('substep_%d_%d.dat'%(i,t))
            self.dens[i] = dat[:s].reshape(ny,nx)
            self.vy[i] = dat[s:2*s].reshape(ny,nx)
            self.vx[i] = dat[2*s:3*s].reshape(ny,nx)
            self.pot[i] = dat[3*s:4*s].reshape(ny,nx)
        return

class Snap():

    def __init__(self,directory='temp_files/',nx=384,ny=128+6):

        dat = np.loadtxt('param_file.txt')
        self.alpha = dat[3]
        self.mp = dat[4]
        self.q = dat[5]
        self.omf = dat[6]
        self.flaringindex = dat[7]
        self.soft = dat[9]
        self.nx = nx
        self.ny = ny
        s = nx*ny
        dat = np.fromfile(directory+'output.dat')
        self.ymin = dat[:(ny+1)]
        dat=dat[(ny+1):]
        self.xmin = dat[:(nx+1)]
        dat=dat[(nx+1):]
        self.dens = dat[:s].reshape(ny,nx)
        self.vy = dat[s:2*s].reshape(ny,nx)
        self.vx= dat[2*s:3*s].reshape(ny,nx)
        self.pot= dat[3*s:4*s].reshape(ny,nx)

        self.ymed = .5*(self.ymin[1:]+self.ymin[:-1])
        self.xmed = .5*(self.xmin[1:] + self.xmin[:-1])

        dat = np.fromfile(directory+'output_init.dat')
        self.dens0 = dat[:s].reshape(ny,nx)
        self.vy0 = dat[s:2*s].reshape(ny,nx)
        self.vx0= dat[2*s:3*s].reshape(ny,nx)
        return
class Fargo():

    def __init__(self,i,directory='',nx=384,ny=128):
        dat = np.loadtxt('param_file.txt')
        self.alpha = dat[3]
        self.mp = dat[4]
        self.q = dat[5]
        self.omf = dat[6]
        self.flaringindex = dat[7]
        self.soft = dat[9]
        self.dt = dat[11]


        self.nx = nx
        self.ny = ny
        self.dens = np.fromfile('gasdens%d.dat'%i).reshape(ny,nx)
        self.vy = np.fromfile('gasvy%d.dat'%i).reshape(ny,nx)
        self.vx = np.fromfile('gasvx%d.dat'%i).reshape(ny,nx)

        self.ymin = np.loadtxt('domain_y.dat')[3:-3]
        self.xmin = np.loadtxt('domain_x.dat')
        self.xmed = .5*(self.xmin[1:]+self.xmin[:-1])
        self.ymed = .5*(self.ymin[1:]+self.ymin[:-1])
        self.L = self.dens *( .5*( np.roll(self.vx,1,axis=1)+self.vx) + self.omf*self.ymed[:,np.newaxis])*self.ymed[:,np.newaxis]

        return
class Torque():

    def __init__(self,directory='temp_files/'):

        dat = np.loadtxt('param_file.txt')
        nx = dat[0]
        ny = dat[1]
        self.nx = nx
        self.ny = ny
        self.nz = dat[2]
        self.alpha = dat[3]
        self.mp = dat[4]
        self.a = dat[5]
        self.h = dat[6]
      #  self.omf = dat[6]
        self.flaringindex = dat[7]
        self.soft = dat[9]
        self.indirectterm = bool(dat[10])
      #  self.dt = dat[11]
        dat = np.fromfile(directory+'torque.dat')
        attrs = ['y','dbar','Lt','Ltn','Ld','Ldn','Lw','Lwn',
                'drFt','drFd','drFw',
                'Lamex','indLam','Lamdep',
                'dtLt','dtLd','dtLw']

        for i,a in enumerate(attrs):
            setattr(self,a,dat[i*ny:(i+1)*ny])

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

def compare(dat,fld,i=-1,relerror=False):
    fig,axes=plt.subplots(1,4,sharey=True,figsize=(20,10))
    if type(fld.dens) != list:
        dens = fld.dens
        vx = fld.vx
        vy = fld.vy
        pot = fld.pot
    else:
        dens = fld.dens[i]
        vx = fld.vx[i]
        vy = fld.vy[i]
        pot = fld.pot[i]

    try:
        errd = (dat.dens[3:-3,:] - dens)
        errvy = (dat.vy[3:-3,:] - vy)
        errvx = (dat.vx[3:-3,:] - vx)
        errp = (dat.pot[3:-3,:] - pot)
    except ValueError:
        dens = dens[3:-3,:]
        vx = vx[3:-3,:]
        vy = vy[3:-3,:]
        pot = pot[3:-3,:]
        errd = (dat.dens[3:-3,:] - dens)
        errvy = (dat.vy[3:-3,:] - vy)
        errvx = (dat.vx[3:-3,:] - vx)
        errp = (dat.pot[3:-3,:] - pot)

    if relerror:
        errd /= dens
        errvy /= vy
        errvx /= vx
        errp /= pot


    ims = [None]*len(axes)
    ims[0]=axes[0].imshow(np.log10(abs(errd)),origin='lower' ,aspect='auto')
    ims[1]=axes[1].imshow(np.log10(abs(errvy)),origin='lower',aspect='auto')
    ims[2]=axes[2].imshow(np.log10(abs(errvx)),origin='lower',aspect='auto')
    ims[3]=axes[3].imshow(np.log10(abs(errp)),origin='lower',aspect='auto')
    axes[0].set_title('$\\Delta\\Sigma$',fontsize=15,rotation=0)
    axes[1].set_title('$\\Delta v_y$',fontsize=15,rotation=0)
    axes[2].set_title('$\Delta v_x$',fontsize=15,rotation=0)
    axes[3].set_title('$\Delta \Phi$',fontsize=15,rotation=0)
    dividers = [make_axes_locatable(ax) for ax in axes]
    caxes = [d.append_axes("right", size="5%", pad=0.05) for d in dividers]
    cbars = [plt.colorbar(im, cax=cax) for im,cax in zip(ims,caxes)]
    for cb,ax in zip(cbars,axes):
        ax.minorticks_on()
        #cb.ax.xaxis.set_ticks_position('top')
        #ax.yaxis.set_label_position('right')
        #ax.yaxis.labelpad = 20


