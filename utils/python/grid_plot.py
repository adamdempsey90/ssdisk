import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy.interpolate import interp2d
from matplotlib.mlab import griddata

class grid:
    """Clase para manejar outputs de FARGO. """
    def __init__(self, input_file, nrad, nsec, rmin, rmax,zmin,zmax):

        self.input = input_file
        self.nrad = nrad
        self.nsec = nsec
        self.rmin = rmin
        self.rmax = rmax
        self.zmin = zmin
        self.zmax = zmax
        try:
            file = open(self.input, 'r')
        except:
            print self.input, 'No encontrado'
        self.data = np.fromfile(input_file, "d", self.nrad*self.nsec)
        self.data = np.reshape(self.data, [self.nrad,self.nsec])
        self.data = np.concatenate((np.array(self.data[:,int(nsec/2):nsec]), \
                                     np.array(self.data[:,0:int(nsec/2)])),1)
        self.data = np.log(self.data)

        self.maxval = np.max(self.data)
        self.minval = np.min(self.data)

    def plot(self):
        fig = plt.figure()
        ax1 = fig.add_subplot(111, axisbg='k')
        ax1.set_xlabel("Theta")
        ax1.set_ylabel("R")
        ax1.imshow(self.data, origin="lower", \
                       extent=[0,2*np.pi,self.rmin,self.rmax])
        plt.show()

    def polar(self, n, file='None'):
        r = np.arange(self.rmin,self.rmax,(self.rmax-self.rmin)/self.nrad)
        t = np.arange(0.,2.*np.pi,2.*np.pi/self.nsec)
        x = np.ndarray([self.nrad*self.nsec], dtype = float)
        y = np.ndarray([self.nrad*self.nsec], dtype = float)
        z = np.ndarray([self.nrad*self.nsec], dtype = float)
        k = 0
        for i in range(self.nrad):
            for j in range(self.nsec):
                x[k] = r[i]*np.cos(t[j])
                y[k] = r[i]*np.sin(t[j])
                z[k] = self.data[i,j]
                k +=1
        xx = np.arange(-self.rmax, self.rmax, (self.rmax-self.rmin)/n)
        yy = np.arange(-self.rmax, self.rmax, (self.rmax-self.rmin)/n)
        zz = griddata(x,y,z,xx,yy)
        fig_pol = plt.figure()
        ax1 = fig_pol.add_subplot(111, axisbg='k')
        ax1.set_xlabel("X")
        ax1.set_ylabel("Y")
        if(self.zmax!='None' and self.zmin!='None'):
            ax1.imshow(zz, cmap=cm.hot, origin="lower", \
                           extent=[-self.rmax,self.rmax, \
                                        -self.rmax,self.rmax])
        else:
            ax1.imshow(zz, cmap=cm.hot, origin="lower", \
                           extent=[-self.rmax,self.rmax, \
                                        -self.rmax,self.rmax])
        if(file!="None"):
            plt.savefig(file+".png",dpi=70, format="png" )
            print file+".png done"
        else:
            plt.show()

    def dual_plot(self, n):
        fig = plt.figure()
        ax1 = fig.add_subplot(121, axisbg='k')
        ax1.set_xlabel("Theta")
        ax1.set_ylabel("R")
        ax1.imshow(self.data, cmap=cm.hot,origin="lower", \
                       extent=[0,2*np.pi,self.rmin,self.rmax])
        r = np.arange(self.rmin,self.rmax,(self.rmax-self.rmin)/self.nrad)
        t = np.arange(0.,2.*np.pi,2.*np.pi/self.nsec)
        x = np.ndarray([self.nrad*self.nsec], dtype = float)
        y = np.ndarray([self.nrad*self.nsec], dtype = float)
        z = np.ndarray([self.nrad*self.nsec], dtype = float)
        k = 0
        for i in range(self.nrad):
            for j in range(self.nsec):
                x[k] = r[i]*np.cos(t[j])
                y[k] = r[i]*np.sin(t[j])
                z[k] = self.data[i,j]
                k +=1
        xx = np.arange(-self.rmax, self.rmax, (self.rmax-self.rmin)/n)
        yy = np.arange(-self.rmax, self.rmax, (self.rmax-self.rmin)/n)
        zz = griddata(x,y,z,xx,yy)
        ax2 = fig.add_subplot(122, axisbg='k')
        ax2.set_xlabel("X")
        ax2.set_ylabel("Y")
        ax2.imshow(zz, cmap=cm.hot,origin="lower" \
                       , extent=[-self.rmax,self.rmax,-self.rmax,self.rmax])
        plt.show()

if (__name__ == "__main__"):
    grilla_ref = grid("gasdens228.dat",306,384, 0.1,3.5)
    zmin = grilla_ref.minval
    zmax = grilla_ref.maxval
    for i in range(10):
        name = 'gasdens'+str(i)+'.dat'
        if(i<1000):
            name_out = 'gasdens'+str(i)
        if(i<100):
            name_out = 'gasdens0'+str(i)
        if(i<10):
            name_out = 'gasdens00'+str(i)
        grilla = grid(name,306,384, 0.1,3.5,zmin,zmax)
        grilla.polar(200,name_out)
