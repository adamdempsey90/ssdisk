import matplotlib.pyplot as plt
import numpy as np
from matplotlib.mlab import griddata


#object.__add__(self, other)
#object.__sub__(self, other)
#object.__mul__(self, other)
#object.__floordiv__(self, other)
#object.__mod__(self, other)
#object.__divmod__(self, other)
#object.__pow__(self, other[, modulo])
#object.__lshift__(self, other)
#object.__rshift__(self, other)
#object.__and__(self, other)
#object.__xor__(self, other)
#object.__or__(self, other)


class animate(object):
    def __init__(self, init = 0, end = 0, i = None, j = None, k = None, \
                   imout = None, dpi=72, colorbox = False, \
                   projection = 'standard', n = 500, \
                   zmin = None, zmax = None, palette = None):

        self.init = init
        self.end  = end

        for l in range(init, end+1):
            a = SnapShot(l)            
            filename = "snap" + "{:06d}".format(l) + '.png'
            a.rho.plot2d(i = i, j = j, k = k, \
                   imout = filename, dpi = dpi, colorbox = colorbox, \
                   projection = projection, n = 500, \
                   zmin = zmin, zmax = zmax, palette = palette)

class Field(object):
    def __add__(self, other):
        data = self.data+other.data
        return Field(self.dim, self.dom, data)

    def __sub__(self, other):
        data = self.data-other.data
        return Field(self.dim, self.dom, data)

    def __mul__(self, other):
        data = self.data*other.data
        return Field(self.dim, self.dom, data)

    def __div__(self, other):
        data = self.data/other.data
        return Field(self.dim, self.dom, data)

    def __init__ (self, dim, dom, data = None):
        if(data == None):
            self.data = np.ndarray([dim[2],dim[1],dim[0]], dtype = "d")
        else:
            data = np.reshape(data, dim)
            self.data = data
        self.dom = dom
        self.dim = dim

    def project_data (self, data, dim1, dim2, xmin, xmax, ymin, ymax, n, \
                          projection = 'meridian'):
        if(projection == 'meridian'):
            r = np.arange(xmin, xmax, (xmax-xmin)/dim1)
            t = np.arange(ymin, ymax, (ymax-ymin)/dim2)
            x = np.ndarray([dim1*dim2], dtype = 'd')
            y = np.ndarray([dim1*dim2], dtype = 'd')
            z = np.ndarray([dim1*dim2], dtype = 'd')
            k = 0
            for i in range(dim1):
                for j in range(dim2):
                    x[k] = r[i]*np.sin(t[j])
                    y[k] = r[i]*np.cos(t[j])
                    z[k] = data[j,i]
                    k += 1
            xx = np.arange(xmin, xmax, (xmax-xmin)/n)
            yy = np.arange(xmax*np.cos(ymax), xmax*np.cos(ymin),\
                               (xmax*(np.cos(ymin)-np.cos(ymax))/n))
            data = griddata(x,y,z,xx,yy)

        if(projection == 'polar'):
            newdim = (dim2+int(ymin*dim2/(ymax-ymin))) #Correction for black centre
            r = np.arange(0.0,ymax,ymax/newdim)
            t = np.arange(xmin,xmax,(xmax-xmin)/dim1)
            x = np.ndarray([newdim*dim1], dtype = float)
            y = np.ndarray([newdim*dim1], dtype = float)
            z = np.ndarray([newdim*dim1], dtype = float)
            k = 0
            for i in range(newdim):
                for j in range(dim1):
                    x[k] = r[i]*np.cos(t[j])
                    y[k] = r[i]*np.sin(t[j])
                    if(r[i]>ymin):
                        z[k] = data[i-(newdim-dim2),j]
                    else:
                        z[k] = 0.0
                    k += 1
            xx = np.arange(-ymax,ymax,2.0*ymax/n)
            yy = np.arange(-ymax,ymax,2.0*ymax/n)
            data = griddata(x,y,z,xx,yy)
        return data

    def plot2d(self, i = None, j = None, k = None, \
                   imout = None, dpi=72, colorbox = False, \
                   projection = 'standard', n = 500, \
                   zmin = None, zmax = None, palette = None):
        
        fig  = plt.figure()
        axis = fig.add_subplot(111)

        if(i !=None ):
            dim1 = self.dim[1]
            dim2 = self.dim[0]
            data = np.ndarray([dim1,dim2], dtype = 'd')
            data = self.data[:,:,i]
            domain = [self.dom[2], self.dom[3], \
                      self.dom[4], self.dom[5]]
            self.x_axis = 'Y'
            self.y_axis = 'Z'
        elif(j !=None ):
            dim1 = self.dim[0]
            dim2 = self.dim[2]
            data = np.ndarray([dim1,dim2], dtype = 'd')
            data = self.data[:,j,:]
            domain = [self.dom[0], self.dom[1], \
                      self.dom[4], self.dom[5]]
            self.x_axis = 'X'
            self.y_axis = 'Z'
        elif(k !=None ):
            dim1 = self.dim[2]
            dim2 = self.dim[1]
            data = np.ndarray([dim1,dim2], dtype = "d")
            data = self.data[k,:,:]
            domain = [self.dom[0], self.dom[1], \
                      self.dom[2], self.dom[3]]
            self.x_axis = 'X'
            self.y_axis = 'Y'
        else: #default plane
            if(self.dim[2] == 1):
                dim1 = self.dim[1]
                dim2 = self.dim[0]
                data = np.ndarray([dim1,dim2], dtype = "d")
                data = self.data[:,:,0]
                self.x_axis = 'Y'
                self.y_axis = 'Z'
                domain = [self.dom[2], self.dom[3], \
                          self.dom[4], self.dom[5]]
            elif(self.dim[1] == 1):
                dim1 = self.dim[2]
                dim2 = self.dim[0]
                data = np.ndarray([dim1,dim2], dtype = "d")
                data = self.data[:,0,:]
                self.x_axis = 'X'
                self.y_axis = 'Z'
                domain = [self.dom[0], self.dom[1], \
                          self.dom[4], self.dom[5]]
            elif(self.dim[0] == 1):
                dim1 = self.dim[2]
                dim2 = self.dim[1]
                data = np.ndarray([dim1,dim2], dtype = "d")
                data = self.data[0,:,:]
                self.x_axis = 'X'
                self.y_axis = 'Y'
                domain = [self.dom[0], self.dom[1], \
                          self.dom[2], self.dom[3]]
            else:
                dim1 = self.dim[2]
                dim2 = self.dim[1]
                data = np.ndarray([dim1,dim2], dtype = "d")
                data = self.data[int(self.dim[0]/2),:,:]
                self.x_axis = 'X'
                self.y_axis = 'Y'
                domain = [self.dom[0], self.dom[1], \
                          self.dom[2], self.dom[3]]

        axis.set_xlabel(self.x_axis)
        axis.set_ylabel(self.y_axis)        

        if(projection == 'meridian'):
            newdata = self.project_data(data, dim1, dim2, \
                                            self.dom[2], self.dom[3], \
                                            self.dom[4], self.dom[5], \
                                            n, projection)
            print domain[2], domain[3]
            temp = domain[2]
            domain[2] = domain[3]
            domain[3] = temp
            domain[2] = np.cos(domain[2])*domain[1]
            domain[3] = np.cos(domain[3])*domain[1]
            im  = plt.imshow(newdata, extent = domain, \
                                 origin = 'lower', vmin = zmin, vmax = zmax)
        elif (projection == 'polar'):
            newdata = self.project_data(data, dim1, dim2, \
                                            self.dom[0], self.dom[1], \
                                            self.dom[2], self.dom[3], \
                                            n, projection)
            im  = plt.imshow(newdata, extent = [-self.dom[3], self.dom[3],-self.dom[3],self.dom[3]],\
                                 origin = 'lower', vmin = zmin, vmax = zmax)
        else:
            im  = plt.imshow(data, extent = domain, origin = 'lower', vmin = zmin, vmax = zmax)
        if(palette != None):
            im.set_cmap(palette)
        if(colorbox == True):
            cbar = fig.colorbar(im)
        if(imout == None):            
            plt.show()
        else:
            plt.savefig(imout, dpi =dpi, format = 'png')
            print imout
        plt.close(fig)
        
class SnapShot(object):
    """
    This class is the core of fargo3d class. This class charge all material
    necessary for the analysis of a simulation. It uses dims.dat information.
    """

    def CompVort(self, k = 0, coord = 'spherical'):

        if(coord == 'cylindrical'):
            DrvpDr = np.ndarray([self.dim[2],self.dim[1]-1],  dtype = "d")
            DvrDp  = np.ndarray([self.dim[2],self.dim[1]-1],  dtype = "d")
            vorticity = np.ndarray([self.dim[2]*self.dim[1]],  dtype = "d")
            rho = self.rho.data[:,:,:]
            try:
                vp = self.vx.data[:,:,:]
            except:
                print "Vortencity error! Vx needed"
            try:
                vr = self.vy.data[:,:,:]
            except:
                print "Vortencity error! Vy needed"
            
            for j in range(self.dim[1]-1):
                for i in range(self.dim[2]):
                    if(i == self.dim[2]-1):
                        ip = i-(self.dim[2]-1)
                    else:            
                        ip = i+1
                    l = i+j*self.dim[2]
                    DvrDp[i,j] = (vr[k,j,ip]-vr[k,j,i]) \
                        /(self.X[i]-self.X[i-1]) # k,j,i
                    DrvpDr[i,j] = (self.Y[j+1]*vp[k,j+1,i]-self.Y[j]*vp[k,j,i]) \
                        /(self.Y[j+1]-self.Y[j])
                    vorticity[l] = 1./self.Y[j]*(DvrDp[i,j]-DrvpDr[i,j])
                    vorticity[l] /= rho[k,j,i]

        if(coord == 'spherical'):
            DrvpDr = np.ndarray([self.dim[2],self.dim[1]-1],  dtype = "d")
            DvrDp  = np.ndarray([self.dim[2],self.dim[1]-1],  dtype = "d")
            vorticity = np.ndarray([self.dim[2]*self.dim[1]],  dtype = "d")
            rho = self.rho.data[:,:,:]
            try:
                vp = self.vx.data[:,:,:]
            except:
                print "Vortencity error! Vx needed"
            try:
                vr = self.vy.data[:,:,:]
            except:
                print "Vortencity error! Vy needed"
            
            for j in range(self.dim[1]-1):
                for i in range(self.dim[2]):
                    if(i == self.dim[2]-1):
                        ip = i-(self.dim[2]-1)
                    else:
                        ip = i+1
                    l = i+j*self.dim[2]
                    DvrDp[i,j] = (vr[k,j,ip]-vr[k,j,i]) \
                        /(self.X[i]-self.X[i-1]) # k,j,i
                    DrvpDr[i,j] = (self.Y[j+1]*vp[k,j+1,i]-self.Y[j]*vp[k,j,i]) \
                        /(self.Y[j+1]-self.Y[j])
                    vorticity[l] = 1./self.Y[j+1] * \
                        ((1./np.sin(self.Z[k])*DvrDp[i,j]-DrvpDr[i,j]))
                    vorticity[l] /= rho[k,j,i]
                    #I do not interpolate density to edge
                    #because derivatives are evaluated in different edges.
                    #If I want to a best aproximation, 
                    #I think that I need a double
                    #interpolation (each one for each border)

        if(coord == 'cartesian'):
            DvxDy  = np.ndarray([self.dim[2],self.dim[1]-1],  dtype = "d")
            DvyDx  = np.ndarray([self.dim[2],self.dim[1]-1],  dtype = "d")
            vorticity = np.ndarray([self.dim[2]*self.dim[1]],  dtype = "d")
            rho = self.rho.data[:,:,:]
            try:
                vx = self.vx.data[:,:,:]
            except:
                print "Vortencity error! Vx needed"
            try:
                vy = self.vy.data[:,:,:]
            except:
                print "Vortencity error! Vy needed"
            
            for j in range(self.dim[1]-1):
                for i in range(self.dim[2]-1):
                    ip = i+1
                    l = i+j*self.dim[2]
                    DvxDy[i,j] = (vx[k,j+1,i] - vx[k,j,i])/(self.Y[j+1]-self.Y[j])
                    DvyDx[i,j] = (vx[k,j,i+1] - vx[k,j,i])/(self.X[i+1]-self.X[i])
                    vorticity[l] = (DvxDy[i,j] - DvyDx[i,j])#*(DvxDy[i,j] + DvyDx[i,j])
                    vorticity[l] /= rho[k,j,i]

        self.vortencity = Field([1,self.dim[1],self.dim[2]],self.dom, data = vorticity)

    def __open_fields(self, n):
        try:
            filename = "Density" + "{:06d}".format(n) + '.dat'
            data = np.fromfile(filename, "d")
            self.rho = Field(self.dim, self.dom, data = data)
        except:
            print filename, "not found"
        try:
            filename = "Vx" + "{:06d}".format(n) + '.dat'
            data = np.fromfile(filename, "d")
            self.vx = Field(self.dim, self.dom, data = data)
        except:
            print filename, "not found"
        try:
            filename = "Vy" + "{:06d}".format(n) + '.dat'
            data = np.fromfile(filename, "d")
            self.vy = Field(self.dim, self.dom, data = data)
        except:
            print filename, "not found"
        try:
            filename = "Vz" + "{:06d}".format(n) + '.dat'
            data = np.fromfile(filename, "d")
            self.vz = Field(self.dim, self.dom, data = data)
        except:
            print filename, "not found"
        try:
            filename = "Energy" + "{:06d}".format(n) + '.dat'
            data = np.fromfile(filename, "d")
            self.e = Field(self.dim, self.dom, data = data)
        except:
            print filename, "not found"
        try:
            filename = "Bx" + "{:06d}".format(n) + '.dat'
            data = np.fromfile(filename, "d")
            self.bx = Field(self.dim, self.dom, data = data)
        except:
            print filename, "not found"
        try:
            filename = "By" + "{:06d}".format(n) + '.dat'
            data = np.fromfile(filename, "d")
            self.by = Field(self.dim, self.dom, data = data)
        except:
            print filename, "not found"
        try:
            filename = "Bz" + "{:06d}".format(n) + '.dat'
            data = np.fromfile(filename, "d")
            self.bz = Field(self.dim, self.dom, data = data)
        except:
            print filename, "not found"

        self.n = n

    def __read_info(self):
        try:
            dims = open("dims.dat",'r')
        except:
            print 'Error!!!', 'dims.dat cannot be opened.'

        DIMS = dims.readlines() #Read all dims file
        dims.close()
        DIMS = DIMS[1].split()  #Splitting line 1 of dims
        dim1 = int(DIMS[6])
        dim2 = int(DIMS[7])
        dim3 = int(DIMS[8])
        
        self.dim = [dim3,dim2,dim1] # k,j,i
            
        if(dim1!=1):
            try:
                domain_x = open("domain_x.dat",'r')
                DOMAIN_X = domain_x.readlines()
                self.X = np.ndarray([dim1], dtype = float)
                for i in range(dim1):
                    self.X[i] = float(DOMAIN_X[i])
            except:
                print 'Error!!!', 'domain_x cannot be opened.'
        if(dim2!=1):
            try:
                domain_y = open("domain_y.dat",'r')
                DOMAIN_Y = domain_y.readlines()
                domain_y.close()
                self.Y = np.ndarray([dim2], dtype = float)
                for i in range(4,dim2+4):
                    self.Y[i-4] = float(DOMAIN_Y[i])
            except:
                print 'Error!!!', 'domain_y cannot be opened.'
        if(dim3!=1):
            try:
                domain_z = open("domain_z.dat",'r')
                DOMAIN_Z = domain_z.readlines()
                domain_z.close()
                self.Z = np.ndarray([dim3], dtype = float)
                for i in range(4,dim3+4):
                    self.Z[i-4] = float(DOMAIN_Z[i])
            except:
                print 'Error!!!', 'domain_z cannot be opened.'

        xmin = float(DIMS[0])
        xmax = float(DIMS[1])
        ymin = float(DIMS[2])
        ymax = float(DIMS[3])
        zmin = float(DIMS[4])
        zmax = float(DIMS[5])

        self.dom = [xmin, xmax,\
                    ymin, ymax,\
                    zmin, zmax]

    def __init__(self,n):
        self.__read_info()
        self.__open_fields(n)
