#!/usr/bin/env python
import sys
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
import sys
#from matplotlib.ticker import AutoMinorLocator #for minor ticks

def open_field(file, dim, axis = 0):
    field = np.reshape(np.fromfile(file), dim)
    return np.average(field, axis = axis)

def read_dims():
    """
    Reading dims.dat function.
    Output: array with elements -> 
    XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,NX,NY,NZ
    """
    try:
        DIMS = open("dims.dat",'r')
    except:
        print "Error reading dims.dat"
        exit(10)
    for line in DIMS:
        dims = line.split()
    return dims

#Mesh parameters------------------------
dims = read_dims() #Reading dims.dat
xmin = float(dims[0]); xmax = float(dims[1]) #ranges
ymin = float(dims[2]); ymax = float(dims[3])
zmin = float(dims[4]); zmax = float(dims[5])
nx = int(dims[6]); ny = int(dims[7]); nz = int(dims[8]) #dimentions
#Disk parameters------------------------
h = 0.05 #Disk aspect ratio
hi = (xmax-xmin)/h #azimutal extention in units of h
hj = (ymax-ymin)/h #radial extention in units of h
mth = h**3 #Thermal mass
#Planet parameters----------------------
rp = 1.0 #radial coord of the planet
phip = 0.0 #azimutal coord of the planet
mp = 6.0e-6 #planet mass
#Wake parameters------------------------
xx = 1.33 # Slice radius in units of h
jp = int((hj/2.0+xx)*ny/hj) #radial positive index of slice
jm = int((hj/2.0-xx)*ny/hj) #radial negative index of slice
#Input parameters-----------------------
number = 15
input = "Density"+"{:06d}".format(number) + ".dat"
Input = open_field(input,[nz,ny,nx])
print Input.shape
#Input = np.fromfile(input)
#Input = Input.reshape([nz,ny,nx])
#Output parameters----------------------
output = "wake"+"{:06d}".format(number)+".dat"
Output = open(output,"w")

#latex package
rc('text', usetex=True)
rc('font', family='serif')

#Averaging Density (Obtaining Sigma0)
Sigma0 = np.average(Input)

curve_p = (mth/mp)*(Input[jp,:]-Sigma0)/Sigma0
curve_p /= np.sqrt(np.abs(xx))
curve_m = (mth/mp)*(Input[jm,:]-Sigma0)/Sigma0
curve_m /= np.sqrt(np.abs(xx))
phi = np.arange(xmin,xmax,(xmax-xmin)/nx)
phi_p = rp*(phip-phi)/h - 3*np.abs(xx)**2/4.0
phi_m = rp*(-phip+phi)/h - 3*np.abs(xx)**2/4.0

for i in range(np.size(phi_p)):
    line = str(phi_p[i]) + "\t" + \
        str(curve_p[i]) + "\t" + \
        str(curve_m[-i-1]) + "\n"
    if(np.abs(phi_p[i])<=6):
        Output.write(line)

plt.plot(phi_p,curve_p)
plt.plot(phi_m,curve_m)
plt.show()


#for i in range(nsec):
#    curve_p[i] = (mth/mp)*((shear[i+jp*nsec]-sigma0)/sigma0)/np.sqrt((np.abs(xx)))
#    curve_m[i] = (mth/mp)*((shear[i+jm*nsec]-sigma0)/sigma0)/np.sqrt((np.abs(xx)))
#    phi = xmin + (xmax-xmin)/nsec*i
#    x_p[i] = (rp*(phip-phi))/h-3*(np.abs(xx))*(np.abs(xx))/4.0;
#    x_m[i] = (rp*(phi-phip))/h-3*(np.abs(xx))*(np.abs(xx))/4.0;
#
#fig = plt.figure()
#axis = fig.add_subplot(111)
#
#xlabel = axis.set_xlabel("$(y/h)-3(x/h)^2/4$")
#ylabel = axis.set_ylabel('$(M_{th}/M_p)(\delta\Sigma/\Sigma_0)/(x/h)^{1/2}$')


#rc('text', usetex=False)
#axis.plot(x_p,curve_p_shear,    label = 'shear_p'   , lw = 0.7)
#axis.plot(x_m,curve_m_shear,    label = 'shear_m'   , lw = 0.7)
#axis.plot(x_p,curve_p_nonshear, label = 'nonshear_p', lw = 0.7)
#axis.plot(x_m,curve_m_nonshear, label = 'nonshear_m', lw = 0.7)
#axis.plot(x_p,curve_p_standard, label = 'standard_p', lw = 0.7)
#axis.plot(x_m,curve_m_standard, label = 'standatd_m', lw = 0.7)

#xminorLocator   = AutoMinorLocator()
#axis.xaxis.set_minor_locator( xminorLocator )

#axis.minorticks_on()

#axis.tick_params(which='both', width=1)
#axis.tick_params(which='major', length=4)
#axis.tick_params(which='minor', length=2)

#axis.set_xlim(-6,6)
#axis.legend(loc = 'upper right')

#plt.show()
#plt.savefig("snap.pdf", dpi = 72, format = 'pdf')
#plt.close()
