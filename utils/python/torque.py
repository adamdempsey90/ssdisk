#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

def open_field(file, dim, axis = 0):
    field = np.reshape(np.fromfile(file), dim)
    return np.average(field, axis = axis)

def read_domain():
    t = np.array([float(line) for line in open('domain_x.dat', 'r')])
    r = np.array([float(line) for line in open('domain_y.dat', 'r')])
    rmed = 0.5*(r[5:-4] + r[4:-5]); tmed = 0.5*(t[1:] + t[:-1])
    return r,t,rmed,tmed

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

def surfaces(r,t):
    dr   = 0.5*(r[5:-4]**2-r[4:-5]**2) #Excluyo ghosts
    dphi = t[1:]-t[:-1]
    ds   = np.tensordot(dr,dphi,axes = 0) #diada r1*phi1 r1*phi2 ...
                                          #      r2*phi1 r2*
    return ds

def compute_mass(surf, rho):
    return rho*surf

def roche_lobe(q,rp):
    return rp*(0.33333*q)**0.333333

def cutoff(b,r):
    #http://arxiv.org/pdf/0908.1863.pdf
    #http://www.springerlink.com/content/h550752588284m88/
    #fermi-like cutoff
    cutoff = 1.0 + np.exp(-(r - b)/(b/10.0))
    return 1.0/cutoff

def compute_torque(rmed,tmed,xp,yp,mass,mp,h):
    ttmed,rrmed = np.meshgrid(tmed,rmed)
    xmed = rrmed*np.cos(ttmed) - xp; ymed = rrmed*np.sin(ttmed) - yp
    rp = np.sqrt(xp**2+yp**2)
    cut = cutoff(0.8*roche_lobe(mp,rp),np.sqrt(xmed*xmed+ymed*ymed))
    force_x = -mass*mp*xmed/np.power(xmed*xmed + ymed*ymed, 1.5)
    force_y = -mass*mp*ymed/np.power(xmed*xmed + ymed*ymed, 1.5)
    ref_torque = rp**(-3)*rp**2*mp**2*h**(-4)
    return (yp*force_x-xp*force_y)*cut/ref_torque

def planet_position(n):
    INPUT = open("planet0.dat",'r')
    lineas = INPUT.readlines()
    xp = (lineas[n].split())[1]
    yp = (lineas[n].split())[2]
    return float(xp),float(yp)

#====================================================================
#Variables
mp  = 6.0e-6
h = 0.05
n = 15
#====================================================================

dims = read_dims()
nx = int(dims[6]); ny = int(dims[7]); nz = int(dims[8]) #dimentions

xp,yp = planet_position(n)
r,t,rmed,tmed = read_domain()

surf = surfaces(r,t)

rho = open_field("Density"+"{:06d}".format(n) + ".dat",[nz,ny,nx])

mass = compute_mass(surf, rho)
average_mass   = np.average(mass,axis = 1)*2.0*np.pi/(t[-1]-t[0])

torque = compute_torque(rmed,tmed,xp,yp,mass,mp,h)
average_torque = np.average(torque,axis = 1)

torque_curve   = (average_torque/average_mass)

#plt.plot(rmed,torque_curve)
output = open("torque"+str(n)+".dat", "w")
for i in range(ny):
    linea = str((rmed[i]-1.0)/h) + "\t" +str(torque_curve[i]) + "\n"
    output.write(linea)
#    
plt.plot(rmed,torque_curve)
#    plt.ylim([-0.3,0.3])
#    #    plt.savefig("snap"+"{:06d}".format(n)+".png")
#    #    plt.clf()
plt.show()
