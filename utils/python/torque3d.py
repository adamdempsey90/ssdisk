#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

def open_field(file, dim, axis = 0):
    field = np.reshape(np.fromfile(file), dim)
    return field

def read_domain():
    """
    t = theta
    r = radius
    pol = colatitude
    """
    t = np.array([float(line) for line in open('domain_x.dat', 'r')])
    r = np.array([float(line) for line in open('domain_y.dat', 'r')])
    pol = np.array([float(line) for line in open('domain_z.dat', 'r')])
    rmed = 0.5*(r[5:-4] + r[4:-5]); tmed = 0.5*(t[1:] + t[:-1])
    polmed = 0.5*(pol[5:-4] + pol[4:-5]);
    return r,t,pol,rmed,tmed,polmed

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

def volumes(r,pol,dtheta):
    dr   = dtheta*1./3.*(r[5:-4]**3-r[4:-5]**3) #Excluyo ghosts
    dpol = np.cos(pol[5:-4])-np.cos(pol[4:-5])
    print dr.shape, dpol.shape
    dvol   = np.tensordot(dpol,dr,axes = 0) #meridional slice volumes
    return dvol

def compute_mass(vol, rho):
    print rho[:,:,1].shape, vol.shape
#    exit()
    for i in xrange(np.size(np.size(rho[0,0,:]))):
        rho[:,:,i] *= vol
    return rho

def roche_lobe(q,rp):
    return rp*(0.33333*q)**0.333333

def cutoff(b,r):
    #http://arxiv.org/pdf/0908.1863.pdf
    #http://www.springerlink.com/content/h550752588284m88/
    #fermi-like cutoff
    cutoff = 1.0 + np.exp(-(r - b)/(b/10.0))
    return 1.0/cutoff

def compute_torque(rmed,tmed,polmed,xp,yp,zp,mass,mp,h):
    nx = np.size(tmed)
    ny = np.size(rmed)
    nz = np.size(polmed)
    ppolmed,rrmed,ttmed = np.mgrid[polmed[0]:polmed[-1]:complex(nz), \
                                       rmed[0]:rmed[-1]:complex(ny), \
                                       tmed[0]:tmed[-1]:complex(nx)]
    xmed = rrmed*np.cos(ttmed)*np.sin(ppolmed) - xp
    ymed = rrmed*np.sin(ttmed)*np.sin(ppolmed) - yp
    zmed = rrmed*np.cos(ppolmed) - zp
    print "Coordenadas listas"
    rp = np.sqrt(xp**2+yp**2+zp**2)
    cut = cutoff(0.8*roche_lobe(mp,rp),np.sqrt(xmed**2+ymed**2+zmed**2))
    force   = -mass*mp/np.power(xmed**2 + ymed**2 + zmed**2, 1.5)    
    ref_torque = rp**(-3)*rp**2*mp**2*h**(-4)
    print "Comenzando calculo de torques..."
    #torque_x = (zp*force*ymed-yp*force*zmed)
    print "Torque_x listo"
    #torque_y = (xp*force*zmed)#-zp*force*xmed)
    print "Torque_y listo"
    torque_z = -xp*force*ymed #(yp*force*xmed-xp*force*ymed)
    print "Torque_z listo"
    #Asumo que el torque en y es nulo por simetria (sum(torque_y) ~ 10^{-19})
    return np.sum(torque_z,axis=0)/ref_torque #Torque sumado sobre Z

def planet_position(n):
    INPUT = open("planet0.dat",'r')
    lineas = INPUT.readlines()
    xp = (lineas[n].split())[1]
    yp = (lineas[n].split())[2]
    zp = (lineas[n].split())[3]
    return float(xp),float(yp),float(zp)

#====================================================================
#Variables
mp  = 6.0e-6
h = 0.05
n = 15
#====================================================================

dims = read_dims()
nx = int(dims[6]); ny = int(dims[7]); nz = int(dims[8]) #dimentions

xp,yp,zp = planet_position(n)
r,t,pol,rmed,tmed,polmed = read_domain()

vol = volumes(r,pol,t[1]-t[0])

rho = open_field("Density000015.dat",[nz,ny,nx])

mass = compute_mass(vol, rho)
#Ahora sumo la masa en Z y la promedio azimutalmente. Si el disco es incompleto,
#extrapolo la masa a la contenida en un anillo.
average_mass   = np.average(np.sum(mass,axis = 0),axis=1)*2.0*np.pi/(t[-1]-t[0])
#Computo la suma del torque en Z, para cada radio y azimuth.
torque_z = compute_torque(rmed,tmed,polmed,xp,yp,zp,mass,mp,h)
#Computo el valor medio en azimuth del torque sumado 
average_torque = np.average(torque_z,axis = 1)
#El cociente de promedios cancela los numeros de celdas.
torque_curve   = (average_torque/average_mass)

output = open("torque"+str(n)+".dat", "w")
for i in range(ny):
    linea = str((rmed[i]-1.0)/h) + "\t" +str(torque_curve[i]) + "\n"
    output.write(linea)
    
plt.plot(rmed,torque_curve)
plt.ylim([-0.3,0.3])
##    #    plt.savefig("snap"+"{:06d}".format(n)+".png")
##    #    plt.clf()
plt.show()
