import numpy as np
import sys
from matplotlib.mlab import griddata
import multiprocessing as mp
#import visit_writer as vw


def write_vtk_file((n,option)):
    print 'Converting number %d'%n
    ifile_dens = "gasdens" + "{:d}".format(n) + '.dat'
    ifile_ener = "gasenergy" + "{:d}".format(n) + '.dat'
    ifile_mass= "gasmass" + "{:d}".format(n) + '.dat'

#    ifile_vx   = "gasvx" + "{:d}".format(n) + '.dat'
#    ifile_vy   = "gasvy" + "{:d}".format(n) + '.dat'
#    ifile_vz   = "gasvz" + "{:d}".format(n) + '.dat'

    density = np.fromfile(ifile_dens, "d")
    energy  = np.fromfile(ifile_ener, "d")
    mass = density.copy()
#    vx      = np.fromfile(ifile_vx  , "d")
#    vy      = np.fromfile(ifile_vy  , "d")
#    vz      = np.fromfile(ifile_vz  , "d")




    NGHY = 3
    NGHZ = 3

    x_arr = np.loadtxt('domain_x.dat')
    y_arr = np.loadtxt('domain_y.dat')[NGHY:-NGHY]
    z_arr = np.loadtxt('domain_z.dat')[NGHZ:-NGHZ]

    dx = x_arr[1]-x_arr[0];

    Nx = len(x_arr)-1
    Ny = len(y_arr)-1
    Nz = len(z_arr)-1


    print Nx,Ny,Nz
    print '%e total points'%(Nx*Ny*Nz)


    mesh  = np.ndarray([3*Nx*Ny*Nz], dtype = float)
    var1  = np.ndarray([Nx*Ny*Nz], dtype = float)


    l = 0
    stride = Nx*Ny
    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                if(option == 0):
                    mass[i + Nx*j +stride*k] *= (1./3.)*dx*(np.cos(z_arr[k+1])-np.cos(z_arr[k]))*(y_arr[j+1]**3  - y_arr[j]**3 )
                    mesh[l] = np.cos(x_arr[i])*np.sin(z_arr[k])*y_arr[j]
                    mesh[l+1] = np.sin(x_arr[i])*np.sin(z_arr[k])*y_arr[j]
                    mesh[l+2] = y_arr[j]*np.cos(z_arr[k])
                if(option == 1):
                    mass[i + Nx*j +stride*k] *= .5*dx*(y_arr[j+1]**2  - y_arr[j]**2 )
                    mesh[l] = np.cos(x_arr[i])*y_arr[j];
                    mesh[l+1] = np.sin(x_arr[i])*y_arr[j];
                    mesh[l+2] = z_arr[k];
                l+=3;
    mass.tofile(ifile_mass)

    #Mesh construction--------------------------------------
    connectivity = [] #Unstructured Mesh!, avoid periodic problem!
    stride = Nx*Ny
    for k in range(Nz-1):
        for j in range(Ny-1):
            for i in range(Nx):
                l = i+j*Nx+k*stride
                if(i<Nx-1):
                    lxp = l+1
                else:
                    lxp = l-(Nx-1) #WARNINR Nx-1?
                connectivity += [[12,
                    l, #  i,j,k
                    lxp, # i+1, j, k
                    lxp+Nx,  # i+1, j+1, k
                    l+Nx, # i, j+1, k
                    l+stride, # i, j, k+1
                    lxp+stride, # i+1, j, k+1
                    lxp+Nx+stride, #  i+1, j+1, k+1
                    l+Nx+stride]] # i, j+1, k+1
    #Mesh construction--------------------------------------
    #print connectivity[0], connectivity[Nx-1]

    dims = (Nx,Ny,Nz)
    var = (("Density", 1, 1, tuple(density)),
           ("Energy" , 1, 1, tuple(energy )),
           ("Mass", 1, 1, tuple(mass)))
           #("Vtheta" , 1, 1, tuple(vx     )),
           #("Vrad"   , 1, 1, tuple(vy     )),
           #("Vpolar" , 1, 1, tuple(vz     ))) #ZONAL

    #vw.WriteCurvilinearMesh("vwcurv3d.vtk", 1, dims, tuple(mesh), var)
    vw.WriteUnstructuredMesh("vwucd3d{:d}.vtk".format(n), 1, tuple(mesh),connectivity, var)


if __name__ == "__main__":
    try:
        import visit_writer as vw
    except ImportError:
        sys.path.append('/home/amd616/lib/python')
        import visit_writer as vw

    nstart = int(sys.argv[1])  #input from term
    nend = int(sys.argv[2]) + 1  #input from term
    option = int(sys.argv[3]) #1 = Cilindrical, 0 = Spherical

    print 'Converting outputs %d to %d'%(nstart,nend-1)
    if len(sys.argv) > 3:
        numprocs = int(sys.argv[4])
        print 'Using %d procs'%numprocs
        if numprocs > 1:
            p = mp.Pool(numprocs)
            args = [(i,option) for i in range(nstart,nend)]
            p.map(write_vtk_file,args)
        else:
            args = [(i,option) for i in range(nstart,nend)]
            for a in args:
                write_vtk_file(a)

    else:
        args = [(i,option) for i in range(nstart,nend)]
        for a in args:
            write_vtk_file(a)

