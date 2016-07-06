#!/usr/bin/python

import os
import os.path
import numpy as np
import matplotlib.pyplot as plt

def adjustdir (dir='.',fo=0,jo=0):
    res = dir.rstrip('\n')
    if (fo == 1):
        cmd=os.popen("tail -n 1 ~/.fargo3drc/lastout")
        res=cmd.read()
        cmd.close()
    if (jo == 1):
        cmd=os.popen("tail -n 1 ~/.jupiterrc/lastout")
        res=cmd.read()
        cmd.close()
    if (res != dir):
        print "Seeking data in "+res.rstrip('\n')
    return res.rstrip('\n')

def guesscode (dir='.'):
    dir = dir.rstrip('\n')
    testjup = dir+'/testendian.dat'
    testfar = dir+'/dimensions.dat'
    #We try to guess here whether we deal with FARGO3D or JUPITER
    if (os.path.exists(testjup)):
        return 'j' #This seems to be a JUPITER output directory
    if (os.path.exists(testfar)):
        return 'f' #This seems to be a FARGO3D output directory
    print "I do not recognize the type of output directory"
    return 'u'

def getreadfilename (dir='.',radix='dens'):
    code = guesscode(dir=dir)
    fargo   = {'dens':['gasdens'], 'vx':['gasvx'], 'vy':['gasvy'], 'vz':['gasvz'],
                'ener':['gasenergy'], 'erad':['energyrad'],
                'temp':['gasenergy','gasdens'],'tau':['Tau'],'stheat':['srad']}
    jupiter = {'dens':['gasdensity'], 'vx':['gasvelocity'], 'vy':['gasvelocity'],
               'vz':['gasvelocity'], 'ener':['gasenergy'], 'erad':['gaserad'],
               'temp':['gasenergy','gasdensity'],'tau':['gastau'],'stheat':['gasstheat']}
    comp = 1
    nbcomp = 1
    if ((code == 'j') and (radix == 'vx')):
        comp=1
        nbcomp = 3 # works with 3D version only
    if ((code == 'j') and (radix == 'vy')):
        comp=2
        nbcomp = 3 # works with 3D version only
    if ((code == 'j') and (radix == 'vz')):
        comp=3
        nbcomp = 3 # works with 3D version only
    if (code == 'f'):
        return [fargo[radix],comp,nbcomp]
    if (code == 'j'):
        return [jupiter[radix],comp,nbcomp]

def getreadpath (dir='.',nb=0,code='f'):
    if (code == 'f'):
        return dir+'/'
    if (code == 'j'):
        return dir+"/output%05d"%nb+'/'

def getvaluefromfile (string, file):
    cmd=os.popen("grep -i "+string+" "+file+"| tr '\t' '\n'|tail -n 1")
    a=cmd.read()
    cmd.close()
    return int(float(a))
    
def floatvaluefromfile (string, file):
    cmd=os.popen("grep -i "+string+" "+file+"| tr '\t' '\n'|tail -n 1")
    a=cmd.read()
    cmd.close()
    return float(a)
    
def getsize (dir='.'):
    code = guesscode (dir=dir)
    if (code == 'j'):
        nx = getvaluefromfile ("SIZE1", dir+'/params.dat')
        ny = getvaluefromfile ("SIZE2", dir+'/params.dat')
        nz = getvaluefromfile ("SIZE3", dir+'/params.dat')
    if (code == 'f'):
        nx = getvaluefromfile ("NX", dir+'/variables.par')
        ny = getvaluefromfile ("NY", dir+'/variables.par')
        nz = getvaluefromfile ("NZ", dir+'/variables.par')
    return [nx,ny,nz]

#low level function
def readfilename (name, dir='.', nb=0):
    code = guesscode (dir=dir)
    path = getreadpath (code=code,nb=nb,dir=dir)
    if (code == 'j'):
        suffix = str(nb)+'_0_0.dat'
    if (code == 'f'):
        suffix = str(nb)+'.dat'
    try:
        m = np.fromfile(path+name+suffix)
    except:
        if (code == 'f'):
            suffix = str(nb)+'_0.dat'
            m = np.fromfile(path+name+suffix)
    return m

#intermediate level function
def readfile (dir='.', nb=0, radix='dens'):
    namelist = getreadfilename(dir=dir,radix=radix)
    count = 0
    for nm in namelist[0]:
        m = readfilename (nm, dir=dir,nb=nb)
        size = m.size/namelist[2]
        m = m[(namelist[1]-1)*size:(namelist[1]*size)]
        if (count == 1):
            result = result/m
        else:
            result = m
        count = count+1
    return result

#high level functions
def r (dir='.',nb=0,radix='dens',fo=0,jo=0):
    dir = adjustdir (dir=dir,fo=fo,jo=jo)
    code = guesscode (dir=dir)
    if (code == 'f'):
        print ("Reading FARGO3D output")
    if (code == 'j'):
        print ("Reading JUPITER output")
    if (code == 'u'):
        print ("Unknown output format in directory "+dir)
        return [-1]
    m = readfile (dir=dir,nb=nb,radix=radix)
    s = getsize (dir=dir)
    print ("Simulation size is "+str(s[0])+" * "+str(s[1])+" * "+str(s[2]))
    try:
        m = m.reshape(s[2],s[1],s[0])
    except ValueError:
        if (code == 'f'):
            try:
                m = m.reshape(s[2]+6,s[1]+6,s[0])
                print ("Data size is "+str(s[0])+" * "+str(s[1]+6)+" * "+str(s[2]+6))
                print ("Ghost values are output")
            except:
                print "Data size mismatch"
                return [-1]
        if (code == 'j'):
            try:
                m = m.reshape(s[2]+4,s[1]+4,s[0]+4)
                print ("Data size is "+str(s[0]+4)+" * "+str(s[1]+4)+" * "+str(s[2]+4))
                print ("Ghost values are output")
            except:
                print "Data size mismatch"
                return [-1]
    return m

def v (m):
    plt.clf()
    cmap = 'gist_stern'
    print ("Min/Max of array: "+str(m.min())+" / "+str(m.max()))
    try:
        plt.imshow(m,origin='lower',interpolation='nearest',aspect='auto',cmap=cmap)
    except TypeError:
        try:
            plt.imshow(m[:,:,3],origin='lower',interpolation='nearest',aspect='auto',cmap=cmap)
            print ("Showing vertical slice near origin in x (outside of ghost)")
            print ("Min/Max of slice displayed: "+str(m[:,:,3].min())+" / "+str(m[:,:,3].max()))
        except IndexError:
            plt.imshow(m[:,:,0],origin='lower',interpolation='nearest',aspect='auto',cmap=cmap)
            print ("Showing vertical slice near origin in x")
            print ("Min/Max of slice displayed: "+str(m[:,:,0].min())+" / "+str(m[:,:,0].max()))
    plt.colorbar()

def rad (dir='.',fo=0,jo=0):
    dir = adjustdir (dir=dir,fo=fo,jo=jo)
    code = guesscode (dir=dir)
    if (code == 'j'):
        rmin = floatvaluefromfile ('RANGE2LOW',dir+'/params.dat')
        rmax = floatvaluefromfile ('RANGE2HIGH',dir+'/params.dat')
        nr   = getvaluefromfile   ('SIZE2',dir+'/params.dat')
    if (code == 'f'):
        rmin = floatvaluefromfile ('YMIN',dir+'/variables.par')
        rmax = floatvaluefromfile ('YMAX',dir+'/variables.par')
        nr   = getvaluefromfile   ('NY',  dir+'/variables.par')
    step = (rmax-rmin)/nr
    radii = np.linspace(rmin+.5*step,rmax-.5*step,num=nr,endpoint=True)
    return radii
        
def ar (dir='.',nb=0,gamma=1.43,fo=0,jo=0):
    dir = adjustdir (dir=dir,fo=fo,jo=jo)
    code = guesscode (dir=dir)
    t = r(radix='temp',dir=dir,nb=nb)
    cs_iso = np.sqrt((gamma-1.0)*t)
    # We need the value of the radius to get V_k
    vk = 1.0/np.sqrt(rad(dir=dir))
    cs = cs_iso[-1,:,:].mean(axis=1)
    if (vk.size != cs.size):
        if (code == 'f'):
            cs = cs[3:-3]
        if (code == 'j'):
            cs = cs[2:-2]
    return cs/vk

def ar2 (dir='.',nb=0,gamma=1.43,fo=0,jo=0):
    dir = adjustdir (dir=dir,fo=fo,jo=jo)
    code = guesscode (dir=dir)
    rho = r(radix='dens',dir=dir,nb=nb)
    rho = rho.mean(axis=2)
    radius = rad(dir=dir)
    radsize = (rho.shape)[1]
    asp_ratio = np.zeros(radsize)
    q=np.loadtxt(dir+"/domain_z.dat")
    dz=q[1]-q[0]
    z = np.arange((rho.shape)[0])*dz
    for i in np.arange(radsize):
        asp_ratio[i] = np.sqrt(np.pi/2.)*np.sum(z*rho[::-1,i])/np.sum(rho[:,i])
    return asp_ratio

