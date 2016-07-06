#!/usr/bin/env python
import os
from sys import argv
import multiprocessing
import matplotlib
matplotlib.use('Agg')
#save_images(q,filelist,dirname=imgdir,fnamebase=q,imgext=ext);
import numpy as np
import matplotlib.pyplot as plt
from numpy import *
from matplotlib.pyplot import *

util_dir = '/projects/b1002/amd616/fargo3d/utils/python/'
execfile(util_dir + 'advanced.py')
directory = argv[1]
n = int(argv[2])
proc = int(argv[3])

irange = np.arange(n)
args = [(i,j,directory) for i,j in enumerate(irange)]

def save_img((i,j,directory)):
    fig=plt.figure()
    ax = fig.add_subplot(111)
    dens=Field('gasdens{0:d}.dat'.format(i),directory=directory)
    dens.plot(ax=ax,cartesian=True,log=True)
    fig.savefig(directory + 'plots/dens%04d.png'%j)
    plt.close(fig)


if proc > 1:
    pool = multiprocessing.Pool(proc)
    print 'Starting'
    pool.map(save_img,args)
    print 'Finished'
else:
    for x in args:
        save_img(x)

