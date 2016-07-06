import matplotlib.pyplot as plt
import numpy as np
import sys

n = int(sys.argv[1])

input_file = "Density" + "{:06d}".format(n) + '.dat'
density = np.fromfile(input_file, "d")

resolucion = 32
hj = 12
hi = 64
nsec = hi*resolucion
nrad = hj*resolucion
xmin = -0.32
xmax = 0.32
ymin = 0.94
ymax = 1.06
h = 0.01
xx = 1.33
j = int((hj/2.0+xx)*nrad/hj)
rp = 1.0
phip = 0.0
sigma0 = 6.3661977237e-4
mth = 1e-6
mp = 1e-7

curve = np.ndarray([nsec])
x     = np.ndarray([nsec])

for i in range(nsec):
    curve[i] = (mth/mp)*((density[i+j*nsec]-sigma0)/sigma0)/np.sqrt((np.abs(xx)))
    phi = xmin + (xmax-xmin)/nsec*i
    x[i] = (rp*(phip-phi))/h-3*(np.abs(xx))*(np.abs(xx))/4.0;

fig = plt.figure()
axis = fig.add_subplot(111)
plt.plot(x,curve)
plt.show()
plt.savefig("snap.png", dpi = 200, format = 'png')
plt.close()
