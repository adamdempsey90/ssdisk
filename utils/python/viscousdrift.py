from scipy.special import iv
import numpy as np
import matplotlib.pyplot as plt

"""
Note: tau = 12*\nu*t
The output of this script should be compared with the 
output of the VISCOUSDRIFT test.
"""

n = 5; nr = 1000

t = np.arange(0.01,0.01*n,0.01)
x = np.arange(0.5+0.5*1.0/nr,1.5+0.5*1.0/nr,1.0/nr)

t0 = 0.01
i = 0
for tau in t:
#    rho = np.fromfile("fargo3d/src/trunk/viscousdrift/Density{:06d}.dat".format(i))#.reshape([100,100])[:,50]
#    plt.plot(rho,'r-')
    y = 1.0/tau*x**(-0.25)*iv(0.25,2.0*x/tau)*np.exp(-(1.0+x**2)/tau)
    plt.plot(y,'b-')
    i += 1
plt.show()
