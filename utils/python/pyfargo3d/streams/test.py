import numpy as np
import matplotlib.pyplot as plt
import stream as s
import pyfargo3d as f

def stream(S,x0,y0,nmax=10000,sign=1):
    """
    A SnapShot Wrapper for the stream module.
    """

    lista = s.stream(S.gasvx.data, S.gasvx.x, S.gasvx.y,
                     S.gasvy.data, S.gasvy.x, S.gasvy.y,
                     x0,y0,nmax,sign)

    return lista

S = f.SnapShot("/home/pablo/Tesis/data/guilet_base/",50)
print "\nExecuting Stream\n"

for i in range(1000):
    x,y = stream(S,-2.0,0.6+0.8*i/100,sign=1)
    plt.plot(x,y,'k')
    x,y = stream(S,-2.0,0.6+0.8*i/100,sign=-1)
    plt.plot(x,y,'k')

S.gasdens.plot()
plt.show()
