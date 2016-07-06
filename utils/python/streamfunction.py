import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate

def streamfunction(u,v,dx,dy):
    term1 = np.ndarray([u.shape[0],u.shape[1]])
    term2 = np.ndarray([u.shape[0],u.shape[1]])
    for i in range(u.shape[1]):
        for j in range(u.shape[0]):
#integral en x para y0 fijo
            term2[j,i] = integrate.trapz(v[j,0:i],dx=dx) 
            term1[j,i] = integrate.trapz(u[0:j,i],dx=dy)
    print np.shape(term1-term2)
    return term1 - term2
#las curvas de nivel de esta funcion son las lineas de flujo....
