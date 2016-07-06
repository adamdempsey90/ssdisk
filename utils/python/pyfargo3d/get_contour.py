import numpy as np
import matplotlib.pyplot as plt
import matplotlib._cntr as cntr

def get_contour(x,y,f,level):
    """
    (C) Pablo Benitez Llambay 2013
    This routine extracts a contour line, at the level 
    specified by 'level'. It depends on routines defined in
    the file contour.py inside of matplotlib package.
    Tested with matplotlib 1.1.1
    
    Parameters:
    -----------
    
    x,y: array_like
         2-D regular domain array.
    z:   array_like
         2-D regular value of the function z=f(x,y)
    level: float
           level of the desired contour
    
    Output:
    -------

    list of contours in the form [(x1,y1),(x2,y2),...(xn,yn)] 
    where xi,yi are 1-d array_like. If there are not contours,
    the function returns None

    Example:
    --------
    
    x = np.linspace(-1,1,128)
    X,Y = np.meshgrid(x,x)
    F = np.sin(X)**2+np.cos(Y)**2
    plt.imshow(F,origin='lower',extent=[-1,1,-1,1],interpolation='nearest')
    contour = get_contour(X,Y,F,0.9)
    for c in contour:
        plt.plot(c[0],c[1],'k')
    plt.show()    
    """

    c = cntr.Cntr(x,y,f)
    contours = c.trace(level)
    conts = []
    c = []
    for contour in contours:
        if contour.dtype != 'uint8':
            x = []
            y = []
            for point in contour:
                x.append(point[0])
                y.append(point[1])
            c.append((x,y))
    if len(c) == 0:
        return None
    else:
        return c

#Example
if __name__ == '__main__':
    x = np.linspace(-1,1,128)
    X,Y = np.meshgrid(x,x)
    F = np.sin(X)**2+np.cos(Y)**2
    plt.imshow(F,origin='lower',extent=[-1,1,-1,1],interpolation='nearest')
    contour = get_contour(X,Y,F,0.9)
    for c in contour:
        plt.plot(c[0],c[1],'k')
    plt.show()
