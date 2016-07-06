import numpy as np
import matplotlib.pyplot as plt
from intersect1d import intersect1d
from monotonize_curve import monotonize_curve

def intersect_curves(x1,y1,x2,y2,n=1):

    if type(x1) == list:
        x1 = np.array(x1)
    if type(y1) == list:
        y1 = np.array(y1)
    if type(x2) == list:
        x2 = np.array(x2)
    if type(y2) == list:
        y2 = np.array(y2)

    mono1 = monotonize_curve(x1,y1)
    mono2 = monotonize_curve(x2,y2)

    common_points = []
    for segment1 in mono1:
        for segment2 in mono2:
            intersection = intersect1d(segment1[0],segment1[1],
                                       segment2[0],segment2[1],n)
            if intersection != None:
                for point in intersection:
                    common_points.append(point)
    return common_points

if __name__ == '__main__':

    n = 1000
    
    t1 = np.linspace(0,2*np.pi,n)*10
    r1 = np.linspace(0.1,1.0,n)
    t2 = np.linspace(0,2*np.pi,n)*5
    r2 = np.linspace(0.05,1.0,n)
    
    x1 = r1*np.cos(t1)
    y1 = r1*np.sin(t1)

    x2 = r2*np.cos(t2)+0.3
    y2 = r2*np.sin(t2)
    
    plt.plot(x1,y1,x2,y2)
    
    intersections = intersect_curves(x1,y1,x2,y2)
    
    for point in intersections:
        plt.plot(point[0],point[1],'ro')
    plt.show()
