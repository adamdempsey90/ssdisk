from get_contour import get_contour
from intersect_curves import intersect_curves

def get_stagnation_points(v1, v2, x0, y0, w=0.1, plot=False, **karg):

    """
    Getting the zero velocity point for a staggered mesh
    Don't forget to use a corotating reference frame.

    Parameters: 
    ----------- 

        v1, v2:  Field class.
        x0, y0:  float
                 The center of the searching box.
        w = 0.1: float.  
                 Is the width of the target zone. Is
                 in units of the width of the main box
        plot = False
                 If it is true, you will see the zero velocity curves. For debuggin only.
                 **karg are the arguments supported by the matplotlib plot function.

     Output: 
     ------- 
        A list with all the stagnation points around
        the planet and contours [(s0x,s0y),...,(snx,sny)].

    """

    nx1 = v1.Params.nx
    ny1 = v1.Params.ny

    x1min = v1.xmin
    x1max = v1.xmax
    x2min = v2.xmin
    x2max = v2.xmax

    y1min = v1.ymin
    y1max = v1.ymax
    y2min = v2.ymin
    y2max = v2.ymax

    nx2 = v2.Params.nx
    ny2 = v2.Params.ny
    
    i1 = int((x0-x1min)/(x1max-x1min)*nx1)
    j1 = int((y0-y1min)/(y1max-y1min)*ny1)

    i2 = int((x0-x2min)/(x1max-x1min)*nx2)
    j2 = int((y0-y2min)/(y2max-y2min)*ny2)

    w *= 0.5

    w1xp  = int(i1+w*nx1)
    w1xm  = int(i1-w*nx1)
    w1yp  = int(j1+w*ny1)
    w1ym  = int(j1-w*ny1)

    w2xp  = int(i2+w*nx2)
    w2xm  = int(i2-w*nx2)
    w2yp  = int(j2+w*ny2)
    w2ym  = int(j2-w*ny2)

    X1, Y1 = np.meshgrid(v1.x,v1.y)
    contour1 = get_contour(X1[w1ym:w1yp,w1xm:w1xp],
                           Y1[w1ym:w1yp,w1xm:w1xp],
                           v1.data[w1ym:w1yp,w1xm:w1xp],0.0)

    X2, Y2 = np.meshgrid(v2.x,v2.y)
    contour2 = get_contour(X2[w2ym:w2yp,w2xm:w2xp],
                           Y2[w2ym:w2yp,w2xm:w2xp],
                           v2.data[w2ym:w2yp,w2xm:w2xp],0.0)
    
    if contour1 == None or contour2 == None:
        print "I can't find zero velocities. Try to play with w parameter..."
        return None

    if np.shape(contour1)[0] > 1 or np.shape(contour2)[0] > 1:
        print "Error: I found more than one contour of zero velocity," \
            " set a small w and try to avoid this problem."
        if plot:
            for c in contour1:
                plt.plot(c[0],c[1],'k')
            for c in contour2:
                plt.plot(c[0],c[1],'r')
        return None

    x1 = np.array(contour1[0][0])
    y1 = np.array(contour1[0][1])
    x2 = np.array(contour2[0][0])
    y2 = np.array(contour2[0][1])

    stagnations = intersect_curves(x1,y1,x2,y2,1)
    
    if plot:
        for c1 in contour1:
            plt.plot(c1[0],c1[1],**karg)
            for c2 in contour2:
                plt.plot(c2[0],c2[1],**karg)
        for s in stagnations:
            plt.plot(s[0],s[1],'ko',**karg)
        plt.show()

    return stagnations#,contour1,contour1]
