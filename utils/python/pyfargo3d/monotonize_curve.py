import numpy as np
import matplotlib.pyplot as plt

"""
Splitting a curve in monotonically segments!
"""

class SplitCurve():
    def __init__(self,x,y):
        """
        x,y are general arrays, the curve is any curve!
        """
        self.x = x
        self.y = y
        self.segments = self.build_segments()

    def scan_segments(self):
        dif  = self.x[1:]-self.x[:-1] + 1e-16 #Avoiding symetry problems!
        sign = np.sign(dif[1:]*dif[:-1])

        index = (np.where(sign == -1)[0]+1).tolist()
        index.insert(0,0)
        index.append(len(self.x))

        index = np.array(index)
        half_index = (0.5*(index[1:]+index[:-1])).astype(np.int)
        signs       = np.sign(dif[half_index])

        separators = []
        for i in range(len(index)-1):
            separators.append((index[i],index[i+1]))
        
        return np.array(separators),signs

    def build_segments(self):
        separators,signs =  self.scan_segments()
        segments = []
        for s,d in zip(separators,signs):
            if d == -1:
                sort = np.argsort(self.x[s[0]:s[1]+1])
                segments.append([self.x[s[0]:s[1]+1][sort],
                                self.y[s[0]:s[1]+1][sort]])
            else:
                segments.append([self.x[s[0]:s[1]+1],
                                self.y[s[0]:s[1]+1]])
        return segments


def monotonize_curve(x,y):
    """
    Funtion that make x-monotonycally increasing segments from
    a general curve.
    
    Parameters:
    -----------
    x,y: 1d array_like
         
    Output:
    -------
    list of segments in the form [(x1,y1),(x2,y2),...(xn,yn)] where
    xi,yi are 1d array_like with the coordinates of each segment.
    "xi is monotonically increasing."

    Example:
    --------

    r = np.linspace(0.3,1.3,1000)
    t = np.linspace(0,2*np.pi,1000)

    x = r*np.cos(5*t)
    y = r*np.sin(5*t)

    segments = monotonize_curve(x,y)

    """

    segments = SplitCurve(x,y).segments
    return segments

if __name__ == '__main__':

    r = np.linspace(0.3,1.3,1000)
    t = np.linspace(0,2*np.pi,1000)
    x = r*np.cos(5*t)
    y = r*np.sin(5*t)

    segments = monotonize_curve(x,y)

    for s in segments:
        plt.plot(s[0],s[1],linewidth=5)
    
    plt.show()
