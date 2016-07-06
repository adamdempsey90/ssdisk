#!/usr/bin/python

import f3djup as f
import matplotlib.pyplot as plt

d1=f.r(nb=0,dir='test1',radix='tau')
d2=f.r(nb=0,dir='test2',radix='tau')

f.v(d2/d1)

plt.show()


