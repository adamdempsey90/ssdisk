#!/usr/bin/env  python
import sys
import os

"""
Script for make a video with plot2dfield.gpl
call form : python plot2dfields.py <nx> <ny> <n_outputs>
"""

nx = sys.argv[1]
ny = sys.argv[2]
n_outputs = sys.argv[3]

n_outputs = int(n_outputs)

for i in range(n_outputs):
    gnuplot_line = 'gnuplot -e "m=' + str(i) + \
        '; nx=' + nx + '; ny=' + ny +\
        '" plot2dfields.gpl'
    print gnuplot_line
    os.system(gnuplot_line)
