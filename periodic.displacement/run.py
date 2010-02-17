#!/usr/bin/env python
import mpl_toolkits.mplot3d as p3
import matplotlib.pyplot as p
from numpy import *
from cell import *

# Parameters:
frames = 10
rho_min = 0
rho_max = .400
wrap = False
mu = 1
eta = 1
dim = 10
radius = 0.4*(1/float(dim))
res = 7
dpi = 300
ext = "png"

count = 0
fig = p.figure()
for i in arange(rho_min, rho_max, rho_max/float(frames)):
    print "Generating frame %d of %d"%(count+1,frames)
    cell = Cell(dim)
    ax = p3.Axes3D(fig)
    cell.displace(mu,eta,i,wrap)
    cell.plotCellAtoms(ax, radius, res)
#    p.suptitle("$\\rho$ = %.3f, $\mu$=$\eta$=1"%i)
    ax.set_xlim3d((0,1))
    ax.set_ylim3d((0,1))
    ax.set_zlim3d((0,1))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    p.savefig("movie/frame-%05d.%s"%(count+1,ext), dpi=dpi, bbox_inches="tight")
    p.savefig("movie/frame-%05d.%s"%(frames + (frames - count),ext), dpi=dpi, bbox_inches="tight")
    count += 1
    p.clf()
