#!/usr/bin/env python
import mpl_toolkits.mplot3d as p3
import matplotlib.pyplot as p
from numpy import *


def generateSphere(center, radius, res = 10):
    u = linspace(0, 2 * pi, res)
    v = linspace(0, pi, res)
    x = radius * outer(cos(u), sin(v)) + center.x
    y = radius * outer(sin(u), sin(v)) + center.y
    z = radius * outer(ones(size(u)), cos(v)) + center.z
    return x,y,z

class Atom:
    def __init__(self, x=0,y=0,z=0, color='b'):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.color = color

    def getTuple():
        return (self.x, self.y, self.z)

class Cell:
    def __init__(self, dim=8, a=1, b=1, c=1):
        self.coords = []
        self.a = a
        self.b = b
        self.c = c
        self.dim = dim

        xoffset = a/float(2*dim)
        yoffset = b/float(2*dim)
        zoffset = c/float(2*dim)
        #c1 = '0.25'
        #c2 = '1.00'
        c1 = 'r'
        c2 = 'aquamarine'
        for x in arange(0,a,a/float(dim)):
            for y in arange(0,b,b/float(dim)):
                color = c1
                for z in arange(0,c,c/float(dim)):
                    if (color == c2): color = c1
                    else: color = c2
                    self.coords.append(Atom(x + xoffset,
                                            y + yoffset,
                                            z + zoffset,
                                            color))

    def printInfo(self):
        print "a =", self.a
        print "b =", self.b
        print "c =", self.c
        print ""
        print "Atoms:"
        for p in self.coords:
            print p.getTuple(), p.color

    def plotCellPoints(self, ax):
        xarray = []
        yarray = []
        zarray = []
        for p in self.coords:
            xarray.append(p.x)
            yarray.append(p.y)
            zarray.append(p.z)
        ax.scatter(xarray,yarray,zarray)

    def plotCellAtoms(self, ax, radius, res=10):
        for p in self.coords:
            x,y,z = generateSphere(p, radius, res)
            ax.plot_surface(x, y, z, rstride=1, cstride=1, color=p.color)

    def wrapToCell(self):
        for p in self.coords:
            while (p.x < 0.0):
                p.x += self.a
            while (p.x > self.a):
                p.x -= self.a
            while (p.y < 0.0):
                p.y += self.b
            while (p.y > self.b):
                p.y -= self.b
            while (p.z < 0.0):
                p.z += self.c
            while (p.z > self.c):
                p.z -= self.c

    def displace(self, mu = None, eta = None, rho = .10, wrap = True):
        if (mu == None): mu = self.dim
        if (eta == None): eta = self.dim
        xarray = []
        yarray = []
        zarray = []
        for p in self.coords:
            xarray.append(p.x)
            yarray.append(p.y)
            zarray.append(p.z)
        new_zarray = rho * ( 
              cos(mu * array(xarray)/float(self.a) * 2*pi)
              ) * (
              cos(eta * array(yarray)/float(self.b) * 2*pi)
              )
        for i in range(len(self.coords)):
            self.coords[i].z = zarray[i] + new_zarray[i]
        if (wrap): self.wrapToCell()
