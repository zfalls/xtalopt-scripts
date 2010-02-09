#!/usr/bin/env python
from dataAnalysis import *
import sys,os

sys.argv.pop(0)

#
# Configuration:
#
# Average energy from a random analysis:
ave_energy = -622.52807
Emin = -636.806
# DPI for images
ext = "pdf"
dpi = 150

#
# Check to see if an update is really necessary
#
if sys.argv[0] != "force":
    if os.path.isfile("halflives.%s"%ext):
        lastRun = os.stat("halflives.%s"%ext).st_mtime
        uptodate = True
        for f in sys.argv:
            if os.stat(f).st_mtime > lastRun:
                uptodate = False
        if uptodate:
            print "Everything update -- not running"
            exit(0)
        else:
            print "Not up to date -- rechecking"
    else:
        print "No files from previous run. Generating data."
    
else:
    sys.argv.pop(0)
    print "Checking forced..."
    
#
# Read data in from files on command line
#

runs = []
halflives = []
rsqs = []

xs = []
ys = []
fs = []
nRuns = 0

for file in sys.argv:
    print "Reading file", file
    nRuns += 1
    n = getCsvArray(file, 0) + 1
    m = getCsvArray(file, 3)
    gen = getCsvArray(file, 1)

    x = []
    y = []

    # Remove bad values, e.g. errored xtals
    for i in range(len(n)):
        if gen[i] != 1:
            if abs(m[i]) > 10:
                x.append(n[i])
                y.append(m[i])

    t = y[0]
    f = []

    for i in y:
        if i < t:
            t = i
        f.append(t)
            
    xs.append(x)
    ys.append(y)
    fs.append(f)

    
#
# Prepare the data for a Hartke plot
#
for numruns in arange(len(xs))+1:
    print "Generating data for %d runs..."%(numruns)
    newxs = xs[:numruns]
    avefs = []

    minlen = len(newxs[0])
    for i in newxs:
        if len(i) < minlen: minlen = len(i)

    x = arange(minlen)+1
    for point in x:
        avef = 0
        count = 0
        for run in range(len(newxs)):
            tf = fs[run][point-1]
            avef += tf
            count += 1
        avef /= float(count)
        avefs.append(avef)

    E_0 = ave_energy-Emin

    def bestFitFunction(c,x,y,get_y=False):
        y_model = E_0*(e**(c[0]*x**c[1])) + Emin
        if get_y == True: return y_model
        error = y-y_model
        return error

    reg,err,rsq = plotFitData(x,avefs, "Structure number", "Enthalpy (eV)", "Hartke Plot", 
                              "custom", plotData=False, color='k', fitFormat=':', 
                              guess=array([-1,1]), func=bestFitFunction, fitLabel="regEqu", 
                              customFitLabel="$%se^{(%%s)x^{%%s}}+%s$"%(E_0,Emin))
    print "R^2 for hartke fit:",rsq

    # Solve for halflife with Newton-Raphson:
    tol = 1e-10
    guess = 100
    diff = 1e-5

    const = -log(1/2.0)
    def bestHalfLife(x):
        return reg[0]*x**reg[1] + const
    x = guess
    val = bestHalfLife(x)
    print E_0, Emin, const, reg, x, val
    while (abs(val) > tol):
        dx = (bestHalfLife(x+diff) - val)/diff
        x = x - val/dx
        val = bestHalfLife(x)
        print "Halflife: %.5f, value: %.6f, dx: %.6f"%(x,val,dx)

    y = bestFitFunction(reg,x,0,True)
    halflife = x
    calchalfE = y
    acthalfE = E_0/2.0+Emin
    print "At structure %.5f, the energy should be at %.5f (compare to %.5f)"%(x,
                                                                               calchalfE,
                                                                               acthalfE)

    halflives.append(halflife)
    rsqs.append(rsq)
    runs.append(numruns)

halflife = halflives[len(halflives)-1]
cla()
plot(runs,halflives, 'k-', label="Halflives")
gca().set_xlabel("Number of Runs")
gca().set_ylabel("Halflife")
gca().set_ylim([halflife-5,halflife+5])
legend(loc=3)
twinx()
plot(runs,rsqs, 'k:', label="$R^2$")
gca().set_ylabel("$R^2$ of Hartke plot")
legend(loc=4)
savefig("halflives.%s"%ext, dpi=dpi, bbox_inches="tight")
cla()
