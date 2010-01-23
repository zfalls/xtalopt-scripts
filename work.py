#!/usr/bin/env python
from dataAnalysis import *
import sys

#
# Average energy from a random analysis:
#
ave_energy = -622.52807

#
# Read data in from files on command line
#
xs = []
ys = []
fs = []

sys.argv.pop(0)

for file in sys.argv:
    print "Reading file", file
    n = getCsvArray(file, 0) + 1
    m = getCsvArray(file, 3)

    x = []
    y = []

    # Remove bad values, e.g. errored xtals
    for i in range(len(n)):
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
    # Create a summary plot for each data set
    scatter(x,y)
    plot(x,f)
    savefig("%s-summary.png"%file)
    cla()

#
# Prepare the data for a Hartke plot
#
minfs = []
maxfs = []
avefs = []

minlen = len(xs[0])
for i in xs:
    if len(i) < minlen: minlen = len(i)

x = arange(minlen)+1
for point in x:
    avef = 0
    minf = 10000000
    maxf = -10000000
    count = 0
    for run in range(len(xs)):
        tf = fs[run][point-1]
        avef += tf
        count += 1
        if minf > tf: minf = tf
        if maxf < tf: maxf = tf
    avef /= count
    minfs.append(minf)
    maxfs.append(maxf)
    avefs.append(avef)

#
# Create a plot of the percent of runs completed by structure
#
percents = []
for point in range(minlen):
    done = 0
    total = 0
    for run in fs:
        if run[point] - minf < 1e-3:
            done += 1
        total += 1
    percents.append(done/float(total)*100)

plotFitData(x, percents, "Structure number", "Percent of runs with lowest energy structure", "Percent complete by structure", 
            plotData=False, regType="connect")
savefig("percents.png")
cla()

#
# Generate the Hartke plot
#
plot(x,maxfs, color='b', label="Worst-best structure")
plot(x,avefs, color='r', label="Average best structure")
plot(x,minfs, color='g', label="Best-best structure")

# Fit average f function
E_min = min(minfs)
E_0 = ave_energy-E_min
def bestFitFunction(c,x,y,get_y=False):
    y_model = E_0*(e**(c[0]*x**c[1] + c[2]*x)) + E_min
    if get_y == True: return y_model
    error = y-y_model
    return error

reg,err,rsq = plotFitData(x,avefs, "Structure number", "Enthalpy (eV)", "Hartke Plot", 
                          "custom", plotData=False, color='r', fitFormat=':', 
                          guess=array([-1,1,1]), func=bestFitFunction, fitLabel="regEqu", 
                          customFitLabel="$%se^{(%%s)x^{%%s}+(%%s)x}+%s$"%(E_0,E_min))

# Solve for halflife with Newton-Raphson:
tol = 1e-10
guess = 100
diff = 1e-5

const = -log(1/2.0)
def bestHalfLife(x):
    return reg[0]*x**reg[1] + reg[2]*x + const
x = guess
val = bestHalfLife(x)
print E_0, E_min, const, reg, x, val
while (abs(val) > tol):
    dx = (bestHalfLife(x+diff) - val)/diff
    x = x - val/dx
    val = bestHalfLife(x)
    print "Halflife: %.5f, value: %.6f, dx: %.6f"%(x,val,dx)

y = bestFitFunction(reg,x,0,True)
print "At structure %.5f, the energy should be at %.5f (compare to %.5f)"%(x,
                                                                           y,
                                                                           E_0/2.0+E_min)

plot(x,y,'x', color='k', label="Halflife of average (x=%.5f)"%x)

prop = matplotlib.font_manager.FontProperties(size=10)
legend(loc=1, prop = prop)
savefig("hartke.png", dpi=150)
cla()
