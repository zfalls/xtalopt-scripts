#!/usr/bin/env python
from dataAnalysis import *
import sys

#
# Configuration:
#
# Average energy from a random analysis:
ave_energy = -622.52807
Emin = -636.806
# DPI for images
dpi = 150

#
# Read data in from files on command line
#
xs = []
ys = []
fs = []
killed = 0
duplicate = 0
optimized = 0
nRuns = 0

sys.argv.pop(0)

for file in sys.argv:
    print "Reading file", file
    nRuns += 1
    n = getCsvArray(file, 0) + 1
    m = getCsvArray(file, 3)
    gen = getCsvArray(file, 1)
    status = getCsvStrings(file, 7, True)

    x = []
    y = []

    # Remove bad values, e.g. errored xtals
    for i in range(len(n)):
        if gen[i] != 1:
            if abs(m[i]) > 10:
                x.append(n[i])
                y.append(m[i])
            if status[i] == "Killed": killed+=1
            if status[i] == "Duplicate": duplicate+=1
            if status[i] == "Optimized": optimized+=1

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
#    scatter(x,y)
#    plotFitData(x,f, "Structure number", "Enthalpy (eV)", "Summary for %s"%file, "connect", plotData = False, color='r')
#    savefig("%s-summary.png"%file, dpi=dpi)
#    cla()

#
# Prepare the data for a Hartke plot
#
minfs = []
maxfs = []
avefs = []
stdfs = []

minlen = len(xs[0])
for i in xs:
    if len(i) < minlen: minlen = len(i)

x = arange(minlen)+1
for point in x:
    avef = 0
    stdf = 0
    minf = 10000000
    maxf = -10000000
    count = 0
    for run in range(len(xs)):
        tf = fs[run][point-1]
        avef += tf
        count += 1
        if minf > tf: minf = tf
        if maxf < tf: maxf = tf
    avef /= float(count)
    minfs.append(minf)
    maxfs.append(maxf)
    avefs.append(avef)
    for run in range(len(xs)):
        tf = fs[run][point-1]
        stdf += pow(tf - avef, 2)
    stdfs.append(sqrt(stdf/float(count)))

#
# Generate the error region around fs
#
higherr = []
lowerr = []
for i in range(len(avefs)):
    higherr.append(avefs[i] + stdfs[i])
    low = avefs[i] - stdfs[i]
    if low < Emin: low = Emin
    lowerr.append(low)

#
# Create a plot of the percent of runs completed by structure
#
percents = []
for point in range(minlen):
    done = 0
    total = 0
    for run in fs:
        if run[point] < Emin + 1e-3:
            done += 1
        total += 1
    percents.append(done/float(total)*100)

plotFitData(x, percents, "Structure number", "Percent of runs with lowest energy structure", "Percent complete by structure", 
            plotData=False, regType="connect")
savefig("percents.png", dpi=dpi)
cla()
percentDone = percents[len(percents)-1]

#
# Generate the Hartke plot
#
plot(x,maxfs, color='b', label="Worst-best structure")
plot(x,avefs, color='r', label="Average best structure")
plot(x,minfs, color='g', label="Best-best structure")
#fill_between(x, lowerr, higherr, alpha=0.2, color='r')

# Fit average f function
if Emin > min(minfs):
    print "Warning: Specified Emin (%f) > found Emin (%f)"%(Emin, min(minfs))

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

#
# Fit error function
#
higherr = []
lowerr = []
for i in range(len(avefs)):
    higherr.append(bestFitFunction(reg + err, x[i], 0, True))
    low = bestFitFunction(reg - err, x[i], 0, True)
    if low < Emin: low = Emin
    lowerr.append(low)
#fill_between(x, lowerr, higherr, alpha=0.2, color='k')

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

plot(halflife,calchalfE,'x', color='k', label="Halflife of average (x=%.5f)"%x)

prop = matplotlib.font_manager.FontProperties(size=10)
legend(loc=1, prop = prop)
savefig("hartke.png", dpi=dpi)
cla()

#
# Write summary file:
#
str = ""
str += "nRuns: 		%d"%nRuns + "\n"
str += "fit-a: 		%s"%numErrorText(reg[0], err[0], 5) + "\n"
str += "fit-b: 		%s"%numErrorText(reg[1], err[1], 5) + "\n"
str += "rsq: 		%f"%rsq + "\n"
str += "halflife:	%0.1f"%halflife + "\n"
str += "percentDone:	%0.1f"%percentDone + "\n"
str += "calchalfE: 	%f"%calchalfE + "\n"
str += "acthalfE:	%f"%acthalfE + "\n"
str += "ave_energy:	%f"%ave_energy + "\n"
str += "Emin:		%f"%Emin + "\n"
str += "killed:		%0.1f (%d)"%(killed/float(duplicate+optimized+killed)*100,killed) + "\n"
str += "duplicate:	%0.1f (%d)"%(duplicate/float(duplicate+optimized)*100,duplicate) + "\n"
str += "optimized:	%0.1f (%d)"%(optimized/float(duplicate+optimized+killed)*100,optimized) + "\n"

f = open("summary", 'w')
f.write(str)
f.close()
print "Summary file written:\n" + str