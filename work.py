#!/usr/bin/env python
from dataAnalysis import *

xs = []
ys = []
fs = []

for run in arange(30)+1:

    n = getCsvArray("run%d-results.txt"%run, 0) + 1
    m = getCsvArray("run%d-results.txt"%run, 3)

    x = []
    y = []

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
    scatter(x,y)
    plot(x,f)
    savefig("summary%d.png"%run)
    cla()

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

plot(x,minfs)
plot(x,maxfs)
plot(x,avefs)

minf = min(minfs)

def fitFunction2(c,x,y,get_y=False):
#    y_model = -minf*(e**(-c[0]*x+c[1])) + minf
    y_model = (e**(c[0]*x**c[1] + c[2]*x)) + minf
    if get_y == True: return y_model
    error = y-y_model
    return error

plotFitData(x,avefs,"Structure number", "Enthalpy (eV)", "Hartke Plot", "custom", plotData=False, color='g', fitFormat=':', 
            guess=array([-1,1,-1]), func=fitFunction2, fitLabel="regEqu", customFitLabel="$e^{(%%s)x^{%%s}+(%%s)x}+%s$"%(minf))

def fitFunction(c,x,y,get_y=False):
#    y_model = -minf*(e**(-c[0]*x+c[1])) + minf
    y_model = -minf*(e**(c[0]*x**c[1] + c[2]*x)) + minf
    if get_y == True: return y_model
    error = y-y_model
    return error

plotFitData(x,avefs,"Structure number", "Enthalpy (eV)", "Hartke Plot", "custom", plotData=False, color='r', fitFormat=':', 
            guess=array([-1,1,-1]), func=fitFunction, fitLabel="regEqu", customFitLabel="$%se^{(%%s)x^{%%s}+(%%s)x}+%s$"%(-minf,minf))

def fitFunction1(c,x,y,get_y=False):
#    y_model = -minf*(e**(-c[0]*x+c[1])) + minf
    y_model = -minf*(e**(c[0]*x**c[1])) + minf
    if get_y == True: return y_model
    error = y-y_model
    return error

plotFitData(x,avefs, "Structure number", "Enthalpy (eV)", "Hartke Plot", "custom", plotData=False, color='b', fitFormat=':', 
            guess=array([-1,1]), func=fitFunction1, fitLabel="regEqu", customFitLabel="$%se^{(%%s)x^{%%s}}+%s$"%(-minf,minf))

prop = matplotlib.font_manager.FontProperties(size=15)
legend(loc=1, prop = prop)
savefig("hartke.png")

cla()

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
