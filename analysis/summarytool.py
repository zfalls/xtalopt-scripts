#!/usr/bin/env python
from dataAnalysis import *
import sys,os

def generateSummary(path, files, force=False):

    #
    # Configuration:
    #
    # Average energy from a random analysis:
    ave_energy = -622.52807
    Emin = -636.806
    Etol = 1e-3
    # DPI for images
    dpi = 150
    # Extension for images
    ext = "pdf"

    #
    # Check to see if an update is really necessary
    #
    if not force:
        if os.path.isfile("%s/summary"%path):
            lastRun = os.stat("%s/summary"%path).st_mtime
            uptodate = True
            for f in files:
                if os.stat("%s/%s"%(path,f)).st_mtime > lastRun:
                    uptodate = False
            if uptodate:
                print "Everything update -- not running"
                return
            else:
                print "Not up to date -- rechecking"
        else:
            print "No files from previous run. Generating data."
    else:
        print "Checking forced..."
    
    #
    # Read data in from files on command line
    #
    xs = []
    ys = []
    fs = []
    firstDone = []
    killed = 0
    duplicate = 0
    optimized = 0
    nRuns = 0

    for f in files:
        file = "%s/%s"%(path,f)
        sys.stdout.write("Reading file %s%s\r"%(file," "*10))
        sys.stdout.flush()
        nRuns += 1
        n = getCsvArray(file, 0) + 1
        m = getCsvArray(file, 3)
        gen = getCsvArray(file, 1)
        status = getCsvStrings(file, 7, True)

        x = []
        y = []
        lowest = 0
        done = False

        # Remove bad values, e.g. errored xtals
        for i in range(len(n)):
            if gen[i] != 1:
                # If value is significantly lower than Emin, there's
                # likely a glitch in the calc. Also remove killed xtals.
                if status[i] != "Killed" and (abs(Emin) - abs(m[i]) > -Etol*10):
                    x.append(n[i])
                    y.append(m[i])
                if not done:
                    if fabs(Emin-m[i]) < Etol:
                        lowest = len(x)
                        done = True
                        optimized+=1
                    else:
                        if status[i] == "Killed": killed+=1
                        if status[i] == "Duplicate": duplicate+=1
                        if status[i] == "Optimized": optimized+=1

        firstDone.append(lowest)                
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
        #    savefig("%s-summary.%s"%(file,ext) dpi=dpi, bbox_inches="tight")
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
            if fabs(run[point] - Emin) < Etol:
                done += 1
            total += 1
        percents.append(done/float(total)*100)
    #plotFitData(x, percents, "Structure number", "Percent of runs with lowest energy structure", "Percent complete by structure", 
                #plotData=False, regType="connect")
    #savefig("%s/percents.%s"%(path,ext), dpi=dpi, bbox_inches="tight")
    #cla()
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

    xmax = max(x)

    reg,err,rsq = plotFitData(x,avefs, "Structure number", "Enthalpy (eV)", "Hartke Plot", 
                              "custom", plotData=False, color='k', fitFormat=':', 
                              guess=array([-1,1]), func=bestFitFunction, fitLabel="Fitted exponential", 
                              customFitLabel="$%se^{(%%s)x^{%%s}}+%s$"%(E_0,Emin))
    #print "R^2 for hartke fit:",rsq

    #
    # Fit error function
    #
    #higherr = []
    #lowerr = []
    #for i in range(len(avefs)):
        #higherr.append(bestFitFunction(reg + err, x[i], 0, True))
        #low = bestFitFunction(reg - err, x[i], 0, True)
        #if low < Emin: low = Emin
        #lowerr.append(low)
    #fill_between(x, lowerr, higherr, alpha=0.2, color='k')

    # Solve for halflife with Newton-Raphson:
    tol = 1e-10
    guess = 100
    diff = 1e-1

    const = -log(1/2.0)
    def bestHalfLife(x):
        return reg[0]*x**reg[1] + const
    x = guess
    val = bestHalfLife(x)
#    print E_0, Emin, const, reg, x, val
    while (abs(val) > tol):
        # Constrain x to be positive
        x = abs(x)
        dx = (bestHalfLife(x+diff) - val)/diff
        x = x - val/dx
        # Constrain x to be positive
        x = abs(x)
        val = bestHalfLife(x)
#        print "x = %f"%x
#        print "val = %f = %f"%(val,bestHalfLife(x))
#        print "val equ: %f*%f**%f+%f = %f"%(reg[0],x,reg[1], const, reg[0]*x**reg[1]+const)
#	print "Halflife: %.5f, value: %.6f, dx: %.6f"%(x,val,dx)

    y = bestFitFunction(reg,x,0,True)
    halflife = x
    calchalfE = y
    acthalfE = E_0/2.0+Emin
    #print "At structure %.5f, the energy should be at %.5f (compare to %.5f)"%(x,
    #calchalfE,
    #acthalfE)

    plot(halflife,calchalfE,'x', color='k', label="Halflife of average (x=%.0f)"%x)

    prop = matplotlib.font_manager.FontProperties(size=10)
    legend(loc=1, prop = prop)
    gca().set_xlim((0,xmax))
    savefig("%s/hartke.%s"%(path,ext), dpi=dpi, bbox_inches="tight")
    cla()

    # First done info
    # Estimate expected finish value for runs that did not complete:
    tol = 1e-10
    guess = 1
    diff = 1e-5
    cutoff = 0.01
    def estimatedFinish(x):
        return bestFitFunction(reg,x,0,True) - (Emin + cutoff)
    x = guess
    val = estimatedFinish(x)
    #print "EstFinish: %.5f, value: %.6f, dx: %.6f"%(x,val,dx)
    while (abs(val) > tol):
        dx = (estimatedFinish(x+diff) - val)/diff
        x = x - val/dx
        # Shouldn't be negative...
        x = abs(x);
        # Shouldn't be infinity
        if x == inf: 
            x = 1e10
            print "Cannot converge, setting expected finished to 1e10 (from inf)."
            break
        val = estimatedFinish(x)
        # Eliminate bad values
        if isnan(val):
            val = 1
        #print "EstFinish: %.5f, value: %.6f, dx: %.6f"%(x,val,dx)

    estFinish = x

    # Replace unfinished runs
    for i in range(len(firstDone)):
        if firstDone[i] == 0: firstDone[i] = estFinish

    fdval = array(firstDone).mean()
    fdstdev = array(firstDone).std()

    #
    # Write summary file:
    #
    str = ""
    str += "nRuns: 		%d"%nRuns + "\n"
    str += "xmax:  		%d"%xmax + "\n"
    str += "fit-a: 		%s"%numErrorText(reg[0], err[0], 5) + "\n"
    str += "fit-b: 		%s"%numErrorText(reg[1], err[1], 5) + "\n"
    str += "rsq: 		%f"%rsq + "\n"
    str += "halflife:		%0.1f"%halflife + "\n"
    str += "percentDone:	%0.1f"%percentDone + "\n"
    str += "calchalfE: 		%f"%calchalfE + "\n"
    str += "acthalfE:		%f"%acthalfE + "\n"
    str += "ave_energy:		%f"%ave_energy + "\n"
    str += "Emin:		%f"%Emin + "\n"
    str += "estFinish:		%d"%estFinish + "\n"
    str += "first-val:		%d"%fdval + "\n"
    str += "first-std:		%d"%fdstdev + "\n"
    str += "killed:		%0.1f (%d)"%(killed/float(duplicate+optimized+killed)*100,killed) + "\n"
    str += "duplicate:		%0.1f (%d)"%(duplicate/float(duplicate+optimized+killed)*100,duplicate) + "\n"
    str += "optimized:		%0.1f (%d)"%(optimized/float(duplicate+optimized+killed)*100,optimized) + "\n"

    f = open("%s/summary"%path, 'w')
    f.write(str)
    f.close()
    print "Summary file written."
    #print str

def generateResults(tests):
    # Extract data from test files
    columns = {
        "numInitial"			: "$N_i$",
        "popSize"			: "$N_\\text{pool}$",
        "genTotal"			: "$N_\\text{cont}$",
        "p_her"				: "$p_\\text{her}$",
        "p_mut"				: "$p_\\text{mut}$",
        "p_perm"			: "$p_\\text{perm}$",
        "her_minimumContribution"	: "$p_\\text{her,min}$",
        "perm_ex"			: "$N_\\text{perm,ex}$",
        "mut_strainStdev"		: "$\\sigma_\\text{mut,max}$",
        "mut_amp"			: "$\\rho_\\text{max}$",
        "mut_per1"			: "$\\eta$",
        "mut_per2"			: "$\\mu$",
        "perm_strainStdev"		: "$\\sigma_\\text{perm,max}$",
        "tol_enthalpy"			: "$\\Delta H$",
        "tol_volume"			: "$\\Delta V$",
        "nRuns"				: "nRuns",
        "xmax"				: "Structures per run",
        "rsq"				: "$R^2$",
        "halflife"			: "$i_\\frac{1}{2}$",
        "percentDone"			: "\\%c",
        "calchalfE"			: "calchalfE",
        "acthalfE"			: "acthalfE",
        "ave_energy"			: "ave_energy",
        "Emin"				: "Emin",
        "first-val"			: "first-val",
        "first-std"			: "first-std",
        "killed"			: "\\%k",
        "duplicate"			: "\\%d",
        "optimized"			: "\\%o"
        }

    data = {}
    for test in tests:
        testdata = {}
        for file in ["%s/summary"%test, "%s/xtalopt.state"%test]:
            f = open(file)
            for line in f.readlines():
                for key in columns.keys():
                    if line.startswith(key):
                        lsplit = line.split()
                        value = "N/a"
                        if len(lsplit) > 1:
                            value = lsplit[1]
                        testdata[key] = value
        data[test] = testdata

    # Clean up data -- e.g. don't report perm_ex if p_perm == 0
    for d in data.keys():
        if 'p_her' in data[d] and int(data[d]["p_her"]) == 0:
            data[d]['her_minimumContribution'] = "X"
        if 'p_mut' in data[d] and int(data[d]["p_mut"]) == 0:
            data[d]['mut_amp'] = "X"
            data[d]['mut_strainStdev'] = "X"
            data[d]['mut_per1'] = "X"
            data[d]['mut_per2'] = "X"
        if 'p_perm' in data[d] and int(data[d]["p_perm"]) == 0:
            data[d]['perm_strainStdev'] = "X"
            data[d]['perm_ex'] = "X"

    usecolumns = [
        "nRuns",
        "popSize",
        "p_her",
        "her_minimumContribution",
        "p_mut",
        "mut_strainStdev",
        "mut_amp",
        "mut_per1",
        "mut_per2",
        "p_perm",
        "perm_strainStdev",
        "perm_ex",
        "tol_enthalpy",
        "tol_volume",
        "halflife",
        "percentDone",
        "duplicate"
        ]

    # Sort searches
    searches = data.keys()
    sortkey = 'halflife'
    for i in range(len(searches)):
        si = searches[i]
        for j in range(i,len(searches)):
            sj = searches[j]
            if sortkey not in data[sj] or data[si][sortkey] > data[sj][sortkey]:
                si, searches[i], searches[j] = sj, sj, si

    # Define format for table
    s  = ""
    s += "\\begin{longtable}{|c|" # One for index
    for i in usecolumns:
        s += "c|" # and one for each column
    s += "}\\hline\n"
    # Headers
    s += "Search \\#"
    for col in usecolumns:
        s += " & "
        s += columns[col]
    s += " \\endhead\\hline\\hline\n"
    # Data
    i = 1
    for search in searches:
        s += "%s"%i
        for col in usecolumns:
            s += " & "
            val = "M"
            if search in data and col in data[search]:
                val = data[search][col]
            s += val
        s += "\\\\\\hline\n"
        i += 1
    s+= "\\end{longtable}\n"
    f = open("../results/tables/table.tex", 'w')
    f.write(s)

    # Use same search ordering to render the average plot
    fig = figure(figsize=(6.5,9))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    i = 1;
    maxval = 600;
    vals = []
    halflives = []
    for search in searches:
        val = 0
        std = 0
        halflife = 0
        if 'first-val' in data[search]:
            val = float(data[search]['first-val'])
        if 'first-std' in data[search]:
            std = float(data[search]['first-std'])
        if 'halflife' in data[search]:
            halflife = float(data[search]['halflife'])
        ax1.errorbar(val, i, xerr=std, ecolor='k')
        halflives.append(halflife)
        vals.append(val)
#        if maxval < val+std: maxval = val+std
        i += 1

    maxval += 0.05*maxval
    for s in range(1,i):
        ax1.plot((0,maxval), (s,s), 'k:')
        ax2.plot((0,maxval), (s,s), 'k:')

    ax1.plot(vals, range(1,i), 'b-')
    ax2.plot(halflives, range(1,i), 'r-')

    ax1.set_xlim(0,maxval)
    ax1.set_ylim(0,i)
    ax1.xaxis.set_ticks(range(0,int(maxval),50))
    ax1.yaxis.set_ticks(range(1,i))

    ax2.set_xlim(0.95*min(halflives), 1.05*max(halflives))
    ax2.set_ylim(0,i)
    ax2.xaxis.set_ticks(range(int(0.95*min(halflives)-1), int(1.05*max(halflives)+1),5))
    ax2.yaxis.set_ticks(range(1,i))

    fontsize = 6
    for ax in [ax1, ax2]:
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)

    ax1.set_title("Average first appearance of rutile", fontsize=fontsize+4)
    ax2.set_title("Halflife", fontsize=fontsize+4)
    ax1.set_xlabel("Number of structures", fontsize=fontsize+1)
    ax2.set_xlabel("Number of structures", fontsize=fontsize+1)
    ax1.set_ylabel("Search #", fontsize=fontsize+1)
    ax2.set_ylabel("Search #", fontsize=fontsize+1)

    savefig("../results/images/dataplot.pdf", bbox_inches='tight')
