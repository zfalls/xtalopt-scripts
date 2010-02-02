import csv
import time
from pylab import *
from numpy import *
from numpy import zeros as ZEROS
from numpy.linalg import norm
from scipy import stats,linspace,polyfit
from scipy.optimize import *
from scipy.interpolate import splprep,splev
from matplotlib import rc
from matplotlib.font_manager import *

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

def getFit(xdata,ydata,func,guess):
	"""
	x- and ydata are the values to be fit to a function func, which must be declared as:

	def func(guess,xdata,ydata)
		...
		error = norm(y - ydata)
		return error

	The return value must be the norm of the error vector to have R^2 calculated.
	"""
	lsq 		= leastsq(func,guess,args=(xdata,ydata),full_output=1)
	fit 		= lsq[0]
	error 		= []
	resid 		= lsq[2]['fvec']
	chisquare	= sum(resid**2)
	try:
		dof	= len(ydata) - len(fit)
	except:
		dof	= len(ydata) - 1
	if lsq[1] == None: 
		try:	error = ZEROS(len(fit))
		except TypeError: error = ZEROS(1)
	else:
		cov 	= lsq[1]
		for i in range(len(cov[0])):
			error.append(sqrt(cov[i][i])*sqrt(chisquare/dof))
	error = array(error)
	RS = norm(resid) ** 2
	SS = ((ydata - ydata.mean()) ** 2).sum()
	rsq = 1 - (RS/SS)
	return fit,error,rsq
	
def getCustomFit(xdata,ydata,func,guess,extFrac = 0,res = 200):
	"""
	func must be a function that has args:

	def func(c,x,y,get_y=False):
		y_model = [function of x]
		if get_y == True: return y_model
		error = [error determination]
		return error
		
	guess is an array of guesses for the fit parameters.
	"""
	fit,error,rsq = getFit(xdata,ydata,func,guess)
	xmin,xmax = xdata.min(),xdata.max()
	ext = (xmax - xmin) * extFrac
	x = linspace(xmin-ext,xmax+ext,res)
	y = func(fit,x,None,get_y=True)
	return x,y,fit,rsq,error

def getLinearFit(xdata,ydata,extFrac = 0,res = 200,guess=None,forceOrigin=False):
	if forceOrigin:
		if guess == None: guess = ones(1)
		def func(c,x,y,get_y=False):
			y_model = c * x
			if get_y == True: return y_model
			error = y-y_model
			return error
	else:
		if guess == None: guess = ones(2)
		def func(c,x,y,get_y=False):
			y_model = c[0] * x + c[1]
			if get_y == True: return y_model
			error = y-y_model
			return error
	fit,error,rsq = getFit(xdata,ydata,func,guess)
	xmin,xmax = xdata.min(),xdata.max()
	ext = (xmax - xmin) * extFrac
	x = linspace(xmin-ext,xmax+ext,res)
	y = func(fit,x,None,get_y=True)
	return x,y,fit,rsq,error

def getQuadraticFit(xdata,ydata,extFrac = 0,res = 200,guess=None,forceOrigin=False):
	if forceOrigin:
		if guess == None: guess = ones(2)
		def func(c,x,y,get_y=False):
			y_model = c[0] * x**2 + c[1] * x
			if get_y == True: return y_model
			error = y-y_model
			return error
	else:
		if guess == None: guess = ones(3)
		def func(c,x,y,get_y=False):
			y_model = c[0] * x**2 + c[1] * x + c[2]
			if get_y == True: return y_model
			error = y-y_model
			return error
	fit,error,rsq = getFit(xdata,ydata,func,guess)
	xmin,xmax = xdata.min(),xdata.max()
	ext = (xmax - xmin) * extFrac
	x = linspace(xmin-ext,xmax+ext,res)
	y = func(fit,x,None,get_y=True)
	return x,y,fit,rsq,error

def getCubicFit(xdata,ydata,extFrac = 0,res = 200,guess=None,forceOrigin=False):
	if forceOrigin:
		if guess == None: guess = ones(3)
		def func(c,x,y,get_y=False):
			y_model = c[0] * x**3 + c[1] * x**2 + c[2] * x
			if get_y == True: return y_model
			error = y-y_model
			return error
	else:
		if guess == None: guess = ones(4)
		def func(c,x,y,get_y=False):
			y_model = c[0] * x**3 + c[1] * x**2 + c[2] * x + c[3]
			if get_y == True: return y_model
			error = y-y_model
			return error
	fit,error,rsq = getFit(xdata,ydata,func,guess)
	xmin,xmax = xdata.min(),xdata.max()
	ext = (xmax - xmin) * extFrac
	x = linspace(xmin-ext,xmax+ext,res)
	y = func(fit,x,None,get_y=True)
	return x,y,fit,rsq,error

def getGaussianFit(xdata,ydata,extFrac = 0, res = 200,guess=None):
	if guess == None: 
		guess = array([ydata.max(),ydata.argmax(),(xdata.max() - xdata.min()) / 3 ])
	def func(c,x,y,get_y=False):
		y_model = c[0]*exp(-(x-c[1])**2/(c[2]**2))
		if get_y == True: return y_model
		error = y-y_model
		return error
	fit,error,rsq = getFit(xdata,ydata,func,guess)
	xmin,xmax = xdata.min(),xdata.max()
	ext = (xmax - xmin) * extFrac
	x = linspace(xmin-ext,xmax+ext,res)
	y = func(fit,x,None,get_y=True)
	return x,y,fit,rsq,error

def getExponentialFit(xdata,ydata,extFrac = 0, res = 200,guess=None):
	if guess == None: guess = array([1,1,0,0])
	def func(c,x,y,get_y=False):
		y_model = c[0]*(e**(-c[1]*x+c[2])) + c[3]
		if get_y == True: return y_model
		error = y-y_model
		return error
	# Remove "unsafe values":  (i.e., y = 0). Don't know if it helps -- try out if you get problems
	"""
	bad = []
	for ind,val in enumerate(ydata):
		if val == 0: bad.append(ind)
	bad.reverse()
	for ind in bad:
		print "Warning! Excluding data point %5s due to incompatibility with fit-type. Point is: (%.5g,%.5g)" % (ind,xdata[ind],ydata[ind])
		xdata.pop(ind)
		ydata.pop(ind)
	"""
	fit,error,rsq = getFit(xdata,ydata,func,guess)
	xmin,xmax = xdata.min(),xdata.max()
	ext = (xmax - xmin) * extFrac
	x = linspace(xmin-ext,xmax+ext,res)
	y = func(fit,x,None,get_y=True)
	return x,y,fit,rsq,error

def getLogFit(xdata,ydata,extFrac = 0, res = 200,guess=None):
	if guess == None: guess = ones(2)
	def func(c,x,y,get_y=False):
		y_model = c[0]*log(x) + c[1]  # NOTE: the log(x) function in python is the inverse of e^(x), not 10^(x).
		if get_y == True: return y_model
		error = y-y_model
		return error
	# Remove "unsafe values":  (i.e., x = 0). Not sure if it helps, try out if you have trouble.
	"""
	bad = []
	for ind,val in enumerate(xdata):
		if val == 0: bad.append(ind)
	bad.reverse()
	for ind in bad:
		print "Warning! Excluding data point %5s due to incompatibility with fit-type. Point is: (%.5g,%.5g)" % (ind,xdata[ind],ydata[ind])
		array(xdata.tolist().pop(ind))
		array(ydata.tolist().pop(ind))
	"""
	fit,error,rsq = getFit(xdata,ydata,func,guess)
	xmin,xmax = xdata.min(),xdata.max()
	ext = (xmax - xmin) * extFrac
	x = linspace(xmin-ext,xmax+ext,res)
	y = func(fit,x,None,get_y=True)
	return x,y,fit,rsq,error

def getHyperbolicFit(xdata,ydata,extFrac = 0, res = 200,guess=None):
	if guess == None: 
		guess = array([1])
	def func(c,x,y,get_y=False):
		y_model = c[0]/x
		if get_y == True: return y_model
		error = y-y_model
		return error
	fit,error,rsq = getFit(xdata,ydata,func,guess)
	xmin,xmax = xdata.min(),xdata.max()
	ext = (xmax - xmin) * extFrac
	x = linspace(xmin-ext,xmax+ext,res)
	y = func([fit],x,None,get_y=True)
	return x,y,fit,rsq,error

def getInterpolation(xdata,ydata,extFrac,res,smoothness=0):
	"Interpolates x,y data"
	xmin,xmax = xdata.min(),xdata.max()
	ext = (xmax - xmin) * extFrac
	x = linspace(xmin-ext,xmax+ext,res)
	s=smoothness # smoothness parameter
	k=2 # spline order
	nest=-1 # estimate of number of knots needed (-1 = maximal)
	# find the knot points
	tckp,u = splprep([xdata,ydata],s=s,k=k,nest=nest)
	# evaluate spline, including interpolated points
	x,y = splev(linspace(0,1,res),tckp)
	return x,y


def plotFitData(xdata, ydata, 
		xlab="x", ylab="y", plotTitle=None, 
		regType=None, 
		forceOrigin=False, 
		dataLabel="Data", fitLabel=None, 
		extFrac=0,
		res=200, 
		plotData=True, 
		guess=None, 
		errordata=None, 
		xerr=None, yerr=None, 
		color='b', errColor='r', 
		dataFormat='o', fitFormat='-', 
		fitAlpha=1, dataAlpha=1, errAlpha=1, 
		linewidth=1,
		smoothness=0, 
		returnFitXY=False, 
		labelsize=10, titlesize=10,
		func=None, 
		customFitLabel=None,
		ax = gca()):
	"""Temp help:

	regTypes: "None", "linear", "exponential", "log", "quadratic", "cubic", "connect", "gaussian", "hyperbolic", "custom", "inter"

	errordata is deprecated -- use yerr instead.

	smoothness is used for spline calculations -- play around with it to get what you want.
	
	If using regType='custom':
		func must be a function that has args:

		def func(c,x,y,get_y=False):
			y_model = [function of x]
			if get_y == True: return y_model
			error = [error determination]
			return error
			
		guess is an array of guesses for the fit parameters.

		customFitLabel should be set if fitLabel="regEqu".
	"""
	ret = 0
	if regType == None:	regType = "None"

	xdata,ydata = sameSize(xdata,ydata,"x","y")
	
	if errordata != None: 
		print "Using 'errordata' is depreciated. Use 'yerr' instead."
		yerr = errordata
	x,y,reg = 0,0,0
	if regType == "custom":				
		if func != None and guess != None:
			x,y,reg,rsq,error	 = getCustomFit(xdata,ydata,func,guess,extFrac,res)
		else:
			print "Psychic determination of fit function and initial guess not yet implemented. Please supply."
			raise Exception, "Need function and/or guess"
	if regType == "linear":				
		x,y,reg,rsq,error	 = getLinearFit(xdata,ydata,extFrac,res,guess,forceOrigin)
	if regType == "exponential":			
		x,y,reg,rsq,error	 = getExponentialFit(xdata,ydata,extFrac,res,guess)
	if regType == "hyperbolic":			
		x,y,reg,rsq,error	 = getHyperbolicFit(xdata,ydata,extFrac,res,guess)
	if regType == "log":				
		x,y,reg,rsq,error	 = getLogFit(xdata,ydata,extFrac,res,guess)
	if regType == "gaussian":			
		x,y,reg,rsq,error	 = getGaussianFit(xdata,ydata,extFrac,res,guess)
	if regType == "quadratic":			
		x,y,reg,rsq,error	 = getQuadraticFit(xdata,ydata,extFrac,res,guess,forceOrigin)
	if regType == "cubic":				
		x,y,reg,rsq,error	 = getCubicFit(xdata,ydata,extFrac,res,guess,forceOrigin)
	if regType == "inter":
		x,y = getInterpolation(xdata,ydata,extFrac,res,smoothness)

	if regType != "None" and regType != "connect" and regType != None and regType != "inter" and fitLabel == "regEqu":
		if len(error) == 1:
			bestVals = numErrorTex(reg,error,5)
		else:
			bestValsArray = array([reg,error]).transpose()
			bestVals = []
			for cur in bestValsArray:
				bestVals.append(numErrorTex(cur[0],cur[1],5))
			bestVals = tuple(bestVals)
	
	if fitLabel == "regEqu":
		if forceOrigin:
			if regType == "linear":			fitLabel = "$y=(%s)x$" % bestVals
			if regType == "quadratic":		fitLabel = "$y=%sx^2+%sx$" % bestVals 
			if regType == "cubic":			fitLabel = "$y=%sx^3+%sx^2+%sx$" % bestVals
		else:
			if regType == "linear":			fitLabel = "$y=(%s)x+%s$" % bestVals
			if regType == "quadratic":		fitLabel = "$y=%sx^2+%sx+%s$" % bestVals 
			if regType == "cubic":			fitLabel = "$y=%sx^3+%sx^2+%sx+%s$" % bestVals
		if regType == "exponential":			fitLabel = "$y=(%s)e^{-(%s)x+%s}+%s$" % bestVals 
		if regType == "hyperbolic":			fitLabel = "$y=\\frac{%s}{x}$" % bestVals 
		if regType == "log":				fitLabel = "$y=%s\\ln(x)+%s$" % bestVals 
		if regType == "gaussian":			fitLabel = "$y=%se^{-\\frac{(x-%s)^2}{%s^2}}$" % bestVals
		if regType == "inter":
			fitLabel = "Interpolated values"
		if customFitLabel != None:			fitLabel = customFitLabel % bestVals 

	fitFullFormat 	= '%s%s' % (color,fitFormat)
	
	if regType != "None" and regType != "connect":	ax.plot(x,y,fitFullFormat,label = fitLabel, alpha=fitAlpha, linewidth=linewidth)
	if regType != "connect" and plotData:		ax.scatter(xdata,ydata,c = color, marker = dataFormat,label = dataLabel, alpha=fitAlpha)
	elif regType == "connect" and plotData != True:	ax.plot(xdata,ydata,"%s%s" % (str(color),str(fitFormat)),label = dataLabel, alpha=dataAlpha, linewidth=linewidth)
	elif plotData:					ax.plot(xdata,ydata,"%s%s" % (str(color),str(fitFormat)),marker=dataFormat,label = dataLabel, alpha=dataAlpha, linewidth=linewidth)

	if xerr != None or yerr != None:		ax.errorbar(xdata,ydata,xerr=xerr,yerr=yerr,color=errColor,linestyle='None',alpha=errAlpha)

	ax.set_xlabel(xlab, size=labelsize)
	ax.set_ylabel(ylab, size=labelsize)
	ax.set_title(plotTitle, size=titlesize)
	if regType != "connect" and regType != "None" and regType != "inter":
		ret = [reg,error,rsq]
	if returnFitXY == True:
		if ret:
			ret = ret[0],ret[1],ret[2],(x,y)
		else:
			ret = (x,y)
	if ret:
		return ret

def roundTex(x,sigfigs):
	format = "%%#.%dg" % sigfigs
	x = format % x
	x = latexSciNot(x)
	return x

def numErrorTex(val,e,maxSigFigs = None):
	"""Given an number and the error, rounds off each appropriately and add \pm in between them.
	
	ex. numErrorTex(45.3223,.04567) returns:
	"45.32 \pm 0.05"
	"""
	val,e 	= numError(val, e, maxSigFigs)
	return "%s \\pm %s" % (val,e)

def numErrorText(val,e,maxSigFigs = None):
	"""Given an number and the error, rounds off each appropriately and adds "+/-" in between them.
	
	ex. numErrorTex(45.3223,.04567) returns:
	"45.32 \pm 0.05"
	"""
	val,e 	= numError(val, e, maxSigFigs)
	return "%s +/- %s" % (val,e)

def numError(val,e,maxSigFigs = None):
	"""Given an number and the error, rounds off each appropriately and returns a tuple:
	
	ex. numError(4532.23,04.567) returns:
	("4.532 \\times 10^{3}", "5.")
	"""
	if e == inf:
		val 	= roundTex(val,5)
		e  	= "\\infty "
	else:
		e	= abs(e)
		errf 	= "%.0e" % abs(e)
		err 	= int(errf[0])
		dec 	= int(errf[2:5])
		valf	= "%.0e" % abs(val)
		vdec	= int(valf[2:5])
		if vdec >= 0 and dec >= 0:	sigfigs	= vdec - dec + 1
		elif vdec <= 0 and dec <= 0:	sigfigs	= abs(vdec - dec) + 1
		else:
			sigfigs	= abs(vdec - dec) + 1
			if sigfigs < 1: sigfigs = 1
			if maxSigFigs != None:	
				if sigfigs > maxSigFigs: 	
					sigfigs = maxSigFigs
					print "sigfigs corrected to %s" % maxSigFigs
		if sigfigs < 0:
			sigfigs = maxSigFigs
			print "sigfigs corrected to %s" % maxSigFigs
		val 	= roundTex(val,sigfigs)
		e 	= roundTex(e,1)
	return (val,e)

def latexSciNot(x):
	"""Pass a string that is formatted using 
	
	'%%.Xg' %% number

	and the 'e+0' or 'e+' will be converted into a LaTeX friendly format."""
	b = x # Default
	try:
		c = x.find("e+0")
		if c != -1: 	x = x[0:c] + "\\times 10^{" + x[c+3:] + "}"
		c = x.find("e+")
		if c != -1:     x = x[0:c] + "\\times 10^{" + x[c+2:] + "}"
		c = x.find("e-0")
		if c != -1:     x = x[0:c] + "\\times 10^{-" + x[c+3:] + "}"
		c = x.find("e-")
		if c != -1:     x = x[0:c] + "\\times 10^{-" + x[c+2:] + "}"
	except:
		print "Something went wrong. Did you check the docstring, dumbass?"
	return x

def sameSize(x,y,xlab="first",ylab="second"):
	x,y = array(x),array(y)
	# check to see if x and y's length can be checked:
	try: 
		if len(x) != len(y):pass
	except: 
		print "Cannot check length of arrays. Probably only 1 value"
		return x,y
	if len(x) > len(y):
		print "Dimension mismatch -- truncating %s array."%xlab
		x = x[0:len(y)]
	if len(y) > len(x):
		print "Dimension mismatch -- truncating %s array."%ylab
		y = y[0:len(x)]
	return x,y

def texTable(data, dataFormat,headers=None,headerFormat="bold",label=None,caption=None):
	"""Creates a latex table from the data provided:
	(Note that data must be a list)

	Example:
	
	size = array(["small",	"medium",	"large"])
	cost = array([2.10, 	3.10, 		5.43])
	volume=array([16,	32,		64])

	headers	= ["Size", "Cost", "Volume"]
	headerFormat = "bold"
	data 	= [size, cost, volume]
	format	= "\\textbf(%s) & %.2f & %d"
	caption	= "Drink sizes"
	label 	= "table:drink"
	
	texTable(headers,data,format,label=label, caption=caption)
	"""
	dataList = []
	for rnum in range(len(data[0])):
		row = []
		for cnum in range(len(data)):
			row.append(data[cnum][rnum])
		dataList.append(tuple(row))

	width = len(dataList[0])
	if width == 0:
		print "No data for the table? Try again!"
		return
	just = ""
	for i in range(width):
		just = just + "|c"
	just = just + "|"
	
	if headerFormat in ["bold"]:
		head = ""
		for num,value in enumerate(headers):
			head = head + "\\textbf{%s} " % value
			if num < len(headers)-1:	head = head + "& "
			else:				head = head + "\\\\\\hline\\hline"
	
	print "\\begin{tabular}{%s}\\hline" % just
	if headers != None:
		print "  " + head
	for i in range(len(dataList)):
		print "  " + dataFormat % dataList[i] + " \\\\\\hline"
	print "\\end{tabular}"
	if caption != None:
		print "\\caption{%s}" % caption
	if label != None:
		print "\\label{%s}" % label

def texTableFile(filename, data, dataFormat,headers=None,headerFormat="bold",label=None,caption=None):
	"""Creates a latex table from the data provided:
	(Note that data must be a list)

	Example:
	
	size = array(["small",	"medium",	"large"])
	cost = array([2.10, 	3.10, 		5.43])
	volume=array([16,	32,		64])

	headers	= ["Size", "Cost", "Volume"]
	headerFormat = "bold"
	data 	= [size, cost, volume]
	format	= "\\textbf(%s) & %.2f & %d"
	caption	= "Drink sizes"
	label 	= "table:drink"
	
	texTable(headers,data,format,label=label, caption=caption)
	"""
        out = file(filename+".tex",'w')
	dataList = []
	for rnum in range(len(data[0])):
		row = []
		for cnum in range(len(data)):
			row.append(data[cnum][rnum])
		dataList.append(tuple(row))

	width = len(dataList[0])
	if width == 0:
		print "No data for the table? Try again!"
		return
	just = ""
	for i in range(width):
		just = just + "|c"
	just = just + "|"
	
	if headerFormat in ["bold"]:
		head = ""
		for num,value in enumerate(headers):
			head = head + "\\textbf{%s} " % value
			if num < len(headers)-1:	head = head + "& "
			else:				head = head + "\\\\\\hline\\hline\n"
	
	out.write( "\\begin{tabular}{%s}\\hline\n" % just)
	if headers != None:
		out.write( "  " + head)
	for i in range(len(dataList)):
		out.write( "  " + dataFormat % dataList[i] + " \\\\\\hline\n")
	out.write( "\\end{tabular}\n")
	if caption != None:
		out.write( "\\caption{%s}\n" % caption )
	if label != None:
		out.write( "\\label{%s}\n" % label)

def getHeaders(file):
    datafile = csv.reader(open(file,"rb"),delimiter="\t",skipinitialspace="true")
    headers = datafile.next()
    return headers
	
def getCsvArray(file,r = 0, delim="\t"):
    datafile = csv.reader(open(file,"rb"),delimiter=delim,skipinitialspace="true")
    data = []
    for row in datafile:
	    if len(row) - 1 < r:
		    print "Warning: Length too small (len=%d, ind=%d)"%(len(row),r)
	    else:
		    try:
			    datum = float(row[r])
		    except:
			    pass
		    else:
			    data.append(datum)
    return array(data)

def getCsvStrings(file,r = 0, quiet=False):
    datafile = csv.reader(open(file,"rb"),delimiter="\t",skipinitialspace="true")
    data = []
    for row in datafile:
	    if len(row) - 1 < r:
		    if not quiet: print "Warning: Length too small (len=%d, ind=%d)"%(len(row),r)
	    else:
		    datum = row[r]
		    data.append(str(datum))
    return array(data)

def writeCsvFile(file,data,headers = None):
	"""Pass data as an array of array, i.e:
	
	headers = ["x","y"]
	x = linspace(0,40,100)
	y = sin(x*pi)

	writeCsvArray(file,array(x,y),headers)
	"""
	datafile = csv.writer(open(file,"wb"),delimiter="\t")
	data = array(data)
	if headers != None: datafile.writerow(headers)
	
	if len(data.shape) == 0:	
					datafile.writerow([data])
					return
	if len(data.shape) == 1:	cols,rows = 1,data.shape[0]
	if len(data.shape) == 2:	cols,rows = data.shape[0],data.shape[1]

	for row in arange(rows):
		rowdata = []
		
		if cols > 1:
			for col in arange(cols):
				rowdata.append(data[col][row])
		else:
			rowdata.append(data[row])

		datafile.writerow(rowdata)

def uStdev(ar):
	"""Unbiased standard deviation"""
	if len(ar) == 1:
		print "Cannot take uStdev of a single measurement. Returning 0."
		return 0
	d = []
	l = len(ar)
	m = ar.mean()
	ar = array(ar)
	for ea in ar:
		d.append(abs(ea - m))
	d = array(d) ** 2
	ustdev = sqrt(d.sum()/(l-1))
	return ustdev

def uStdevOfMean(ar):
	"""Unbiased standard deviation of the mean"""
	ustdev = uStdev(ar)
	return ustdev/sqrt(len(ar))

def interpolateQuad(xdata,ydata,values):
	"Given x and y data along with a list of y values, a regression is performed and an array of corresponding x values is returned. Rather unnecessary, indeed."
	newData = []
	(a,b,c) = polyfit(xdata,ydata,2)
	for x in values:
		newData.append( ( - b + sqrt(b ** 2 - (4 * a * ( c - x ) ) ) ) / (2 * a)   )
	return array(newData)

def averageXY(x,y,xWidth=0):
	"""Condenses large data sets from:
	
	x,y ---> x,y,xUStdev,yUSstdev
	
	by combining like values of x. i.e.,
	x 	= array([0,0,1,1,2,2,3,3])
	y 	= array([5,5,3,3,2,2,1,1])

	into

	x 	= array([0,1,2,3])
	y 	= array([5,3,2,1])
	ustdev 	= array([0,0,0,0])

	xWidth allows some variation in the x data.
	"""
	if len(x) != len(y):
		print "Dimensions of x and y must match!"
		return
	
	x = array(x).tolist()
	y = array(y).tolist()
	
	newx = []
	newy = []
	newxUStdevs = []
	newyUStdevs = []

	while len(x) > 0:
		curx = [x.pop(0)]
		cury = [y.pop(0)]
		matches = []
		for ind in range(len(x)):
			if x[ind] >= mean(curx)-xWidth and x[ind] <= mean(curx)+xWidth:
				matches.append(ind)

		matches.reverse() # reverses values so they can be popped correctly
		for ind in matches:
			curx.append(x.pop(ind))
			cury.append(y.pop(ind))

		newx.append(array(curx).mean())
		newy.append(array(cury).mean())
		newxUStdevs.append(uStdev(array(curx)))
		newyUStdevs.append(uStdev(array(cury)))
	
	return array(newx),array(newy),array(newxUStdevs),array(newyUStdevs)

def findNextInflectionPoint(ar):
	"""
	Performs a second derivative approximation to find the next inflection point in an array. Since an inflection point will occur between two discrete data points, it should be noted that this method returns the index of the point following the change in inflection.
	"""
	ar = array(ar)
	if len(ar) <= 3:
		raise Exception, "Array too short to find inflection points."
	cur		= 0
        prevSlope 	= ar[cur+1] - ar[cur]
       	cur 		+= 1
        slope 		= ar[cur+1] - ar[cur]
        prev_d2y_dx2 	= slope - prevSlope 
        prevSlope 	= slope
	infPoint	= None
	while infPoint == None:
		cur 		+= 1
		if cur >= len(ar)-1:
			raise Exception, "No inflection point in array"
		slope		= ar[cur+1] - ar[cur]
		d2y_dx2		= slope - prevSlope
		prevSign	= prev_d2y_dx2 / abs(prev_d2y_dx2)
		sign 		= d2y_dx2 / abs(d2y_dx2)
		if prev_d2y_dx2 != 0 and d2y_dx2 != 0 and sign != prevSign:
			infPoint = cur
	return infPoint
	
def findInflectionPoints(ar,minData = 3):
	done 	= False
	points 	= []
	prev	= 0
	while done == False:
		try:
			point 	= findNextInflectionPoint(ar[prev:]) + prev
			points.append(point)
			prev 	= point
		except:
			done	= True

	# Enforces minData restriction:
	x 	= 1
	comp 	= points[0]
	while x < len(points):
		cur = points[x]
		if cur - comp < minData:
			a = points.pop(x) # Pop two inflection points to get back to maintain concavity.
			b = points.pop(x)
		else:
			comp = cur
		x += 1
	return array(points)

def decomposeGaussian(xdata,ydata,minData=3,thresholdHeight=0,thresholdRsq=0.9,infPoints=None,plotCurves=False):
	xdata = array(xdata)
	ydata = array(ydata)
	pos	= []
	height	= []
	if infPoints == None:
		infPoints = findInflectionPoints(ydata,minData=minData)
	for peak in arange(len(infPoints)-1):
		ind1			= infPoints[peak]-1
		ind2			= infPoints[peak+1]
		if ydata[ind1:ind2].max() >= thresholdHeight:
			x,y,fit,rsq,error	= getGaussianFit(xdata[ind1:ind2],ydata[ind1:ind2],guess=array([ydata.max(),ydata.argmax(),(xdata.max() - xdata.min()) * 1 ]),extFrac=3)
			if rsq >= thresholdRsq and fit[0] >= thresholdHeight and fit[1] >= xdata.min() and fit[1] <= xdata.max():
				pos.append(fit[1])
				height.append(fit[0])
				if plotCurves == True:	plot(x,y)
	
	return array(pos), array(height)
		
def findNextMinimum(a, start, end):
                """
                Finds the first minima after start and before end in a. Returns the index of said minimum.

                returns minIndex
                """
                cur = start
                prevSlope = a[cur+1] - a[cur]
                minIndex = None
                while minIndex == None:
                        cur += 1
                        if cur == end-1:
                                raise Exception, "No minimum encountered."
                        slope = a[cur+1] - a[cur]
                        if slope >= 0.0 and prevSlope < 0.0: # If slope is positive and previous slope is negative:
                                minIndex = cur
                        prevSlope = slope

                return minIndex

def findNextMaximum(a, start, end):
                """
                Finds the closest maximum after start and before end in a. Returns the index of said maximum.

                returns maxIndex
                """
                cur = start
                prevSlope = a[cur+1] - a[cur]
                maxIndex = None
                while maxIndex == None:
                        cur += 1
                        if cur == end-1:
                                raise Exception, "No maximum encountered."
                        slope = a[cur+1] - a[cur]
                        if slope < 0.0 and prevSlope >= 0.0: # If slope is negative and previous slope is positive:
                                maxIndex = cur
                        prevSlope = slope

                return maxIndex

def printDict(values,subdir=".",ext=".tex"):
	"""Pass a dictionary of ("filename":value) pairs and they will be saved in the subdirectory specified (default is ./) with the specified extentions (default ext=".text") """
	for key,item in values.items():
	        print "%-20s %s" % (key,item)
	        writeCsvFile("%s/%s.tex"%(subdir,key),item)

def valueToIndex(array, value):
	for ind,val in enumerate(array):
		if val >= value: return ind
	return 0

def funcTime(func,len,*args):
	"""Times a function for comparison."""
	start = time.time()
	for i in range(len):
		func(*args)
	end = time.time()
	return end - start

def respace(x,y,res=1):
	"Takes a set of x,y data and changes the resolution of the x set to res"
	xmin	= x.min()
	xmax	= x.max()
	newX	= []
	newY	= []
	curX	= xmin
	curY	= 0
	count	= 0
	for ind in range(len(x)):
		if x[ind] >= curX and x[ind] < curX+res:
			curY += y[ind]
			count+= 1
		elif x[ind] >= curX+res:
			newX.append(curX)
			newY.append(curY/count)
			curX += res
			curY = 0
			count= 0
	return array(newX),array(newY)
