#!/usr/bin/env python

import sys

filename = sys.argv[1]

class convTable:
    def __init__(self):
        self.m_lookup = {}
        self.m_bools = []
        self.addNew("optType", "edit\\optType")
        self.addNew("numInitial", "opt\\opt\\numInitial")
        self.addNew("popSize", "opt\\opt\\popSize")
        self.addNew("genTotal", "opt\\opt\\contStructs")
        self.addNew("p_her", "opt\\opt\\p_cross")
        self.addNew("p_mut", "opt\\opt\\p_strip")
        self.addNew("p_perm", "opt\\opt\\p_perm")
        self.addNew("her_minimumContribution", "opt\\opt\\cross_minimumContribution")
        self.addNew("perm_ex", "opt\\opt\\perm_ex")
        self.addNew("perm_strainStdev_max", "opt\\opt\\perm_strainStdev_max")
        self.addNew("mut_strainStdev_min", "opt\\opt\\strip_strainStdev_min")
        self.addNew("mut_strainStdev_max", "opt\\opt\\strip_strainStdev_max")
        self.addNew("mut_amp_min", "opt\\opt\\strip_amp_min")
        self.addNew("mut_amp_max", "opt\\opt\\strip_amp_max")
        self.addNew("mut_per1", "opt\\opt\\strip_per1")
        self.addNew("mut_per2", "opt\\opt\\strip_per2")
        self.addNew("A_min", "init\\limits\\a\\min")
        self.addNew("B_min", "init\\limits\\b\\min")
        self.addNew("C_min", "init\\limits\\c\\min")
        self.addNew("alpha_min", "init\\limits\\alpha\\min")
        self.addNew("beta_min", "init\\limits\\beta\\min")
        self.addNew("gamma_min", "init\\limits\\gamma\\min")
        self.addNew("A_max", "init\\limits\\a\\max")
        self.addNew("B_max", "init\\limits\\b\\max")
        self.addNew("C_max", "init\\limits\\c\\max")
        self.addNew("alpha_max", "init\\limits\\alpha\\max")
        self.addNew("beta_max", "init\\limits\\beta\\max")
        self.addNew("gamma_max", "init\\limits\\gamma\\max")
        self.addNew("vol_min", "init\\limits\\volume\\min")
        self.addNew("vol_max", "init\\limits\\volume\\max")
        self.addNew("vol_fixed", "init\\limits\\volume\\fixed")
        self.addNew("shortestInteratomicDistance", "init\\limits\\shortestInteratomicDistance")
        self.addNew("tol_enthalpy", "opt\\tol\\enthalpy")
        self.addNew("tol_volume", "opt\\tol\\volume")
        self.addNew("using_fixed_volume", "init\\using\\fixedVolume")
        self.addNew("using_shortestInteratomicDistance", "init\\using\\shortestInteratomicDistance")
        self.addNew("using_remote", "sys\\remote")
        self.addBool("using_fixed_volume")
        self.addBool("using_shortestInteratomicDistance")
        self.addBool("using_remote")
        self.addNew("limitRunningJobs", "opt\\opt\\limitRunningJobs")
        self.addNew("runningJobLimit", "opt\\opt\\runningJobLimit")
        self.addNew("failLimit", "opt\\opt\\failLimit")
        self.addNew("failAction", "opt\\opt\\failAction")
        self.addNew("filePath", "sys\\file\\path")
        self.addNew("fileBase", "sys\\file\\base")
        self.addNew("launchCommand", "sys\\queue\\qsub")
        self.addNew("queueCheck", "sys\\queue\\qstat")
        self.addNew("queueDelete", "sys\\queue\\qdel")
        self.addNew("host", "sys\\remote\\host")
        self.addNew("username", "sys\\remote\\username")
        self.addNew("rempath", "sys\\remote\\rempath")

    def addNew(self,oldkey,newkey):
        self.m_lookup[oldkey] = newkey

    def addBool(self,oldkey):
        self.m_bools.append(oldkey)

    def convert(self,oldkey,value):
        if oldkey not in self.m_lookup.keys():
            print "Unknown key: %s"%oldkey
            return "%s=%s\n"%(oldkey,value.strip())
        if oldkey in self.m_bools:
            if int(value) == 1:	value = "true"
            else: 		value = "false"
        return "%s=%s\n"%(self.m_lookup[oldkey],value.strip())

    def has(self, key):
        return (key in self.m_lookup.keys())

t = convTable()

old = file (filename, 'r')
new = file (filename + ".converted", 'w')
opts =file (filename + ".extra", 'w')
new.write("[xtalopt]\n")
for line in old:
    # Read and convert keys
    key,junk,val = line.partition(":")
    if t.has(key):
        new.write(t.convert(key,val))
    # Copy templates to new file
    else:
        opts.write(line)
    # Set composition
    if key == "Composition":
        comp = val.split()
        size = 0
        for i in range(len(comp)/2):
            size +=1
            new.write("init\composition\%d\\atomicNumber=%d\n"%(i+1,int(comp[2*i])))
            new.write("init\composition\%d\quantity=%d\n"%(i+1,int(comp[2*i+1])))
        new.write("init\composition\size=%d\n"%(size))

print "%s converted!"%filename
