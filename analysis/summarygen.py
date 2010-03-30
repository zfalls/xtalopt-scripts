#!/usr/bin/env python
from dataAnalysis import *
from summarytool import *
import sys,os

curdir = os.environ['PWD']
tests = os.listdir(curdir)

force = False
if len(sys.argv) > 1 and sys.argv[1] == "force":
    force = True

# find and summarize results files
for test in tests:
    files = os.listdir(test)
    # remove all non results files from list
    i = 0
    badidx = []
    for file in files:
        if not file.startswith("run") or not file.endswith("results.txt"):
            badidx.append(i)
        i+=1
    badidx.reverse()
    for idx in badidx:
        files.pop(idx)

    print test
    generateSummary(test, files, force)

# make table and plot
#generateResults(tests)
