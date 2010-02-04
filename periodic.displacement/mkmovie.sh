#!/bin/bash
mencoder "mf://*.png" -mf fps=10 -o test.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=800
mencoder test.avi -mf fps=24 -o animation.wmv -ovc lavc -lavcopts vcodec=wmv2:vbitrate=800
