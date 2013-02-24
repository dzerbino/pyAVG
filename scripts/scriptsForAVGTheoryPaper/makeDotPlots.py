import time
import os
import sys
import subprocess

"""Generates all the the graph viz plots
"""

def system(command):
    sts = subprocess.call(command, shell=True, bufsize=-1, stdout=sys.stdout, stderr=sys.stderr)
    if sts != 0:
        raise RuntimeError("Command: %s exited with non-zero status %i" % (command, sts))
    return sts

def main():
    experimentNumber=5
    outputDir = "results"
    print "Making graphviz plots"

    #Write the dot files of the H, G and G_min, G_max and G_med for each experiment
    for experiment in range(experimentNumber):
        dirName = os.path.join(outputDir, "%s_graphViz" % experiment)
        def writeDot(fileName):
            fileName = os.path.join(dirName, fileName)
            system("dot %s -Tpdf > %s.pdf" % (fileName, fileName))
        writeDot("history")
        writeDot("g")
        writeDot("gPRRMin")
        writeDot("gPRRMax")
        writeDot("gPRSMin")
        writeDot("gPRSMax")
     
if __name__ == '__main__':
    main()