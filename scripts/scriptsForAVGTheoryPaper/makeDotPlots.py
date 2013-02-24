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
    print "Making graphviz plots"

    #Write the dot files of the H, G and G_min, G_max and G_med for each experiment
    for experiment in range(experimentNumber):
        historyRow, gRow, gPRows = getRows() 
        dirName = os.path.join(outputDir, "%s_graphViz" % experiment)
        def writeDotPDF(fileName, row):
            fileName = os.path.join(dirName, fileName)
            system("dot %s -Tpdf > %s.pdf" % (fileName, fileName))
        writeDot("history", historyRow)
        writeDot("g", gRow)
        writeDot("gPRRMin" , fn(min, "ubrc"))
        writeDot("gPRRMax", fn(max, "ubrc"))
        writeDot("gPRSMin" , fn(min, "ubsc"))
        writeDot("gPRSMax", fn(max, "ubsc"))
        
    print "Plots took %s seconds" % (time.time() - startTime)

if __name__ == '__main__':
    main()