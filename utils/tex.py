"""Lib for generating latex tables.
"""

def formatFloat(string, decimals=3):
    f = float(string)
    if f == 2147483647:
        return "NaN"
    return (("%." + str(decimals) + "f") % f)[1:]

def writeDocumentPreliminaries(fileHandle):
    fileHandle.write("\\documentclass[10pt]{article}\n")
    
    fileHandle.write("\\pagenumbering{arabic}\n")
    fileHandle.write("\\pagestyle{plain}\n")
    
    fileHandle.write("\\usepackage{epsfig}\n")
    fileHandle.write("\\usepackage{url}\n")
    fileHandle.write("\\usepackage{rotating}\n")
    
    fileHandle.write("\\usepackage{multirow}\n")
    
    fileHandle.write("\\usepackage{color}\n")
    fileHandle.write("\\usepackage[table]{xcolor}\n")
    
    fileHandle.write("\\setlength{\\evensidemargin}{0in}\n")
    fileHandle.write("\\setlength{\\oddsidemargin}{0in}\n")
    fileHandle.write("\\setlength{\\marginparwidth}{1in}\n")
    fileHandle.write("\\setlength{\\textwidth}{6.5in}\n")
    fileHandle.write("\\setlength{\\topmargin}{-0.5in}\n")
    fileHandle.write("\\setlength{\\textheight}{9in}\n")
    fileHandle.write("\\usepackage[table]{xcolor}\n")
    fileHandle.write("\\definecolor{lightgray}{gray}{0.9}\n")

    fileHandle.write("\\begin{document}\n")
    
 
def writeDocumentEnd(fileHandle):
    fileHandle.write("\\end{document}\n")
 
def writePreliminaries(columnNumber, fileHandle):
    fileHandle.write("\\begin{table}\n\\centering\n")
    fileHandle.write("\\rowcolors{1}{}{lightgray}\n")
    fileHandle.write("\\begin{tabular}{" + ("c"*columnNumber) + "}\n")

def writeEnd(fileHandle, tableLabel, caption):
    fileHandle.write("\\end{tabular}\n")
    fileHandle.write("\caption{%s}\n" % caption)
    fileHandle.write("\label{%s}\n" % tableLabel)
    fileHandle.write("\end{table}\n")

def writeRow(entries, fileHandle):
    fileHandle.write("%s \\\\\n" % " & ".join(entries))
