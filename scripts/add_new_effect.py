#!/usr/bin/python
# Call this after adding a new effect source file, e.g., effect.c.  It takes all the .c files in reboundx/src and adds them to all the Makefiles where needed

import glob
import os
import sys

sources = [each for each in os.listdir('../src/') if each.endswith('.c')] 
print("Adding {0}".format(sources))
sourcestring = "SOURCES="
sourcesExamples = "SOURCES="
sourcesSetup = "        sources = ["
for source in sources:
    sourcestring += source + " "
    sourcesExamples += "$(REBOUNDXDIR)" + source + " "
    sourcesSetup += ' \'src/' + source + '\','
    
sourcestring += "\n"
sourcesExamples += "\n"
sourcesSetup = sourcesSetup[:-1] + '],\n' # remove comma after last source

with open("../src/Makefile") as f:
    lines = f.readlines()
    for i,l in enumerate(lines):
        if "modify_orbits_direct.c" in l:
            lines[i] = sourcestring

    with open("../src/Makefile", "w") as f:
        f.writelines(lines)

with open("../examples/Makefile") as f:
    lines = f.readlines()
    for i,l in enumerate(lines):
        if "modify_orbits_direct.c" in l:
            lines[i] = sourcesExamples

    with open("../examples/Makefile", "w") as f:
        f.writelines(lines)
'''
for example in os.listdir('../examples/'):
    try:
        with open("../examples/" + example + "/Makefile") as f:
            lines = f.readlines()
            for i,l in enumerate(lines):
                if "modify_orbits_direct.c" in l:
                    print('here')
                    lines[i] = sourcesExamples

            with open("../examples/" + example + "/Makefile", "w") as f:
                f.writelines(lines)
    except:
        pass
'''
with open("../setup.py") as f:
    setuplines = f.readlines()
    for i,l in enumerate(setuplines):
        if 'src/' in l:
            setuplines[i] = sourcesSetup

    with open("../setup.py", "w") as f:
        f.writelines(setuplines)

