#!/usr/bin/python
# This script automatically creates a list of examples by reading the header in all problem.c files.
import glob

with open("version.txt") as f:
    reboundxversion = f.readlines()[0].strip()
    print "Updating version to "+reboundxversion

with open("src/reboundx.c") as f:
    reboundxlines = f.readlines()
    for i,l in enumerate(reboundxlines):
        if "**VERSIONLINE**" in l:
            reboundxlines[i] = "const char* rebx_version_str = \""+reboundxversion+"\";			// **VERSIONLINE** This line gets updated automatically. Do not edit manually.\n"

    with open("src/reboundx.c", "w") as f:
        f.writelines(reboundxlines)

with open("setup.py") as f:
    setuplines = f.readlines()
    for i,l in enumerate(setuplines):
        if "version='" in l:
            setuplines[i] = "    version='"+reboundxversion+"',\n"

    with open("setup.py", "w") as f:
        f.writelines(setuplines)

