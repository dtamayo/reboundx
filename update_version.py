#!/usr/bin/python
# This script automatically creates a list of examples by reading the header in all problem.c files.
import glob
import subprocess
ghash = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("ascii").strip()

with open("version.txt") as f:
    reboundxversion = f.readlines()[0].strip()
    print("Updating version to "+reboundxversion)

with open("doc/index.rst") as f:
    index = f.readlines()

with open("doc/index.rst","w") as f:
    for i in range(0,len(index)):
        if "Welcome to REBOUNDx" in index[i]:
            index[i] = "Welcome to REBOUNDx ("+reboundxversion+")\n"
            underline = ""
            for j in range(len(index[i])-1):
                underline += "="
            underline += "\n"
            index[i+1] = underline
        f.write(index[i])

with open("README.rst") as f:
    readme = f.readlines()

with open("README.rst","w") as f:
    for i in range(0,len(readme)):
        if "badge/REBOUNDx-v" in readme[i]:
            readme[i] = ".. image:: https://img.shields.io/badge/REBOUNDx-v"+reboundxversion+"-green.svg?style=flat\n"
        f.write(readme[i])

with open("src/core.c") as f:
    reboundxlines = f.readlines()
    for i,l in enumerate(reboundxlines):
        if "**VERSIONLINE**" in l:
            reboundxlines[i] = "const char* rebx_version_str = \""+reboundxversion+"\";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.\n"

    with open("src/core.c", "w") as f:
        f.writelines(reboundxlines)

with open("setup.py") as f:
    setuplines = f.readlines()
    for i,l in enumerate(setuplines):
        if "version='" in l:
            setuplines[i] = "    version='"+reboundxversion+"',\n"
        if "GITHASHAUTOUPDATE" in l:
            setuplines[i] = "    ghash_arg = \"-DREBXGITHASH="+ghash+"\" #GITHASHAUTOUPDATE\n"

    with open("setup.py", "w") as f:
        f.writelines(setuplines)

shortversion = reboundxversion
while shortversion[-1] != '.':
    shortversion = shortversion[:-1]
    
shortversion = shortversion[:-1]

with open("doc/conf.py") as f:
    conflines = f.readlines()
    for i,l  in enumerate(conflines):
        if "version =" in l:
            conflines[i] = "version = '"+shortversion+"'\n"
        if "release =" in l:
            conflines[i] = "release = '"+reboundxversion+"'\n"

    with open("doc/conf.py", "w") as f:
        f.writelines(conflines)
