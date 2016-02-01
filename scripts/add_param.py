#!/usr/local/bin/python

import sys

if len(sys.argv) is not 3:
    print("Error\n*****\nMust pass name of parameter and description, e.g. python add_param.py tau_a 'semimajor axis decay timescale'")
    sys.exit()

param_name = str(sys.argv[1])
param_desc = str(sys.argv[2])

tab = "    " # 4 space tab
with open("../src/reboundx.h") as f:
    rebxh = f.readlines()

with open("../src/reboundx.h", "w") as f:
    for i in range(len(rebxh)):
        f.write(rebxh[i])
        if "enum REBX_PARAMS{" in rebxh[i]:
            f.write(tab+"{0},".format(param_name.upper())+tab+tab+tab+tab+tab+tab+tab+tab+"// {0}\n".format(param_desc))
        if "Getter setter landmark" in rebxh[i]:
            f.write("void rebx_set_{0}(struct reb_particle* p, double value);\n".format(param_name))
            f.write("double rebx_get_{0}(struct reb_particle* p);\n".format(param_name))

with open("../src/reboundx.c") as f:
    rebxc = f.readlines()

with open("../src/reboundx.c", "w") as f:
    for i in range(len(rebxc)):
        f.write(rebxc[i])
        if "Getter setter landmark" in rebxc[i]:
            f.write("void rebx_set_{0}(struct reb_particle* p, double value){{\n".format(param_name))
            f.write(tab+"double* {0}Ptr = rebx_search_param(p, {1});\n".format(param_name, param_name.upper()))
            f.write(tab+"if({0}Ptr == NULL){{\n".format(param_name))
            f.write(tab+tab+"rebx_add_param_double(p, {0}, value);\n".format(param_name.upper()))
            f.write(tab+"}\n")
            f.write(tab+"else{\n")
            f.write(tab+tab+"*{0}Ptr = value;\n".format(param_name))
            f.write(tab+"}\n")
            f.write("}\n\n")
        
            f.write("double rebx_get_{0}(struct reb_particle* p){{\n".format(param_name))
            f.write(tab+"double* {0}Ptr = rebx_search_param(p, {1});\n".format(param_name, param_name.upper()))
            f.write(tab+"if({0}Ptr == NULL){{\n".format(param_name))
            f.write(tab+tab+"return 0.;\n")
            f.write(tab+"}\n")
            f.write(tab+"else{\n")
            f.write(tab+tab+"return *{0}Ptr;\n".format(param_name))
            f.write(tab+"}\n")
            f.write("}\n\n")

with open("../reboundx/extras.py") as f:
    extras = f.readlines()

with open("../reboundx/extras.py", "w") as f:
    for i in range(len(extras)):
        f.write(extras[i])
        if "def add_Particle_props" in extras[i]:
            f.write(tab+tab+"@property\n")
            f.write(tab+tab+"def {0}(self):\n".format(param_name))
            f.write(tab+tab+"tab+clibreboundx.rebx_get_{0}.restype = c_double\n".format(param_name))
            f.write(tab+tab+"tab+return clibreboundx.rebx_get_{0}(byref(self))\n".format(param_name))
            f.write(tab+tab+"@{0}.setter\n".format(param_name))
            f.write(tab+tab+"def {0}(self, value):\n".format(param_name))
            f.write(tab+tab+tab+"clibreboundx.rebx_set_{0}(byref(self), c_double(value))\n".format(param_name))
        if "Monkeypatch landmark" in extras[i]:
            f.write(tab+tab+"rebound.Particle.{0} = {1}\n".format(param_name, param_name))

with open("../doc/modules.rst") as f:
    mod = f.readlines()

with open("../doc/modules.rst", "w") as f:
    for i in range(len(mod)):
        f.write(mod[i])
        if ".. Parameters (marker for add_param.py)" in mod[i]:
           pass 
