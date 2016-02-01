#!/usr/local/bin/python

import sys

if len(sys.argv) is not 3:
    print("Error\n*****\nMust pass name of parameter and description, e.g. python add_param.py tau_a 'semimajor axis decay timescale'")
    sys.exit()

param_name = str(sys.argv[1])
param_desc = str(sys.argv[2])

with open("../src/reboundx.h") as f:
    rebxh = f.readlines()

with open("../src/reboundx.h", "w") as f:
    for i in range(len(rebxh)):
        f.write(rebxh[i])
        if "enum REBX_PARAMS{" in rebxh[i]:
            f.write("\t{0},\t\t\t\t\t\t\t\t// {1}\n".format(param_name.upper(), param_desc))
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
            f.write("\tdouble* {0}Ptr = rebx_search_param(p, {1});\n".format(param_name, param_name.upper()))
            f.write("\tif({0}Ptr == NULL){{\n".format(param_name))
            f.write("\t\trebx_add_param_double(p, {0}, value);\n".format(param_name.upper()))
            f.write("\t}\n")
            f.write("\telse{\n")
            f.write("\t\t*{0}Ptr = value;\n".format(param_name))
            f.write("\t}\n")
            f.write("}\n\n")
        
            f.write("double rebx_get_{0}(struct reb_particle* p){{\n".format(param_name))
            f.write("\tdouble* {0}Ptr = rebx_search_param(p, {1});\n".format(param_name, param_name.upper()))
            f.write("\tif({0}Ptr == NULL){{\n".format(param_name))
            f.write("\t\treturn 0.;\n")
            f.write("\t}\n")
            f.write("\telse{\n")
            f.write("\t\treturn *{0}Ptr;\n".format(param_name))
            f.write("\t}\n")
            f.write("}\n\n")

with open("../reboundx/extras.py") as f:
    extras = f.readlines()

with open("../reboundx/extras.py", "w") as f:
    for i in range(len(extras)):
        f.write(extras[i])
        if "def add_Particle_props" in extras[i]:
            f.write("\t\t@property\n")
            f.write("\t\tdef {0}(self):\n".format(param_name))
            f.write("\t\t\tclibreboundx.rebx_get_{0}.restype = c_double\n".format(param_name))
            f.write("\t\t\treturn clibreboundx.rebx_get_{0}(byref(self))\n".format(param_name))
            f.write("\t\t@{0}.setter\n".format(param_name))
            f.write("\t\tdef {0}(self, value):\n".format(param_name))
            f.write("\t\t\tclibreboundx.rebx_set_{0}(byref(self), c_double(value))\n".format(param_name))
        if "Monkeypatch landmark" in extras[i]:
            f.write("\t\trebound.Particle.{0} = {1}\n".format(param_name, param_name))

with open("../doc/modules.rst") as f:
    mod = f.readlines()

with open("../doc/modules.rst", "w") as f:
    for i in range(len(mod)):
        f.write(mod[i])
        if ".. Parameters (marker for add_param.py)" in mod[i]:
           pass 
