.. _getters:

Getters/Setters available in REBOUNDx
=====================================

There needs to be a separate ``rebx_get_param_type`` and ``rebx_set_param_type`` for each (C) variable type you need to use.  

To add a new one, copy the functions for type ``double`` in ``reboundx/src/core.c`` directly below the ones under the heading `Getters and Setters for particle parameters (need new set for each variable type)`, and modify them for your new variable type.  
Contact me for help with this if needed.

You will also need to decide on a default return value if the parameter is not set for the particular particle, so that different implementations can check for this.
You can find the default value returned when a parameter is not found for each variable type below.

Getters and setters are currently implemented for the following variable types:

=================== =================================
Type                Return Value if Parameter Not Set
=================== =================================
double              nan
=================== =================================
