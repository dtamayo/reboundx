/**
 * @file    inner_disk_edge.c
 * @brief   Inner disk edge implemention at a chosen location, while planets are undergoing migration
 * @author  Kaltrina Kajtazi <1kaltrinakajtazi@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Dan Tamayo, Hanno Rein
 * 
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The section after the dollar signs gets built into the documentation by a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 * 
 * $Orbit Modifications$       // Effect category 
 * 
 * ======================= ================================================================================================================
 * Authors                 Kajtazi, Kaltrina and D. Petit, C. Antoine
 * Implementation Paper    Kajtazi et al. in prep.
 * Based on                `Pichierri et al 2018 <https://ui.adsabs.harvard.edu/abs/2018CeMDA.130...54P/abstract>`_.
 * C example               :ref:`c_example_inner_disk_edge`
 * Python example          `InnerDiskEdge.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/InnerDiskEdge.ipynb>`_.
 * ======================= ================================================================================================================
 * 
 * This applies an inner disk edge that functions as a planet trap. Within its width the planet's migration is reversed by an opposite and roughly equal magnitude torque. Thus, stopping further migration and trapping the planet within the width of the trap. 
 * The functions here provide a way to modify the tau_a timescale in modify_orbits_forces, modify_orbit_direct, and type_I_migration.
 * Note that the present prescription is very useful for simple simulations when an inner trap is needed during the migration but it shouldn't be considered as a realistic model of the inner edge of a disk.
 * 
 * **Effect Parameters**
 * 
 * ============================ =========== ===================================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ===================================================================================
 * ide_position (double)        Yes         The position of the inner disk edge in code units 
 * ide_width (double)           Yes         The disk edge width (planet will stop within ide_width of ide_position)
 * ============================ =========== ===================================================================================
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"

/* A planet trap that is active only at the inner disk edge, to reverse the planetary migration and prevent migration onto the star */
const double rebx_calculate_planet_trap(const double r, const double dedge, const double hedge){
    double tau_a_red;

    if (r > dedge*(1.0 + hedge)){
        tau_a_red = 1.0;
    }

    else if (dedge*(1.0 - hedge) < r){
        tau_a_red =  5.5 * cos( ((dedge * (1.0 + hedge) - r ) * 2 * M_PI) / (4 * hedge * dedge) ) - 4.5;
    }

    else {
        tau_a_red = -10.0;
    }

    return tau_a_red;
}
