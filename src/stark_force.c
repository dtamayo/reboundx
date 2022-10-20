#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

void rebx_stark_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    particles[0].ax += 0.01;
}