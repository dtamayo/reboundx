#include "rebound.h"
#include "reboundx.h"

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    struct reb_particle star = {0};
    star.m     = 1.;
    reb_add(sim, star);

    struct reb_particle planet = {0};  // add a planet on a circular orbit (with default units where G=1)
    planet.x = 1.;
    planet.vy = 1.;
    reb_add(sim, planet);

    struct rebx_extras* rebx = rebx_attach(sim);  // first initialize rebx
    struct rebx_force* stark = rebx_load_force(rebx, "stark_force"); // add our new force
    rebx_add_force(rebx, stark);

    double tmax = 100000.;
    reb_integrate(sim, tmax);
}
