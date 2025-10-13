#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rebound.h"

char TITLE[100] = "cham_simarchive_";

const double min_frag_mass = 1.4e-8; //should not be global
const double rho1 = 1.684e6;  //Msun/AU^3 
const double cstar = 1.8;      //need to add proper comments and defs


#define MIN(a, b) ((a) > (b) ? (b) : (a))    // Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    // Returns the maximum of a and b

//define collision parameters
struct collision_params
{
    int target;
    int projectile;
    double dx;
    double dy;
    double dz;
    double b;
    double dvx;
    double dvy;
    double dvz;
    double Vi;
    double l;
    double mu;
    double QR;
    double QpRD;
    double V_esc;
    double separation_distance;
    double Mlr;
    double Mslr;
    double Q;
    double Mlr_dag;
    double Q_star;
    double vrel;
    double xrel;
}; 

//Get dot product of two vectors
double get_dot(double x1, double y1, double z1, double x2, double y2, double z2){ 
    return (x1*x2)+(y1*y2)+(z1*z2);
}  

//Get magnitude of vector
double get_mag(double x, double y, double z){
    return sqrt(pow(x,2)+pow(y,2)+pow(z,2));
}  

//Get radius based on mass and density
double get_radii(double m, double rho){
    return pow((3.*m)/(4.*M_PI*rho),1./3.);
}   


/* add_fragments is the function to compute number and mass of fragments
   it also adds fragments to the simulation
   and sets their location and velocity
   */
void add_fragments(struct reb_simulation* const r, struct reb_collision c, const struct collision_params *const params){
    struct reb_particle* target = &(r->particles[params->target]);
    struct reb_particle* projectile = &(r->particles[params->projectile]);
    struct reb_particle com = reb_particle_com_of_pair(*target, *projectile);
    double initial_mass = target -> m + projectile -> m;
    double remaining_mass = initial_mass - params->Mlr;
    double rho = target->m/(4./3*M_PI*pow(target ->r, 3));
    double rtot = target -> r + projectile -> r;

    //In hit and run collisions, there is a second largest remnant or a "big fragment"
    //If mass of second largest remnant (Mslr) is bigger than zero, then number of big_frags is 1
    //and fragments will be made of what remains
    int big_frags = 0;
    if (params->Mslr > 0){
        remaining_mass = remaining_mass -  params->Mslr;
        big_frags = 1;
    }

    //Fragments are broken into equal sizes
    int small_frags = remaining_mass/min_frag_mass;
    double frag_mass = remaining_mass/small_frags;

    //New bodies introduced into the simulation,
    //Might be just fragments, or fragments + second largest remnant from a hit and run collision (big frag)
    int n_frags = small_frags + big_frags; //TODO do we need this extra param?

    double mxsumx = 0;
    double mxsumy = 0;
    double mxsumz = 0;

    double mvsumx = 0;
    double mvsumy = 0;
    double mvsumz = 0;


    //target gets mass of Mlr and is assigned COM position and velocity;
    target -> last_collision = r->t;
    target -> m = params->Mlr;
    target -> r = get_radii(params->Mlr, rho);
    target->x = com.x;
    target->y = com.y;
    target->z = com.z;

    target->vx = com.vx;
    target->vy = com.vy;
    target->vz = com.vz;

    //If Mlr is less than the computed fragment mass, then Mlr and fragment are swaped (//TODO not sure if this is okay)
    if (small_frags == 1 && params->Mlr <= frag_mass){
        target->m = frag_mass;
        frag_mass = params->Mlr;
    }

    mxsumx = mxsumx + target->m*target->x;  //adding target mx, to check for com offset later
    mxsumy = mxsumy + target->m*target->y;  //same with y  
    mxsumz = mxsumz + target->m*target->z;  //same with z

    mvsumx = mvsumx + target->m*target->vx; //adding target mv, to check for momentum offset
    mvsumy = mvsumy + target->m*target->vy; //same with y
    mvsumz = mvsumz + target->m*target->vz; //same with z

    double theta_inc = (2.*M_PI)/n_frags; //seperation angle between fragments


    double unit_dvx, unit_dvy, unit_dvz, zx, zy, zz, z, ox, oy, oz, o; //TODO not sure if all of these are needed

    //unit vectors parallel to target velocity
    unit_dvx = params->dvx/params->vrel;
    unit_dvy = params->dvy/params->vrel;
    unit_dvz = params->dvz/params->vrel;

    //vectors normal to the collision plane, vrel cross xrel
    //TODO what is this doing?
    zx = (params->dvy*params->dz - params->dvz*params->dy);                 
    zy = (params->dvz*params->dx - params->dvx*params->dz);
    zz = (params->dvx*params->dy - params->dvy*params->dx);

    //get magnitude of vrel cross vrel
    z = get_mag(zx, zy, zz);

    //TODO make them unit vectors ?
    zx = zx/z;          
    zy = zy/z;
    zz = zz/z;


    ox = (zy*params->dvz - zz*params->dvy);                   //TODO (??) vector normal to target velocity in collision plane; z cross vrel
    oy = (zz*params->dvx - zx*params->dvz);
    oz = (zx*params->dvy - zy*params->dvx);

    o = get_mag(ox, oy, oz);

    ox = ox/o;      //unit vector               //TODO no idea what this part is doing
    oy = oy/o;
    oz = oz/o;

    //Based on Chambers et al. 2013, fragments leave at a velocity %5 higher than mutual escape velocity
    //In such a way to conserve momentum
    //Separation distance is set as 4*R_tot (//TODO didn't find this in Chambers)
    //TODO why is fragment velocity not simply 1.05v_esc? Isn't G*m_init/rtot = vesc^2??
    double fragment_velocity =sqrt(1.1*pow(params->V_esc,2) - 2*r->G*initial_mass*(1./rtot - 1./params->separation_distance));

    if (big_frags == 1){  //assign radii, positions and velocities to second largest remnant, theta=0
        struct reb_particle Slr1 = {0};
        Slr1.m = params->Mslr;
        Slr1.x = com.x + params->separation_distance*unit_dvx;  //TODO I think this makes sense but need to double check
        Slr1.y = com.y + params->separation_distance*unit_dvy;  //TODO why unit v?? 
        Slr1.z = com.z + params->separation_distance*unit_dvz;

        Slr1.vx = com.vx + fragment_velocity*unit_dvx;
        Slr1.vy = com.vy + fragment_velocity*unit_dvy;
        Slr1.vz = com.vz + fragment_velocity*unit_dvz;

        Slr1.r = get_radii(Slr1.m, rho); //TODO this uses target's density, do we want to be able to change densities?

        //Center of mass offset parameter
        mxsumx += Slr1.m*Slr1.x;
        mxsumy += Slr1.m*Slr1.y;    
        mxsumz += Slr1.m*Slr1.z;

        //Momentum offset parameter
        mvsumx += Slr1.m*Slr1.vx;
        mvsumy += Slr1.m*Slr1.vy;   
        mvsumz += Slr1.m*Slr1.vz;

        //Record the collision time as last collision
        Slr1.last_collision = r->t;

        //Add the second largest remnant to simulation
        reb_simulation_add(r, Slr1);
    }

    //I think this is probably for unique hash generation, but isn't there a better way to make hashes?
    for (int j=1; j < small_frags + 1; j++){          //add fragments
        struct reb_particle fragment = {0};

        fragment.m = frag_mass;                  
        fragment.x = com.x + params->separation_distance*(cos(theta_inc*j)*unit_dvx + sin(theta_inc*j)*ox); //TODO I don't fully understand this part
        fragment.y = com.y + params->separation_distance*(cos(theta_inc*j)*unit_dvy + sin(theta_inc*j)*oy);
        fragment.z = com.z + params->separation_distance*(cos(theta_inc*j)*unit_dvz + sin(theta_inc*j)*oz);
        fragment.vx = com.vx + fragment_velocity*(cos(theta_inc*j)*unit_dvx + sin(theta_inc*j)*ox);
        fragment.vy = com.vy + fragment_velocity*(cos(theta_inc*j)*unit_dvy + sin(theta_inc*j)*oy);
        fragment.vz = com.vz + fragment_velocity*(cos(theta_inc*j)*unit_dvz + sin(theta_inc*j)*oz);

        //Computes fragment radius based on its mass and target's density //TODO maybe we want to make an average density?
        fragment.r = get_radii(frag_mass, rho);

        //Records collision
        fragment.last_collision = r->t;


        //COM offset
        mxsumx +=fragment.m*fragment.x;
        mxsumy += fragment.m*fragment.y;    
        mxsumz += fragment.m*fragment.z;

        //Momentum offset
        mvsumx += fragment.m*fragment.vx;
        mvsumy += fragment.m*fragment.vy;    
        mvsumz += fragment.m*fragment.vz;

        reb_simulation_add(r, fragment); 
    }


    //Ensure momentum is conserved (//TODO ?)
    //Compute offset for center of mass
    double xoff[3] = {com.x - mxsumx/initial_mass, com.y - mxsumy/initial_mass, com.z - mxsumz/initial_mass};

    //Compute offset for momentum
    double voff[3] = {com.vx - mvsumx/initial_mass, com.vy - mvsumy/initial_mass, com.vz - mvsumz/initial_mass};

    //Correct for COM and momentum offset
    target -> x +=  xoff[0]*target->m/initial_mass; //TODO fixing for offset of momentum, need to check this
    target -> y += xoff[1]*target->m/initial_mass; 
    target -> z += xoff[2]*target->m/initial_mass; 
    target -> vx += voff[0]*target->m/initial_mass; 
    target -> vy += voff[1]*target->m/initial_mass; 
    target -> vz += voff[2]*target->m/initial_mass; 

    //Add fragments into the simulation
    for (int i = r->N - n_frags; i < r->N; i++){ 
        double mass_fraction = r->particles[i].m/initial_mass; // cm: what is this for?
        r->particles[i].x += xoff[0]*mass_fraction;
        r->particles[i].y += xoff[1]*mass_fraction;
        r->particles[i].z += xoff[2]*mass_fraction;

        r->particles[i].vx += voff[0]*mass_fraction;
        r->particles[i].vy += voff[1]*mass_fraction;
        r->particles[i].vz += voff[2]*mass_fraction;
    }

    return;
}


/* Function to merge two particles
   This happens when v_i < v_esc
   Keeps target, gets rid of projectile
   Conserves mass and momentum
   */
void merge(struct reb_simulation* const r, struct reb_collision c, const struct collision_params *const params){
    struct reb_particle* pi = &(r->particles[params->target]); //sets particle pi as target
    struct reb_particle* pj = &(r->particles[params->projectile]); //sets particle pj as projectile

    double invmass = 1.0/(pi->m + pj->m);
    double targ_rho = pi->m/(4./3*M_PI*pow(pi->r,3));  //new body recieves density of the target // TODO why? I think it should be their mean density
                                                       // Merge by conserving mass, volume and momentum
    pi->vx = (pi->vx*pi->m + pj->vx*pj->m)*invmass;     //TODO why not make this a vector to ease our lives?
    pi->vy = (pi->vy*pi->m + pj->vy*pj->m)*invmass;
    pi->vz = (pi->vz*pi->m + pj->vz*pj->m)*invmass;     //TODO same here
    pi->x  = (pi->x*pi->m + pj->x*pj->m)*invmass;
    pi->y  = (pi->y*pi->m + pj->y*pj->m)*invmass;
    pi->z  = (pi->z*pi->m + pj->z*pj->m)*invmass;
    pi->m  = pi->m + pj->m;
    pi->r  = pow((3*pi->m)/(4*M_PI*targ_rho),1./3.);    //TODO need to make this combined density
    pi->last_collision = r->t;

    return;
}


/* Function to compute collision parameters for hit-and-run regime
   hit-and-run collisions are the ones where sin(theta) < R_t/(R_t + R_p) (or b < b_crit)
   In this regime, if M_lr > M_t, then we will have M_lr, M_slr, and fragments
   Where M_lr is mass of the largest remnant, M_slr is mass of second largest remnant, and fragments are smaller equal mass remnants
   This means that projectile is eroded, giving some of its mass to target, and then breaks into M_slr and fragments
   If M_lr < M_t collision is resolved similar to an erosion.
   */
int hit_and_run(struct reb_simulation* const r, struct reb_collision c, struct collision_params *params){  //also includes partial accretion.  Mlr = M_target.  Projectile is erroded.
    struct reb_particle* target = &(r->particles[params->target]); 
    struct reb_particle* projectile = &(r->particles[params->projectile]);
    double rho_t = target->m/(4./3*M_PI*pow(target ->r, 3)); //density of the target

    //To make sure projectile is removed not the target,
    //need to check for masses of colliding particles and set swap accordingly
    int i = c.p1;   //TODO why is merge not doing this?
    int j = c.p2;   
    struct reb_particle* pi = &(r->particles[i]);
    struct reb_particle* pj = &(r->particles[j]);

    int remove;  
    if (pi->m < pj->m){ //makes sure projectile is removed
        remove = 1;
    }else{
        remove = 2;
    }

    double phi = 2*acos((params->l-projectile->r)/projectile->r); //Helps with finding part of the projectile that is NOT crossing the target
    double A_interact = pow(projectile->r, 2)*((M_PI-(phi-sin(phi))/2.));  //Leinhardt Eq. 46; cross section of projectile interacting with the target
    double L_interact = 2.*pow(pow(target->r,2)-(pow(target->r-params->l/2.,2)), .5);   //Leinhardt Eq. 47, interacting length

    //TODO Changed line below, multiplied by density to get a dimensionless value for beta
    double beta = ((A_interact*L_interact) * rho_t)/target->m;  //Leinhardt Eq. 48, used in Chambers Eq. 11
    double Rc1 = pow(3./(4.*M_PI*rho1)*(beta*target->m + projectile->m), 1./3.); //Based on Chambers Eq. 11
    double Q0 = .8*cstar*M_PI*rho1*r->G*pow(Rc1, 2); //Chambers Eq. 11
    double gamma = (beta*target->m)/projectile->m; //Based on Chambers Eq. 11
    double Q_star = (pow(1+gamma, 2)/4*gamma)* Q0; //Chambers Eq. 10

    double mu = (beta*target->m*projectile->m)/(beta*target->m+projectile->m);  //Chambers Eq. 13
    double Q = .5*(mu*pow(params->Vi,2))/(beta*target->m+projectile->m); //Chambers Eq. 12

    /* If  velocity in the hit-and-run regime is very low, the collision
     * might eventually lead to a merger. Here, we compute the threshhold velocity for this event,
     * called critical velocity. If v < v_crit, then we have a "graze and merge" event.
     */
    double c1 = 2.43; //c1 to c4 are constants used in Chambers Eq. 17
    double c2 = -0.0408;
    double c3 = 1.86;
    double c4 = 1.08;


    double zeta = pow((target->m - projectile->m)/(target->m + projectile->m),2); //TODO in Chambers, this is [(1-gamma)^2/(1+gamma)^2]^2, maybe we should use the same gamma as above?
    double fac = pow(1-params->b/(target->r + projectile->r),2.5); //TODO maybe can re-write this prettier
    double v_crit = params->V_esc*(c1*zeta*fac + c2*zeta +c3*fac + c4); // Velocity boundary between graze and merge, Chambers Eq. 15

    if (params->Vi <= v_crit){             //Graze and merge, impact velocity less than critical velocity         
        merge(r,c,params);
        return remove;
    }else{ //vi > v_crit
        if (params->Mlr < target->m){ //Resolved as an erosion event
            if (target->m + projectile->m <= 2*min_frag_mass){ //not enough mass to produce new fragments
                reb_collision_resolve_hardsphere(r,c);
                return 0; //remove no particle
            }else{
                params->Mlr = MAX(params->Mlr, min_frag_mass); //If Mlr is less than min frag mass, then it is set to min frag mass
                                                               //TODO shouldn't we interpret this also as elastic bounce? seems like we are adding mass to simulation in this way!
                add_fragments(r,c,params);
            }
        }else{ //Mlr > Mt, either a hit and run or an elastic bounce
            double Mlr_dag = (beta*target->m + projectile->m)/10 * pow(Q/(1.8*Q_star), -1.5);
            if (Q < 1.8*Q_star){
                Mlr_dag = (beta*target->m + projectile->m)*(1 - Q/ (2*Q_star));
            }

            double projectile_mass_accreted = params->Mlr - target->m;
            double new_projectile_mass = projectile->m - projectile_mass_accreted;
            Mlr_dag = MAX(Mlr_dag, min_frag_mass); //TODO again, if Mlr_dag is less than min frag mass, 
                                                   // you set it as min frag mass, so introducing additional mass to system
            if (new_projectile_mass-Mlr_dag < min_frag_mass){ //fragment mass lower than resolution
                reb_collision_resolve_hardsphere(r,c);
                return 0; //remove no particle
            }else{ //Hit-and-run event (partial accretion)
                params->Mslr = Mlr_dag;
                add_fragments(r,c,params);
            }
        }
        return remove;
    }
}


int reb_collision_resolve_fragment(struct reb_simulation* const r, struct reb_collision c){
    if (r->particles[c.p1].last_collision==r->t || r->particles[c.p2].last_collision==r->t) return 0;
    if (c.p1 < c.p2) return 0;      //only return one collision callback

    int i;
    int j; 
    int remove;


    if (r->particles[i].m < r->particles[j].m){ //unless remove is redfined as 0, projectile is going to be removed.
        remove =1;
        i = c.p2; //i is the higher mass object (target)
        j = c.p1; //j is the lower mass object (projectile)
    }else{
        remove = 2;
        i = c.p1;
        j = c.p2;
    }

    struct reb_particle* particles = r->particles;
    struct collision_params params = {0};

    double imp_r = particles[j].r; //radius of projectile
    double targ_r = particles[i].r; //radius of target
    double R_tot = imp_r + targ_r; //Sum of radii

    double imp_m = particles[j].m; //Projectile mass
    double targ_m = particles[i].m; //Target mass

    double M_tot = imp_m + targ_m; //Total mass
    double G = r->G;
    double Mlr,dx,dy,dz,dvx,dvy,dvz;
    double x2rel, xrel, v2rel, v2imp, Vi;
    double hx,hy,hz,h2,b;

    dx = particles[i].x - particles[j].x;
    dy = particles[i].y - particles[j].y;
    dz = particles[i].z - particles[j].z;

    x2rel = get_dot(dx, dy, dz, dx, dy, dz); //distance between two particles ^ 2

    dvx = particles[i].vx - particles[j].vx;
    dvy = particles[i].vy - particles[j].vy;
    dvz = particles[i].vz - particles[j].vz;

    v2rel = get_dot(dvx,dvy,dvz,dvx,dvy,dvz); //relative velociry between two particles ^ 2

    xrel = sqrt(x2rel);  //distance between the centers of the projectile and target


    hx = (dy*dvz - dz*dvy);                     //angular momentum vector xrel X Vrel  //TODO ?
    hy = (dz*dvx - dx*dvz);
    hz = (dx*dvy - dy*dvx);

    h2 = get_dot(hx,hy,hz,hx,hy,hz);

    v2imp = v2rel + 2*G*M_tot*(1./R_tot - 1./xrel); //impact velocity with gravitational focusing at time of detected collision

    if (1./R_tot - 1./xrel < 0){v2imp = v2rel;}  //if collision is detected after physical contact

    Vi = sqrt(v2imp);  //magnitude of impact velocity vector
    b = sqrt(h2/v2imp);  //impact parameter, b=R_tot*sin(theta)
    if (isnan(b)){
        reb_simulation_error(r, "b is not a number.");
        return 0;
    }
    //Stewart & Leinhardt 2012 parameters
    double mu = (targ_m*imp_m)/M_tot;  //Chambers Eq. 2, reduced mass
    double l = R_tot-b;  //Leinhardt Eq. 7, the projected length of the projectile overlapping the target
    l = MIN(l, 2*imp_r);
    double alpha = (pow(l,2)*(3*imp_r-l))/(4*pow(imp_r, 3)); //Leinhardt Eq. 11, interacting mass fraction
    alpha = MIN(1., alpha);
    double Q = .5*v2imp*targ_m*imp_m/pow(M_tot,2);  //specific energy per unit mass
    double V_esc = pow(2.*G*M_tot/R_tot, .5); //mutal escape velocity as defined in Wallace et al 2018 and Chambers 2013
    double alphamu = (alpha*targ_m*imp_m)/(alpha*imp_m + targ_m);  //Leinhardt Eq. 12, reduced interacting mass for fraction alpha.
    double gamma = imp_m/targ_m;  //Chambers Eq. 6

    double Rc1 = pow((M_tot*3)/(4.*M_PI*rho1), 1./3.);  //Chambers Eq. 4, combined radius of target and projectile with constant density
    double Q0 = .8*cstar*M_PI*rho1*G*pow(Rc1,2);  //Chambers Eq. 3, critical value of impact energy for head-on collisions
    double Q_star = pow(mu/alphamu, 1.5)*(pow(1+gamma, 2)/ (4*gamma))*Q0;  //Chambers Eq. 5, critical value for oblique or different mass collisons.  
    if (alpha == 0.0){
        reb_simulation_error(r, "alpha = 0");
        return 0;
    }
    if (b == 0 && imp_m == targ_m){
        Q_star = Q0;
    }
    double qratio = Q/Q_star;
    if (qratio < 1.8){
        Mlr = M_tot*(1.0-.5*qratio);
    }else{
        Mlr = .1*M_tot*pow(qratio/1.8, -1.5);  //Chambers Eq.8
    }

    double separation_distance = 4 * R_tot;  //seperation distance of fragments.  Should be userdefined but this is what chambers uses //cm: wasn't this based on number of frags?
                                             ///POPULATE STRUCT OBJECTS
    params.target = i;
    params.projectile = j;
    params.dx = dx;
    params.dy = dy;
    params.dz = dz;
    params.b = b;
    params.dvx = dvx;
    params.dvy = dvy;
    params.dvz = dvz;
    params.Vi = Vi;
    params.l = l;
    params.mu = mu;
    params.Q = Q;
    params.separation_distance = separation_distance;
    params.V_esc = V_esc;
    params.vrel = sqrt(v2rel);
    params.Mslr = 0;
    params.xrel = xrel;
    params.Mlr = Mlr;

    if (Vi <= V_esc){
        merge(r,c, &params);
    }else{  //Vi > V_esc
        if (b<targ_r){ //non-grazing regime
            if (M_tot - params.Mlr < min_frag_mass){
                merge(r,c,&params);
            }else{ // M_tot - params->Mlr >= min_frag_mass; fragments will be produced unless it is a graze and merge or elastic bounce 
                if (params.Mlr < targ_m){
                    if (params.Mlr <= 0.1*targ_m){
                        params.Mlr = MAX(Mlr, min_frag_mass);
                        add_fragments(r,c,&params);
                    }else{
                        params.Mlr = MAX(Mlr, min_frag_mass);
                        add_fragments(r,c,&params);
                    }
                }else{  //(params->Mlr >= targ_m)
                    add_fragments(r,c,&params);
                }
            }
        }else{ // b > b_crit, grazing regime
            return hit_and_run(r,c,&params); //swap gets redefined here as it may be set to 0 in the case of a bounce // cm: why?
        }
    }

    return remove;
}





int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    r->G = 39.476926421373;
    r->dt = 6./365.;
    //r->exit_max_distance = 100.; 

    // The random seed is passed as a command line argument
    if (argc == 2){
        r->rand_seed = atoi(argv[1]);
        strcat(TITLE, argv[1]);
    }
    r->integrator = REB_INTEGRATOR_MERCURIUS;
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_fragment;

    //Assigning mass and number of planetary embryos and planetesimals
    double rho = 5.05e6; //3 g/cm^3
    int n_emb = 14; //number of planetary embryos
    double m_emb = 2.8e-7; //mass of each embryo
    int n_pl = 140; //number of planetesimals
    double m_pl = 2.8e-08;
    int ef = 1;
    double a_pl[n_pl];  // Array to hold the values of planetesimal semi-major axis
    double a_emb[n_emb]; // Same for embryos

    FILE *file_pl = fopen("/home/tajer.1/planetgarten/initial_disk_setup/output_semis/semis_pl_chambers.txt", "r");
    for (int i = 0; i < n_pl; i++) {
        fscanf(file_pl, "%lf", &a_pl[i]);
    }
    // Close the file after reading
    fclose(file_pl);

    FILE *file_emb = fopen("/home/tajer.1/planetgarten/initial_disk_setup/output_semis/semis_emb_chambers.txt", "r");
    for (int i = 0; i < n_emb; i++) {
        fscanf(file_emb, "%lf", &a_emb[i]);
    }
    // Close the file after reading
    fclose(file_emb);


    struct reb_particle star = {0};
    star.m = 1.00;
    star.r = 0.1; 
    star.hash = 0; 
    reb_simulation_add(r, star);

    // Add planetary embryos

    for (int i=0; i<n_emb; i++){
        double a = a_emb[i];      // semi major axis
        double m = m_emb;
        double inc = reb_random_uniform(r,0,0.0175);
        double ecc = reb_random_uniform(r,0,0.01);
        double omega = reb_random_uniform(r,0,2*M_PI);
        double Omega = reb_random_uniform(r,0,2*M_PI);
        double f = reb_random_uniform(r,0,2*M_PI);
        double hash = i + 1;
        //now build particle from orbit
        struct reb_particle emb = reb_particle_from_orbit(r->G, star, m, a, ecc, inc, Omega, omega, f);
        emb.r = get_radii(m, rho)*ef;
        emb.hash = hash;

        reb_simulation_add(r, emb); 
    }
    //add planetesimals
    for (int i=0; i<n_pl; i++){
        double a = a_pl[i];      // semi major axis
        double m = m_pl;
        double inc = reb_random_uniform(r,0,0.0175);
        double ecc = reb_random_uniform(r,0,0.01);
        double omega = reb_random_uniform(r,0,2*M_PI);
        double Omega = reb_random_uniform(r,0,2*M_PI);
        double f = reb_random_uniform(r,0,2*M_PI);
        double hash = n_emb + i + 1;
        //now build particle from orbit
        struct reb_particle pl = reb_particle_from_orbit(r->G, star, m, a, ecc, inc, Omega, omega, f);
        pl.r = get_radii(m, rho)*ef;
        pl.hash = hash;

        reb_simulation_add(r, pl); 
    }

    fclose(of_dbcti);


    double m,a,e,inc,Omega,omega,f; //Omega=longitude of ascending node, omega= argument of pericenter in RADIANS

    //Add Jupiter and Saturn
    struct reb_particle Jup = reb_particle_from_orbit(r->G, r->particles[0], m=9.543e-4, a=5.20349, e=0.048381, inc=0.365*(M_PI/180), Omega=0.0, omega=68.3155*(M_PI/180), f=227.0537*(M_PI/180));
    Jup.r = get_radii(Jup.m, rho); 
    Jup.hash = reb_hash("JUPITER");
    reb_simulation_add(r,Jup);

    struct reb_particle Sat = reb_particle_from_orbit(r->G, r->particles[0], m=0.0002857, a=9.54309, e=0.052519, inc=0.8892*(M_PI/180), Omega=M_PI, omega=324.5263*(M_PI/180),f=256.9188*(M_PI/180));
    Sat.r = get_radii(Sat.m, rho);
    Sat.hash = reb_hash("SATURN");
    reb_simulation_add(r,Sat);


    reb_simulation_move_to_com(r);  // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    double run_time = 10;
    reb_simulation_save_to_file_interval(r,TITLE,1.e5);
    reb_simulation_integrate(r, run_time);
}



