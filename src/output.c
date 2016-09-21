/**
 * @file    output.c
 * @brief   Output functions.
 * @author  Dan Tamayo <tamayo.daniel@gmail.com>
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
 */

#include <stdio.h>
#include <string.h>
#include "reboundx.h"

void rebx_write_effect(struct rebx_effect* effect){
    printf("%s\n", effect->name);
}

void rebx_output_binary(struct rebx_extras* rebx, char* filename){
    struct reb_simulation* sim = rebx->sim;

    FILE* of = fopen(filename,"wb"); 
    if (of==NULL){
        reb_error(sim, "REBOUNDx error: Can not open file passed to rebx_output_binary.");
    }

    // Output header.
    const char str[] = "REBOUNDx Binary File. Version: ";
    size_t lenheader = strlen(str)+strlen(rebx_version_str);
    fwrite(str,sizeof(char),strlen(str),of);
    fwrite(rebx_version_str,sizeof(char), strlen(rebx_version_str),of);
    while (lenheader<64){ //padding
        char space = ' ';
        if (lenheader==63) space = '\0';
        fwrite(&space,sizeof(char),1,of);
        lenheader += 1;
    }

    int neffects=0;
    struct rebx_effect* current = rebx->effects;
    while (current != NULL){
        current = current->next;
        neffects++;
    }
    
    while (neffects > 0){
        current = rebx->effects;
        for(int i=0;i<neffects-1;i++){
            current=current->next;
        }
        rebx_write_effect(current);
        neffects--;
    }
    fclose(of);
}
//fwrite(current, sizeof(struct rebx_effect), 1, of);
/*// Output main simulation structure
 fwrite(r,sizeof(struct reb_simulation),1,of);
 
 // Output particles
 fwrite(r->particles,sizeof(struct reb_particle),r->N,of);
 
 // Output variational configuration structures
 if (r->var_config_N){
 fwrite(r->var_config,sizeof(struct reb_variational_configuration),r->var_config_N,of);
 }
 
 // Output IAS15 temporary arrays (needed for bit-by-bit reproducability)
 if (r->ri_ias15.allocatedN){
 int N3 = r->ri_ias15.allocatedN;
 fwrite(r->ri_ias15.at,sizeof(double),N3,of);
 fwrite(r->ri_ias15.x0,sizeof(double),N3,of);
 fwrite(r->ri_ias15.v0,sizeof(double),N3,of);
 fwrite(r->ri_ias15.a0,sizeof(double),N3,of);
 fwrite(r->ri_ias15.csx,sizeof(double),N3,of);
 fwrite(r->ri_ias15.csv,sizeof(double),N3,of);
 fwrite(r->ri_ias15.csa0,sizeof(double),N3,of);
 reb_save_dp7(&(r->ri_ias15.g)  ,N3,of);
 reb_save_dp7(&(r->ri_ias15.b)  ,N3,of);
 reb_save_dp7(&(r->ri_ias15.csb),N3,of);
 reb_save_dp7(&(r->ri_ias15.e)  ,N3,of);
 reb_save_dp7(&(r->ri_ias15.br) ,N3,of);
 reb_save_dp7(&(r->ri_ias15.er) ,N3,of);
 }*/

/*
void reb_output_binary_positions(struct reb_simulation* r, char* filename){
    const int N = r->N;
#ifdef MPI
    char filename_mpi[1024];
    sprintf(filename_mpi,"%s_%d",filename,r->mpi_id);
    FILE* of = fopen(filename_mpi,"wb"); 
#else // MPI
    FILE* of = fopen(filename,"wb"); 
#endif // MPI
    if (of==NULL){
        reb_exit("Can not open file.");
    }
    for (int i=0;i<N;i++){
        struct reb_vec3d v;
        v.x = r->particles[i].x;
        v.y = r->particles[i].y;
        v.z = r->particles[i].z;
        fwrite(&(v),sizeof(struct reb_vec3d),1,of);
    }
    fclose(of);
}
*/
