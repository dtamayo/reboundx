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
#include "core.h"

void rebx_write_param(struct rebx_param* param, FILE* of){
    long pos_rewrite = ftell(of);
    struct rebx_binary_field field = {.type = REBX_BINARY_FIELD_TYPE_PARAM, .size=0};
    fwrite(&field, sizeof(field), 1, of);
    
    long pos_start_effect = ftell(of);
    size_t namelength = strlen(param->name) + 1; // +1 for \0 at end
    //printf("bef name %lu\n", ftell(of));
    fwrite(&namelength, sizeof(namelength), 1, of);
    printf("before name %lu\n", ftell(of));
    fwrite(param->name, sizeof(*param->name), namelength, of);
    fwrite(&param->param_type, sizeof(param->param_type), 1, of);
    fwrite(&param->ndim, sizeof(param->ndim), 1, of);
    printf("before shape %lu\n", ftell(of));
    fwrite(param->shape, sizeof(*param->shape), param->ndim, of);
    printf("after shape %lu\n", ftell(of));
    fwrite(param->contents, rebx_sizeof(param->param_type), param->size, of);
    printf("after contents %lu\n", ftell(of));
    
    /*double a = 8;
    double b= 6.;
    fwrite(&a, sizeof(a), 1, of);
    fwrite(&b, sizeof(b), 1, of);
    */
    long pos_end_effect = ftell(of);
    field.size = pos_end_effect - pos_start_effect;
    //printf("%lu\t%lu\t%lu\n", pos_start_effect, pos_end_effect, field.size);
    
    fseek(of, pos_rewrite, SEEK_SET);
    fwrite(&field, sizeof(field), 1, of);
    fseek(of, 0, SEEK_END);
    //printf("%lu\n", ftell(of));
}
    
void rebx_write_params(struct rebx_param* ap, FILE* of){
    int nparams=0;
    struct rebx_param* current = ap;
    while (current != NULL){
        current = current->next;
        nparams++;
    }
    
    while (nparams > 0){
        current = ap;
        for(int i=0;i<nparams-1;i++){
            current=current->next;
        }
        rebx_write_param(current, of);
        nparams--;
    }
}

void rebx_write_effect(struct rebx_effect* effect, FILE* of){
    //printf("bef effect %lu\n", ftell(of));
    long pos_rewrite = ftell(of);
    struct rebx_binary_field field = {.type = REBX_BINARY_FIELD_TYPE_EFFECT, .size=0};
    fwrite(&field, sizeof(field), 1, of);
    
    long pos_start_effect = ftell(of);
    size_t namelength = strlen(effect->name) + 1; // +1 for \0 at end
    //printf("bef name %lu\n", ftell(of));
    fwrite(&namelength, sizeof(namelength), 1, of);
    fwrite(effect->name, sizeof(*effect->name), namelength, of);
    rebx_write_params(effect->ap, of);
    
    // update size in binary_field to reflect the size of the effect now that we can calculate how big it is
    long pos_end_effect = ftell(of);
    field.size = pos_end_effect - pos_start_effect;
    //printf("%lu\t%lu\t%lu\n", pos_start_effect, pos_end_effect, field.size);
    fseek(of, pos_rewrite, SEEK_SET);
    fwrite(&field, sizeof(field), 1, of);
    fseek(of, 0, SEEK_END);
    //printf("%lu\n", ftell(of));
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
        rebx_write_effect(current, of);
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
