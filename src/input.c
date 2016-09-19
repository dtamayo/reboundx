/**
 * @file    input.c
 * @brief   Input functions.
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

void rebx_load_effects(struct reb_simulation* sim, char* filename, enum reb_input_binary_messages* warnings){
    FILE* inf = fopen(filename,"rb"); 
    if (inf){
        struct reb_simulation* r = malloc(sizeof(struct reb_simulation));
        long objects = 0;
        // Input header.
        const char str[] = "REBOUND Binary File. Version: ";
        char readbuf[65], curvbuf[65];
        sprintf(curvbuf,"%s%s",str,reb_version_str);
        for(size_t j=strlen(curvbuf);j<63;j++){
            curvbuf[j] = ' ';
        }
        curvbuf[63] = '\0';
        objects += fread(readbuf,sizeof(char),64,inf);
        if(strcmp(readbuf,curvbuf)!=0){
            *warnings |= REB_INPUT_BINARY_WARNING_VERSION;
        }

        // Read main simulation oject.
        objects += fread(r,sizeof(struct reb_simulation),1,inf);
        int ri_ias15_allocatedN = r->ri_ias15.allocatedN;
        if(reb_reset_function_pointers(r)){
            *warnings |= REB_INPUT_BINARY_WARNING_POINTERS;
        }
        reb_reset_temporary_pointers(r);
        r->allocatedN = r->N;
        r->tree_root = NULL;

        // Read particles
        if (r->N>0){
            r->particles = malloc(sizeof(struct reb_particle)*r->N);
            if (r->particles){
                objects = fread(r->particles,sizeof(struct reb_particle),r->N,inf);
                if (objects==r->N){
                    for (int l=0;l<r->N;l++){
                        r->particles[l].c = NULL;
                        r->particles[l].ap = NULL;
                        r->particles[l].sim = r;
                    }
                }else{
                    *warnings |= REB_INPUT_BINARY_WARNING_PARTICLES;
                }
            }else{
                *warnings |= REB_INPUT_BINARY_WARNING_PARTICLES;
            }
        }
        
        // Read variational config structures
        if (r->var_config_N>0){
            r->var_config = malloc(sizeof(struct reb_variational_configuration)*r->var_config_N);
            if (r->var_config){
                objects = fread(r->var_config,sizeof(struct reb_variational_configuration),r->var_config_N,inf);
                if (objects==r->var_config_N){
                    for (int l=0;l<r->var_config_N;l++){
                        r->var_config[l].sim = r;
                    }
                }else{
                    *warnings |= REB_INPUT_BINARY_WARNING_VARCONFIG;
                }
            }else{
                *warnings |= REB_INPUT_BINARY_WARNING_VARCONFIG;
            }
        }

        // Read temporary arrays for IAS15 (needed for bit-by-bit reproducability)
        if (ri_ias15_allocatedN && !(*warnings & REB_INPUT_BINARY_WARNING_PARTICLES)){
            int N3 = ri_ias15_allocatedN;
            r->ri_ias15.allocatedN = N3;
            r->ri_ias15.at = malloc(sizeof(double)*N3);
            fread(r->ri_ias15.at,sizeof(double),N3,inf);
            r->ri_ias15.x0 = malloc(sizeof(double)*N3);
            fread(r->ri_ias15.x0,sizeof(double),N3,inf);
            r->ri_ias15.v0 = malloc(sizeof(double)*N3);
            fread(r->ri_ias15.v0,sizeof(double),N3,inf);
            r->ri_ias15.a0 = malloc(sizeof(double)*N3);
            fread(r->ri_ias15.a0,sizeof(double),N3,inf);
            r->ri_ias15.csx = malloc(sizeof(double)*N3);
            fread(r->ri_ias15.csx,sizeof(double),N3,inf);
            r->ri_ias15.csv = malloc(sizeof(double)*N3);
            fread(r->ri_ias15.csv,sizeof(double),N3,inf);
            r->ri_ias15.csa0 = malloc(sizeof(double)*N3);
            fread(r->ri_ias15.csa0,sizeof(double),N3,inf);
            reb_read_dp7(&(r->ri_ias15.g)  ,N3,inf);
            reb_read_dp7(&(r->ri_ias15.b)  ,N3,inf);
            reb_read_dp7(&(r->ri_ias15.csb),N3,inf);
            reb_read_dp7(&(r->ri_ias15.e)  ,N3,inf);
            reb_read_dp7(&(r->ri_ias15.br) ,N3,inf);
            reb_read_dp7(&(r->ri_ias15.er) ,N3,inf);
        }
        fclose(inf);
        
        return r;
    }
    *warnings |= REB_INPUT_BINARY_ERROR_NOFILE;
    return NULL;
}


