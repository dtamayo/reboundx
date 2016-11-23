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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "rebound.h"
#include "reboundx.h"
#include "core.h"

#define CASE(typename, value) case REB_BINARY_FIELD_TYPE_##typename: \
    {\
    fread(value, field.size,1,inf);\
    }\
    break;

int rebx_load_param(struct rebx_extras* rebx, void* const object, FILE* inf, enum rebx_input_binary_messages* warnings){
    size_t namelength = 0;
    char* name = NULL;
    enum rebx_param_type param_type = INT_MIN;
    int ndim = INT_MIN;
    int* shape = NULL;

    struct rebx_binary_field field;
    int reading_fields = 1;
    while (reading_fields){
        fread(&field, sizeof(field), 1, inf);
        switch (field.type){
            case REBX_BINARY_FIELD_TYPE_PARAM_TYPE:
            {
                fread(&param_type, sizeof(param_type), 1, inf);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_NDIM:
            {
                fread(&ndim, sizeof(ndim), 1, inf);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_NAMELENGTH:
            {
                fread(&namelength, sizeof(namelength), 1, inf);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_NAME:
            {
                if (namelength == 0){
                    fseek(inf,field.size,SEEK_CUR);
                    break;
                }
                name = malloc(namelength);
                fread(name, sizeof(*name), namelength, inf);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_SHAPE:
            {
                if(ndim<0){
                    fseek(inf,field.size,SEEK_CUR);
                    break;
                }
                shape = malloc(ndim*sizeof(*shape));
                fread(shape, sizeof(*shape), ndim, inf);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_CONTENTS:
            {
                if (!name || param_type == INT_MIN || ndim == INT_MIN){ // OK for shape to be NULL if ndim = 0
                    *warnings |= REBX_INPUT_BINARY_WARNING_PARAM_NOT_LOADED;
                    fseek(inf,field.size,SEEK_CUR);
                }
                else{
                    struct rebx_param* param = rebx_add_param_node(object, name, param_type, ndim, shape);
                    fread(param->contents, rebx_sizeof(param->param_type), param->size, inf);
                }
                break;
            }
            case REBX_BINARY_FIELD_TYPE_END:
            {
                reading_fields=0;
                break;
            }
            default:
                *warnings |= REBX_INPUT_BINARY_WARNING_FIELD_UNKOWN;
                fseek(inf,field.size,SEEK_CUR);
                break;
        }
    }
    free(name);
    free(shape);
    return 1;
}

static int rebx_load_effect(struct rebx_extras* rebx, FILE* inf, enum rebx_input_binary_messages* warnings){
    size_t namelength = 0;
    char* name = NULL;
    
    struct rebx_effect* effect = NULL;
    struct rebx_binary_field field;
    int reading_fields = 1;
    while (reading_fields){
        fread(&field, sizeof(field), 1, inf);
        switch (field.type){
            case REBX_BINARY_FIELD_TYPE_NAMELENGTH:
            {
                fread(&namelength, sizeof(namelength), 1, inf);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_NAME:
            {
                if (namelength == 0){
                    fseek(inf,field.size,SEEK_CUR);
                    break;
                }
                name = malloc(namelength);
                fread(name, sizeof(*name), namelength, inf);
                effect = rebx_add(rebx, name);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_PARAM:
            {
                if(!effect){
                    *warnings |= REBX_INPUT_BINARY_WARNING_PARAM_NOT_LOADED;
                    fseek(inf,field.size,SEEK_CUR);
                }
                rebx_load_param(rebx, effect, inf, warnings);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_END:
            {
                reading_fields=0;
                break;
            }
            default:
                *warnings |= REBX_INPUT_BINARY_WARNING_FIELD_UNKOWN;
                fseek(inf,field.size,SEEK_CUR);
                break;
        }
    }
    free(name);
    return 1;
}

static int rebx_load_particle(struct rebx_extras* rebx, FILE* inf, enum rebx_input_binary_messages* warnings){
    struct reb_particle* p = NULL;
    struct rebx_binary_field field;
    int reading_fields = 1;
    while (reading_fields){
        if (!fread(&field, sizeof(field), 1, inf)){
            return 0;
        }
        switch (field.type){
            case REBX_BINARY_FIELD_TYPE_PARTICLE_INDEX:
            {
                int index;
                fread(&index, sizeof(index), 1, inf);
                p = &rebx->sim->particles[index];
                break;
            }
            case REBX_BINARY_FIELD_TYPE_PARAM:
            {
                if(!p){
                    *warnings |= REBX_INPUT_BINARY_WARNING_PARAM_NOT_LOADED;
                    fseek(inf,field.size,SEEK_CUR);
                }
                else if (!rebx_load_param(rebx, p, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_WARNING_PARAM_NOT_LOADED;
                    fseek(inf,field.size,SEEK_CUR);
                }
                break;
            }
            case REBX_BINARY_FIELD_TYPE_END:
            {
                reading_fields=0;
                break;
            }
            default:
                *warnings |= REBX_INPUT_BINARY_WARNING_FIELD_UNKOWN;
                fseek(inf,field.size,SEEK_CUR);
                break;
        }
    }
    return 1;
}


void rebx_create_extras_from_binary_with_messages(struct rebx_extras* rebx, const char* const filename, enum rebx_input_binary_messages* warnings){
    FILE* inf = fopen(filename,"rb");
    if (!inf){
        *warnings |= REBX_INPUT_BINARY_ERROR_NOFILE;
        return;
    }
    
    long objects = 0;
    // Input header.
    const char str[] = "REBOUNDx Binary File. Version: ";
    char readbuf[65], curvbuf[65];
    sprintf(curvbuf,"%s%s",str,rebx_version_str);
    for(size_t j=strlen(curvbuf);j<63;j++){
        curvbuf[j] = ' ';
    }
    curvbuf[63] = '\0';
    objects += fread(readbuf,sizeof(*str),64,inf);
    if(strcmp(readbuf,curvbuf)!=0){
        *warnings |= REBX_INPUT_BINARY_WARNING_VERSION;
    }
    
    struct rebx_binary_field field;
    int reading_fields = 1;
    while (reading_fields){
        if (!fread(&field, sizeof(field), 1, inf)){
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            break;
        }
        switch (field.type){
            case REBX_BINARY_FIELD_TYPE_EFFECT:
            {
                reading_fields = rebx_load_effect(rebx, inf, warnings);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_PARTICLE:
            {
                reading_fields = rebx_load_particle(rebx, inf, warnings);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_END:
            {
                reading_fields=0;
                break;
            }
            default:
                *warnings |= REBX_INPUT_BINARY_WARNING_FIELD_UNKOWN;
                fseek(inf,field.size,SEEK_CUR);
                break;
        }
    }

    fclose(inf);
    return;
}

struct rebx_extras* rebx_create_extras_from_binary(struct reb_simulation* sim, const char* const filename){
    enum rebx_input_binary_messages warnings = REBX_INPUT_BINARY_WARNING_NONE;
    struct rebx_extras* rebx = rebx_init(sim);
    rebx_create_extras_from_binary_with_messages(rebx, filename, &warnings);
    
    if (warnings & REBX_INPUT_BINARY_WARNING_VERSION){
        reb_warning(sim,"REBOUNDx: Binary file was saved with a different version of REBOUNDx. Binary format might have changed and corrupted the loading. Check that effects and parameters are loaded as expected.");
    }
    if (warnings & REBX_INPUT_BINARY_WARNING_FIELD_UNKOWN){
        reb_warning(sim,"REBOUNDx: Unknown field found in binary file. Any unknown fields not loaded.  This can happen if the binary was created with a later version of REBOUNDx than the one used to read it.");
    }
    if (warnings & REBX_INPUT_BINARY_ERROR_NOFILE){
        reb_error(sim,"REBOUNDx: Cannot read binary file. Check filename and file contents.");
        rebx_free(rebx);
        rebx = NULL;
    }
    return rebx;
}
