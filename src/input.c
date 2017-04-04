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

#define CASE(typename, value) case REBX_BINARY_FIELD_TYPE_##typename: \
    {\
    fread(value, field.size,1,inf);\
    }\
    break;

#define CASE_MALLOC(typename, valueref) case REBX_BINARY_FIELD_TYPE_##typename: \
    {\
    valueref = malloc(field.size);\
    fread(valueref, field.size,1,inf);\
    }\
    break;

static int rebx_load_param(struct rebx_extras* rebx, void* const object, FILE* inf, enum rebx_input_binary_messages* warnings){
    struct rebx_param* param = rebx_create_param();
    if(!param){
        return 0;
    }
    struct rebx_binary_field field;
    int reading_fields = 1;
    while (reading_fields){
        if (!fread(&field, sizeof(field), 1, inf)){
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            break;
        }
        switch (field.type){
            CASE(PARAM_TYPE,        &param->param_type);
            CASE(PYTHON_TYPE,       &param->python_type);
            CASE(NDIM,              &param->ndim);
            CASE_MALLOC(CONTENTS,   param->contents);
            case REBX_BINARY_FIELD_TYPE_NAME:
            {
                param->name = malloc(field.size);
                fread(param->name, field.size,1,inf);
                param->hash = reb_hash(param->name);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_SHAPE:
            {
                param->size = 1;
                if (param->ndim > 0){
                    param->shape = malloc(field.size);
                    param->strides = malloc(field.size);
                    fread(param->shape, field.size,1,inf);
                    for(int i=param->ndim-1;i>=0;i--){
                        param->strides[i] = param->size;
                        param->size *= param->shape[i];
                    }
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
    param = rebx_attach_param_node(object, param);
    return 1;
}

static int rebx_load_effect(struct rebx_extras* rebx, FILE* inf, enum rebx_input_binary_messages* warnings){
    char* name = NULL;
    struct rebx_effect* effect = NULL;
    struct rebx_binary_field field;
    int reading_fields = 1;
    while (reading_fields){
        if (!fread(&field, sizeof(field), 1, inf)){
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            break;
        }
        switch (field.type){
            case REBX_BINARY_FIELD_TYPE_NAME:
            {
                name = malloc(field.size);
                fread(name, field.size, 1, inf);
                effect = rebx_add(rebx, name);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_FORCE_AS_OPERATOR:
            {
                if(!effect){
                    return 0;
                }
                fread(&effect->force_as_operator, sizeof(effect->force_as_operator), 1, inf);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_OPERATOR_ORDER:
            {
                if(!effect){
                    return 0;
                }
                fread(&effect->operator_order, sizeof(effect->operator_order), 1, inf);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_PARAM:
            {
                if(!effect){
                    return 0;
                }
                long field_start = ftell(inf);
                if (!rebx_load_param(rebx, effect, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_WARNING_PARAM_NOT_LOADED;
                    fseek(inf,field_start + field.size,SEEK_SET);
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
    return 1;
}

static int rebx_load_particle(struct rebx_extras* rebx, FILE* inf, enum rebx_input_binary_messages* warnings){
    struct reb_particle* p = NULL;
    struct rebx_binary_field field;
    int reading_fields = 1;
    while (reading_fields){
        if (!fread(&field, sizeof(field), 1, inf)){
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            break;
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
                    return 0;
                }
                long field_start = ftell(inf);
                if (!rebx_load_param(rebx, p, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_WARNING_PARAM_NOT_LOADED;
                    fseek(inf,field_start + field.size,SEEK_SET);
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

static int rebx_load_rebx(struct rebx_extras* rebx, FILE* inf, enum rebx_input_binary_messages* warnings){
    struct rebx_binary_field field;
    int reading_fields = 1;
    while (reading_fields){
        if (!fread(&field, sizeof(field), 1, inf)){
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            break;
        }
        switch (field.type){
            case REBX_BINARY_FIELD_TYPE_REBX_INTEGRATOR:
            {
                fread(&rebx->integrator, sizeof(rebx->integrator), 1, inf);
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
    const char zero = '\0';
    char readbuf[65], curvbuf[65];
    sprintf(curvbuf,"%s%s",str,rebx_version_str);
    memcpy(curvbuf+strlen(curvbuf)+1,rebx_githash_str,sizeof(char)*(62-strlen(curvbuf)));
    curvbuf[63] = zero;
    
    objects += fread(readbuf,sizeof(*str),64,inf);
    // Note: following compares version, but ignores githash.
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
            case REBX_BINARY_FIELD_TYPE_REBX_STRUCTURE:
            {
                rebx_load_rebx(rebx, inf, warnings);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_EFFECT:
            {
                long field_start = ftell(inf);
                if (!rebx_load_effect(rebx, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_WARNING_EFFECT_NOT_LOADED;
                    fseek(inf,field_start + field.size,SEEK_SET);
                }
                break;
            }
            case REBX_BINARY_FIELD_TYPE_PARTICLE:
            {
                long field_start = ftell(inf);
                if (!rebx_load_particle(rebx, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_WARNING_PARTICLE_NOT_LOADED;
                    fseek(inf,field_start + field.size,SEEK_SET);
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

    fclose(inf);
    return;
}

struct rebx_extras* rebx_create_extras_from_binary(struct reb_simulation* sim, const char* const filename){
    enum rebx_input_binary_messages warnings = REBX_INPUT_BINARY_WARNING_NONE;
    struct rebx_extras* rebx = rebx_init(sim);
    rebx_create_extras_from_binary_with_messages(rebx, filename, &warnings);
    
    if (warnings & REBX_INPUT_BINARY_ERROR_NOFILE){
        reb_error(sim,"REBOUNDx: Cannot open binary file. Check filename.");
        rebx_free(rebx);
        rebx = NULL;
    }
    if (warnings & REBX_INPUT_BINARY_ERROR_CORRUPT){
        reb_error(sim,"REBOUNDx: Binary file is unreadable. Please open an issue on Github mentioning the version of REBOUND and REBOUNDx you are using and including the binary file.");
        rebx_free(rebx);
        rebx = NULL;
    }
    if (warnings & REBX_INPUT_BINARY_WARNING_VERSION){
        reb_warning(sim,"REBOUNDx: Binary file was saved with a different version of REBOUNDx. Binary format might have changed. Check that effects and parameters are loaded as expected.");
    }
    if (warnings & REBX_INPUT_BINARY_WARNING_PARAM_NOT_LOADED){
        reb_warning(sim,"REBOUNDx: At least one parameter was not loaded from the binary file.");
    }
    if (warnings & REBX_INPUT_BINARY_WARNING_PARTICLE_NOT_LOADED){
        reb_warning(sim,"REBOUNDx: At least one particle's parameters were not loaded from the binary file.");
    }
    if (warnings & REBX_INPUT_BINARY_WARNING_EFFECT_NOT_LOADED){
        reb_warning(sim,"REBOUNDx: At least one effect was not loaded from the binary file. If binary was created with newer version of REBOUNDx, a particular effect may not be implemented in your current version of REBOUNDx.");
    }
    if (warnings & REBX_INPUT_BINARY_WARNING_FIELD_UNKOWN){
        reb_warning(sim,"REBOUNDx: Unknown field found in binary file. Any unknown fields not loaded.  This can happen if the binary was created with a later version of REBOUNDx than the one used to read it.");
    }
    
    return rebx;
}
