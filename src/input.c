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

static int rebx_load_param(struct rebx_extras* rebx, struct rebx_node** ap, FILE* inf, enum rebx_input_binary_messages* warnings){
    struct rebx_param* param = rebx_malloc(rebx, sizeof(*param));
    if (param == NULL){
        *warnings |= REBX_INPUT_BINARY_ERROR_NO_MEMORY;
        return 0;
    }
    param->value = NULL;
    param->name = NULL;
    param->type = REBX_TYPE_NONE;
    
    struct rebx_binary_field field;
    int reading_fields = 1;
    while (reading_fields){
        if (!fread(&field, sizeof(field), 1, inf)){ // means we didn't reach an END field. Corrupt
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            rebx_free_param(param);
            return 0;
        }
        switch (field.type){
            case REBX_BINARY_FIELD_TYPE_PARAM_TYPE:
            {
                fread(&param->type, field.size, 1, inf);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_VALUE:
            {
                param->value = rebx_malloc(rebx, field.size);
                fread(param->value, field.size, 1, inf);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_NAME:
            {
                param->name = malloc(field.size);
                fread(param->name, field.size, 1, inf);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_END:
            {
                reading_fields=0;
                break;
            }
            default: // Might have added new fields, saved with new version and loaded with old version
            {
                *warnings |= REBX_INPUT_BINARY_WARNING_FIELD_UNKOWN; // not necessarily fatal
                fseek(inf, field.size, SEEK_CUR);
                break;
            }
        }
    }
    
    // Necessary for valid param. Don't check value since registered_params will have NULL
    if (param->type == REBX_TYPE_NONE){ // Didn't include type
        *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
        rebx_free_param(param);
        return 0;
    }
    if (param->name == NULL){ // All params must include name
        *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
        rebx_free_param(param);
        return 0;
    }
    
    rebx_add_param(rebx, ap, param);
    return 1;
}

// NEED TO LOAD_FORCE so that function pointers are set
static int rebx_load_force_field(struct rebx_extras* rebx, FILE* inf, enum rebx_input_binary_messages* warnings){
    struct rebx_force* force = rebx_create_force(rebx, NULL);
    struct rebx_binary_field field;
    int reading_fields = 1;
    while (reading_fields){
        if (!fread(&field, sizeof(field), 1, inf)){
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            rebx_remove_force(rebx, force);
            return 0;
        }
        switch (field.type){
            case REBX_BINARY_FIELD_TYPE_NAME:
            {
                force->name = rebx_malloc(rebx, field.size);
                fread(force->name, field.size, 1, inf);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_FORCE_TYPE:
            {
                fread(&force->force_type, field.size, 1, inf);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_PARAM:
            {
                long field_start = ftell(inf);
                if (!rebx_load_param(rebx, &force->ap, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_WARNING_PARAM_NOT_LOADED;
                    fseek(inf, field_start + field.size, SEEK_SET);
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

// Force is already loaded in allocated_forces. Need to get from that list and add to sim
static int rebx_load_additional_force_field(struct rebx_extras* rebx, FILE* inf, enum rebx_input_binary_messages* warnings){
    struct rebx_force* force = NULL;
    struct rebx_binary_field field;
    int reading_fields = 1;
    while (reading_fields){
        if (!fread(&field, sizeof(field), 1, inf)){
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            return 0;
        }
        switch (field.type){
            case REBX_BINARY_FIELD_TYPE_NAME:
            {
                char* name = malloc(field.size);
                fread(name, field.size, 1, inf);
                force = rebx_get_force(rebx, name);
                free(name);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_END:
            {
                reading_fields=0;
                break;
            }
            default:
                *warnings |= REBX_INPUT_BINARY_WARNING_FIELD_UNKOWN;
                fseek(inf, field.size, SEEK_CUR);
                break;
        }
    }
    if(force != NULL){
        int success = rebx_add_force(rebx, force);
        return success;
    }
    else{
        *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
        return 0;
    }
}

static int rebx_load_operator_field(struct rebx_extras* rebx, FILE* inf, enum rebx_input_binary_messages* warnings){
    struct rebx_operator* operator = rebx_create_operator(rebx, NULL);
    struct rebx_binary_field field;
    int reading_fields = 1;
    while (reading_fields){
        if (!fread(&field, sizeof(field), 1, inf)){
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            rebx_remove_operator(rebx, operator);
            return 0;
        }
        switch (field.type){
            case REBX_BINARY_FIELD_TYPE_NAME:
            {
                operator->name = rebx_malloc(rebx, field.size);
                fread(operator->name, field.size, 1, inf);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_OPERATOR_TYPE:
            {
                fread(&operator->operator_type, field.size, 1, inf);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_PARAM:
            {
                long field_start = ftell(inf);
                if (!rebx_load_param(rebx, &operator->ap, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_WARNING_PARAM_NOT_LOADED;
                    fseek(inf, field_start + field.size, SEEK_SET);
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
                fseek(inf, field.size, SEEK_CUR);
                break;
        }
    }
    return 1;
}

// NEED TO LOAD_OPERATOR TO MAKE SURE FUNCTION POINTERS ARE SET
static int rebx_load_step_field(struct rebx_extras* rebx, FILE* inf, enum rebx_input_binary_messages* warnings, enum rebx_timing timing){
    char* stepname = NULL;
    char* operatorname = NULL;
    double dt_fraction = 0.;
    struct rebx_operator* operator = NULL;
    
    struct rebx_binary_field field;
    int reading_fields = 1;
    while (reading_fields){
        if (!fread(&field, sizeof(field), 1, inf)){
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            return 0;
        }
        switch (field.type){
            case REBX_BINARY_FIELD_TYPE_NAME:
            {
                stepname = malloc(field.size);
                fread(stepname, field.size, 1, inf);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_OPERATOR_NAME:
            {
                operatorname = malloc(field.size);
                fread(operatorname, field.size, 1, inf);
                operator = rebx_get_operator(rebx, operatorname);
                free(operatorname);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_DT_FRACTION:
            {
                fread(&dt_fraction, field.size, 1, inf);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_END:
            {
                reading_fields=0;
                break;
            }
            default:
                *warnings |= REBX_INPUT_BINARY_WARNING_FIELD_UNKOWN;
                fseek(inf, field.size, SEEK_CUR);
                break;
        }
    }
    if(operator != NULL){
        int success = rebx_add_force(rebx, force);
        return success;
    }
    else{
        *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
        return 0;
    }
}

static int rebx_load_particle(struct rebx_extras* rebx, FILE* inf, enum rebx_input_binary_messages* warnings){
    struct reb_particle* p = NULL;
    struct rebx_binary_field field;
    int reading_fields = 1;
    while (reading_fields){
        if (!fread(&field, sizeof(field), 1, inf)){
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            return 0;
        }
        switch (field.type){
            case REBX_BINARY_FIELD_TYPE_PARTICLE_INDEX:
            {
                int index;
                fread(&index, field.size, 1, inf);
                p = &rebx->sim->particles[index];
                break;
            }
            case REBX_BINARY_FIELD_TYPE_PARAM:
            {
                if(!p){ // PARTICLE_INDEX should always be first in binary
                    *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                    return 0;
                }
                long field_start = ftell(inf);
                if (!rebx_load_param(rebx, &p->ap, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_WARNING_PARAM_NOT_LOADED;
                    fseek(inf, field_start + field.size, SEEK_SET);
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
                fseek(inf, field.size, SEEK_CUR);
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
