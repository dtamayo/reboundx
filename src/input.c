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

// Macro to read a single field from a binary file.
#define CASE(typename, valueref) case REBX_BINARY_FIELD_TYPE_##typename: \
{\
if(!fread(valueref, field.size, 1, inf)){\
*warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;\
}\
break;\
}\

#define CASE_MALLOC(typename, valueref) case REBX_BINARY_FIELD_TYPE_##typename: \
{\
valueref = malloc(field.size);\
if(valueref == NULL){\
*warnings |= REBX_INPUT_BINARY_ERROR_NO_MEMORY;\
}\
else{\
if(!fread(valueref, field.size, 1, inf)){\
*warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;\
free(valueref);\
}\
}\
break;\
}\

void rebx_input_skip_binary_field(FILE* inf, long field_size){
    fseek(inf, field_size, SEEK_CUR);
}

static int rebx_load_list(struct rebx_extras* rebx, enum rebx_binary_field_type expected_type, struct rebx_node** ap, FILE* inf, enum rebx_input_binary_messages* warnings);

static struct rebx_param* rebx_read_param(struct rebx_extras* rebx, FILE* inf, enum rebx_input_binary_messages* warnings){
    
    struct rebx_param* param = malloc(sizeof(*param));
    if (param == NULL){
        *warnings |= REBX_INPUT_BINARY_ERROR_NO_MEMORY;
        return NULL;
    }
    param->value = NULL;
    param->name = NULL;
    param->type = REBX_TYPE_NONE;
    
    struct rebx_binary_field field;
    int reading_fields = 1;
    while (reading_fields){
        if (!fread(&field, sizeof(field), 1, inf)){ // means we didn't reach an END field. Corrupt
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            break;
        }
        switch (field.type){
            CASE(PARAM_TYPE,                  &param->type);
            CASE_MALLOC(NAME,                 param->name);
            CASE_MALLOC(PARAM_VALUE,          param->value);
            case REBX_BINARY_FIELD_TYPE_END:
                reading_fields=0;
                break;
            default: // Might have added new fields, saved with new version and loaded with old version
            {
                *warnings |= REBX_INPUT_BINARY_WARNING_FIELD_UNKNOWN;
                break;
            }
        }
    }
    // Check type and name after param has been loaded. Check value later (registered params should have value=NULL)
    if (param->type == REBX_TYPE_NONE || param->name == NULL){
        *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
        rebx_free_param(param);
        return NULL;
    }
    return param;
}

static int rebx_load_param(struct rebx_extras* rebx, struct rebx_node** ap, FILE* inf, enum rebx_input_binary_messages* warnings){
    struct rebx_param* param = rebx_read_param(rebx, inf, warnings);
    
    if(param == NULL){
        return 0;
    }
    
    if(param->value == NULL){
        *warnings |= REBX_INPUT_BINARY_WARNING_PARAM_VALUE_NULL;
        rebx_free_param(param);
        return 0;
    }
    
    if(param->type == REBX_TYPE_FORCE){
        struct rebx_force* force = rebx_get_force(rebx, param->value);
        if (force == NULL){
            *warnings |= REBX_INPUT_BINARY_WARNING_FORCE_PARAM_NOT_LOADED;
            rebx_free_param(param);
            return 0;
        }
        param->value = force;
    }
    int success = rebx_add_param(rebx, ap, param);
    if(!success){
        return 0;
    }
    return 1;
    
}

static int rebx_load_registered_param(struct rebx_extras* rebx, FILE* inf, enum rebx_input_binary_messages* warnings){
    struct rebx_param* param = rebx_read_param(rebx, inf, warnings);
    
    if(param == NULL){
        return 0;
    }
    
    int success = rebx_add_param(rebx, &rebx->registered_params, param);
    if(!success){
        return 0;
    }
    return 1;
}

static char* rebx_load_name(FILE* inf, enum rebx_input_binary_messages* warnings){
    struct rebx_binary_field field;
    if (!fread(&field, sizeof(field), 1, inf)){
        *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
        return NULL;
    }
    if (field.type != REBX_BINARY_FIELD_TYPE_NAME){
        *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
        return NULL;
    }
    char* name = malloc(field.size);
    if (name == NULL){
        *warnings |= REBX_INPUT_BINARY_ERROR_NO_MEMORY;
        return NULL;
    }
    if (!fread(name, field.size, 1, inf)){
        *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
        free(name);
        return NULL;
    }
    return name;
}

static int rebx_load_force_field(struct rebx_extras* rebx, FILE* inf, enum rebx_input_binary_messages* warnings){
    
    // Name of force always comes first so that we can load it
    char* name = rebx_load_name(inf, warnings);
    if(name == NULL){
        return 0;
    }
    struct rebx_force* force = rebx_load_force(rebx, name);
    free(name);
    if(force == NULL){
        *warnings |= REBX_INPUT_BINARY_WARNING_FORCE_NOT_LOADED;
        return 0;
    }
    
    struct rebx_binary_field field;
    int reading_fields = 1;
    while (reading_fields){
        if (!fread(&field, sizeof(field), 1, inf)){
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            return 0;
        }
        switch (field.type){
            case REBX_BINARY_FIELD_TYPE_PARAM_LIST:
            {
                if (!rebx_load_list(rebx, REBX_BINARY_FIELD_TYPE_PARAM, &force->ap, inf, warnings)){
                    return 0;
                }
                break;
            }
            case REBX_BINARY_FIELD_TYPE_END:
            {
                reading_fields=0;
                break;
            }
            default:
            {
                *warnings |= REBX_INPUT_BINARY_WARNING_FIELD_UNKNOWN;
                rebx_input_skip_binary_field(inf, field.size);
                break;
            }
        }
    }
        
    return 1;
}

// Force is already loaded in allocated_forces. Need to get from that list and add to sim
static int rebx_load_additional_force_field(struct rebx_extras* rebx, FILE* inf, enum rebx_input_binary_messages* warnings){
    
    char* name = rebx_load_name(inf, warnings);
    if(name == NULL){
        return 0;
    }
    struct rebx_force* force = rebx_get_force(rebx, name);
    free(name);
    if(force == NULL){
        return 0;
    }
    
    // Just catches END for now. This makes it flexible to addition of fields
    struct rebx_binary_field field;
    int reading_fields = 1;
    while (reading_fields){
        if (!fread(&field, sizeof(field), 1, inf)){
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            return 0;
        }
        switch (field.type){
            case REBX_BINARY_FIELD_TYPE_END:
            {
                reading_fields=0;
                break;
            }
            default:
            {
                *warnings |= REBX_INPUT_BINARY_WARNING_FIELD_UNKNOWN;
                rebx_input_skip_binary_field(inf, field.size);
                break;
            }
        }
    }
    
    int success = rebx_add_force(rebx, force); // add to additional_forces
    return success;
}

static int rebx_load_operator_field(struct rebx_extras* rebx, FILE* inf, enum rebx_input_binary_messages* warnings){
    // Name of force always comes first so that we can load it
    char* name = rebx_load_name(inf, warnings);
    if(name == NULL){
        return 0;
    }
    struct rebx_operator* operator = rebx_load_operator(rebx, name);
    free(name);
    if(operator == NULL){
        *warnings |= REBX_INPUT_BINARY_WARNING_OPERATOR_NOT_LOADED;
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
            case REBX_BINARY_FIELD_TYPE_PARAM_LIST:
            {
                if (!rebx_load_list(rebx, REBX_BINARY_FIELD_TYPE_PARAM, &operator->ap, inf, warnings)){
                    return 0;
                }
                break;
            }
            case REBX_BINARY_FIELD_TYPE_END:
            {
                reading_fields=0;
                break;
            }
            default:
            {
                *warnings |= REBX_INPUT_BINARY_WARNING_FIELD_UNKNOWN;
                rebx_input_skip_binary_field(inf, field.size);
                break;
            }
        }
    }
    
    return 1;
}

static int rebx_load_step_field(struct rebx_extras* rebx, FILE* inf, enum rebx_input_binary_messages* warnings, struct rebx_node** ap){
    char* name = rebx_load_name(inf, warnings);
    if(name == NULL){
        return 0;
    }
    struct rebx_operator* operator = rebx_get_operator(rebx, name);
    free(name);
    if(operator == NULL){
        *warnings |= REBX_INPUT_BINARY_WARNING_OPERATOR_NOT_LOADED;
        return 0;
    }
    
    double dt_fraction = 0.;
    struct rebx_binary_field field;
    int reading_fields = 1;
    while (reading_fields){
        if (!fread(&field, sizeof(field), 1, inf)){
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            break;
        }
        switch (field.type){
            CASE(STEP_DT_FRACTION,                 &dt_fraction);
            case REBX_BINARY_FIELD_TYPE_END:
                reading_fields=0;
                break;
            default:
            {
                *warnings |= REBX_INPUT_BINARY_WARNING_FIELD_UNKNOWN;
                break;
            }
        }
    }
    
    if(dt_fraction == 0){
        return 0;
    }
    
    int success = 0;
    if (ap == &rebx->pre_timestep_modifications) {
        success = rebx_add_operator_step(rebx, operator, dt_fraction, REBX_TIMING_PRE);
    }
    if (ap == &rebx->post_timestep_modifications) {
        success = rebx_add_operator_step(rebx, operator, dt_fraction, REBX_TIMING_POST);
    }
    return success;
}

static int rebx_load_particle(struct rebx_extras* rebx, FILE* inf, enum rebx_input_binary_messages* warnings){
    struct reb_particle* p = NULL;
    struct rebx_binary_field field;
    if (!fread(&field, sizeof(field), 1, inf)){
        *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
        return 0;
    }
    
    if(field.type != REBX_BINARY_FIELD_TYPE_PARTICLE_INDEX){
        *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
        return 0;
    }
    int index;
    if(!fread(&index, field.size, 1, inf)){
        *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
        return 0;
    }
    
    p = &rebx->sim->particles[index]; // checked sim is valid in init_from_binary
    
    int reading_fields = 1;
    while (reading_fields){
        if (!fread(&field, sizeof(field), 1, inf)){
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            return 0;
        }
        switch (field.type){
            case REBX_BINARY_FIELD_TYPE_PARAM_LIST:
            {
                if (!rebx_load_list(rebx, REBX_BINARY_FIELD_TYPE_PARAM, &p->ap, inf, warnings)){
                    return 0;
                }
                break;
            }
            case REBX_BINARY_FIELD_TYPE_END:
            {
                reading_fields=0;
                break;
            }
            default:
            {
                *warnings |= REBX_INPUT_BINARY_WARNING_FIELD_UNKNOWN;
                rebx_input_skip_binary_field(inf, field.size);
                break;
            }
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
            case REBX_BINARY_FIELD_TYPE_END:
            {
                reading_fields=0;
                break;
            }
            case REBX_BINARY_FIELD_TYPE_REGISTERED_PARAMETERS:
            {
                if (!rebx_load_list(rebx, REBX_BINARY_FIELD_TYPE_REGISTERED_PARAM, &rebx->registered_params, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                    rebx_input_skip_binary_field(inf, field.size);
                }
                break;
            }
            case REBX_BINARY_FIELD_TYPE_ALLOCATED_FORCES:
            {
                if (!rebx_load_list(rebx, REBX_BINARY_FIELD_TYPE_FORCE, NULL, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                    rebx_input_skip_binary_field(inf, field.size);
                }
                break;
            }
            case REBX_BINARY_FIELD_TYPE_ALLOCATED_OPERATORS:
            {
                if (!rebx_load_list(rebx, REBX_BINARY_FIELD_TYPE_OPERATOR, NULL, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                    rebx_input_skip_binary_field(inf, field.size);
                }
                break;
            }
            case REBX_BINARY_FIELD_TYPE_ADDITIONAL_FORCES:
            {
                if (!rebx_load_list(rebx, REBX_BINARY_FIELD_TYPE_ADDITIONAL_FORCE, NULL, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                    rebx_input_skip_binary_field(inf, field.size);
                }
                break;
            }
            case REBX_BINARY_FIELD_TYPE_PRE_TIMESTEP_MODIFICATIONS:
            {
                if (!rebx_load_list(rebx, REBX_BINARY_FIELD_TYPE_STEP, &rebx->pre_timestep_modifications, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                    rebx_input_skip_binary_field(inf, field.size);
                }
                break;
            }
            case REBX_BINARY_FIELD_TYPE_POST_TIMESTEP_MODIFICATIONS:
            {
                if (!rebx_load_list(rebx, REBX_BINARY_FIELD_TYPE_STEP, &rebx->post_timestep_modifications, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                    rebx_input_skip_binary_field(inf, field.size);
                }
                break;
            }
            default:
            {
                *warnings |= REBX_INPUT_BINARY_WARNING_FIELD_UNKNOWN;
                rebx_input_skip_binary_field(inf, field.size);
                break;
            }
        }
    }
    return 1;
}

static int rebx_load_snapshot(struct rebx_extras* rebx, FILE* inf, enum rebx_input_binary_messages* warnings){
    struct rebx_binary_field field;
    if (!fread(&field, sizeof(field), 1, inf)){
        *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
        return 0;
    }
    if (field.type != REBX_BINARY_FIELD_TYPE_SNAPSHOT){
        *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
        return 0;
    }

    int reading_fields = 1;
    while (reading_fields){
        if (!fread(&field, sizeof(field), 1, inf)){
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            break;
        }
        switch (field.type){
            case REBX_BINARY_FIELD_TYPE_REBX_STRUCTURE:
            {
                if (!rebx_load_rebx(rebx, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_ERROR_REBX_NOT_LOADED;
                    rebx_input_skip_binary_field(inf, field.size);
                }
                break;
            }
            case REBX_BINARY_FIELD_TYPE_PARTICLES:
            {
                if (!rebx_load_list(rebx, REBX_BINARY_FIELD_TYPE_PARTICLE, NULL, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                    rebx_input_skip_binary_field(inf, field.size);
                }
                break;
            }
            case REBX_BINARY_FIELD_TYPE_END:
            {
                reading_fields=0;
                break;
            }
            default:
            {
                *warnings |= REBX_INPUT_BINARY_WARNING_LIST_UNKNOWN;
                rebx_input_skip_binary_field(inf, field.size);
                break;
            }
        }
    }
    
    return 1;
}

// Only fails (returns 0) if binary is in wrong format
static int rebx_load_list(struct rebx_extras* rebx, enum rebx_binary_field_type expected_type, struct rebx_node** ap, FILE* inf, enum rebx_input_binary_messages* warnings){
    struct rebx_binary_field field;
    int reading_fields = 1;
    while (reading_fields){
        if (!fread(&field, sizeof(field), 1, inf)){
            return 0;
        }
        
        // Check whether we've reached end before checking for expected type
        if (field.type == REBX_BINARY_FIELD_TYPE_END){
            break;
        }
        
        if (field.type != expected_type){
            return 0;
        }
        
        // Only will have fields of expected_type, check function to call
        switch (field.type){
            case REBX_BINARY_FIELD_TYPE_PARAM:
            {
                if(!rebx_load_param(rebx, ap, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_WARNING_PARAM_NOT_LOADED;
                    rebx_input_skip_binary_field(inf, field.size);
                }
                break;
            }
            case REBX_BINARY_FIELD_TYPE_REGISTERED_PARAM:
            {
                if(!rebx_load_registered_param(rebx, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_ERROR_REGISTERED_PARAM_NOT_LOADED;
                    rebx_input_skip_binary_field(inf, field.size);
                }
                break;
            }
            case REBX_BINARY_FIELD_TYPE_FORCE:
            {
                if (!rebx_load_force_field(rebx, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_WARNING_FORCE_NOT_LOADED;
                    rebx_input_skip_binary_field(inf, field.size);
                }
                break;
            }
            case REBX_BINARY_FIELD_TYPE_ADDITIONAL_FORCE:
            {
                if (!rebx_load_additional_force_field(rebx, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_WARNING_ADDITIONAL_FORCE_NOT_LOADED;
                    rebx_input_skip_binary_field(inf, field.size);
                }
                break;
            }
            case REBX_BINARY_FIELD_TYPE_OPERATOR:
            {
                if (!rebx_load_operator_field(rebx, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_WARNING_OPERATOR_NOT_LOADED;
                    rebx_input_skip_binary_field(inf, field.size);
                }
                break;
            }
            case REBX_BINARY_FIELD_TYPE_STEP:
            {
                if (!rebx_load_step_field(rebx, inf, warnings, ap)){
                    *warnings |= REBX_INPUT_BINARY_WARNING_STEP_NOT_LOADED;
                    rebx_input_skip_binary_field(inf, field.size);
                }
                break;
            }
            case REBX_BINARY_FIELD_TYPE_PARTICLE:
            {
                if (!rebx_load_particle(rebx, inf, warnings)){
                    *warnings |= REBX_INPUT_BINARY_WARNING_PARTICLE_PARAMS_NOT_LOADED;
                    rebx_input_skip_binary_field(inf, field.size);
                }
                break;
            }
            default:
            {
                rebx_error(rebx, "REBOUNDx Error. Reached default in rebx_load_list reading binary. Should never reach this case. Means we added a list to rebx and didn't add new case to load_list. Please report bug as Github issue.\n");
                return 0;
            }
        }
    }
    return 1;
}

static void rebx_input_read_header(FILE* inf, enum rebx_input_binary_messages* warnings){
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
}

void rebx_init_extras_from_binary(struct rebx_extras* rebx, const char* const filename, enum rebx_input_binary_messages* warnings){
    if (rebx->sim == NULL){
        rebx_error(rebx, ""); // rebx_error gives meaningful err
        return;
    }
    FILE* inf = fopen(filename,"rb");
    if (!inf){
        *warnings |= REBX_INPUT_BINARY_ERROR_NOFILE;
        return;
    }
    
    rebx_input_read_header(inf, warnings);
    rebx_load_snapshot(rebx, inf, warnings);
    
    fclose(inf);
    return;
}

struct rebx_extras* rebx_create_extras_from_binary(struct reb_simulation* sim, const char* const filename){
    if (sim == NULL){
        fprintf(stderr, "REBOUNDx Error: Simulation pointer passed to rebx_create_extras_from_binary was NULL.\n");
        return NULL;
    }
    enum rebx_input_binary_messages warnings = REBX_INPUT_BINARY_WARNING_NONE;
    // create manually so that default registered parameters not loaded
    struct rebx_extras* rebx = malloc(sizeof(*rebx));
    rebx_initialize(sim, rebx);
    rebx_init_extras_from_binary(rebx, filename, &warnings);
    
    if (warnings & REBX_INPUT_BINARY_ERROR_NOFILE){
        reb_error(sim,"REBOUNDx: Cannot open binary file. Check filename.");
    }
    if (warnings & REBX_INPUT_BINARY_ERROR_CORRUPT){
        reb_error(sim,"REBOUNDx: Binary file is unreadable. Please open an issue on Github mentioning the version of REBOUND and REBOUNDx you are using and include the binary file.");
    }
    if (warnings & REBX_INPUT_BINARY_ERROR_NO_MEMORY){
        reb_error(sim,"REBOUNDx: Ran out of system memory.");
    }
    if (warnings & REBX_INPUT_BINARY_ERROR_REBX_NOT_LOADED){
        reb_error(sim,"REBOUNDx: REBOUNDx structure couldn't be loaded.");
    }
    if (warnings & REBX_INPUT_BINARY_ERROR_REGISTERED_PARAM_NOT_LOADED){
        reb_error(sim,"REBOUNDx: At least one registered parameter was not loaded. This typically indicates the binary is corrupt or was saved with an incompatible version to the current one being used.");
    }
    if (warnings & REBX_INPUT_BINARY_WARNING_PARAM_NOT_LOADED){
        reb_warning(sim,"REBOUNDx: At least one force or operator parameter was not loaded from the binary file. This typically indicates the binary is corrupt or was saved with an incompatible version to the current one being used.");
    }
    if (warnings & REBX_INPUT_BINARY_WARNING_PARTICLE_PARAMS_NOT_LOADED){
        reb_warning(sim,"REBOUNDx: At least one particle's parameters were not loaded from the binary file.");
    }
    if (warnings & REBX_INPUT_BINARY_WARNING_FORCE_NOT_LOADED){
        reb_warning(sim,"REBOUNDx: At least one force was not loaded from the binary file. If binary was created with a newer version of REBOUNDx, a particular force may not be implemented in your current version of REBOUNDx.");
    }
    if (warnings & REBX_INPUT_BINARY_WARNING_OPERATOR_NOT_LOADED){
        reb_warning(sim,"REBOUNDx: At least one operator was not loaded from the binary file. If binary was created with a newer version of REBOUNDx, a particular force may not be implemented in your current version of REBOUNDx.");
    }
    if (warnings & REBX_INPUT_BINARY_WARNING_STEP_NOT_LOADED){
        reb_warning(sim,"REBOUNDx: At least one operator step was not loaded from the binary file.");
    }
    if (warnings & REBX_INPUT_BINARY_WARNING_ADDITIONAL_FORCE_NOT_LOADED){
        reb_warning(sim,"REBOUNDx: At least one force was not added to the simulation. If binary was created with a newer version of REBOUNDx, a particular force may not be implemented in your current version of REBOUNDx.");
    }
    if (warnings & REBX_INPUT_BINARY_WARNING_FIELD_UNKNOWN){
        reb_warning(sim,"REBOUNDx: Unknown field found in binary file. Any unknown fields not loaded.  This can happen if the binary was created with a later version of REBOUNDx than the one used to read it.");
    }
    if (warnings & REBX_INPUT_BINARY_WARNING_LIST_UNKNOWN){
        reb_warning(sim,"REBOUNDx: Unknown list in the REBOUNDx structure wasn't loaded. This can happen if the binary was created with a later version of REBOUNDx than the one used to read it.");
    }
    if (warnings & REBX_INPUT_BINARY_WARNING_PARAM_VALUE_NULL){
        reb_warning(sim,"REBOUNDx: The value of at least one parameter was not loaded. This can happen if a custom structure was added by the user as a parameter. See Parameters.ipynb jupyter notebook example.");
    }
    if (warnings & REBX_INPUT_BINARY_WARNING_VERSION){
        reb_warning(sim,"REBOUNDx: Binary file was saved with a different version of REBOUNDx. Binary format might have changed. Check that effects and parameters are loaded as expected.");
    }
    if (warnings & REBX_INPUT_BINARY_WARNING_FORCE_PARAM_NOT_LOADED){
        reb_warning(sim,"REBOUNDx: A force parameter failed to load from the list of REBOUNDx implemented forces. Custom forces can't be saved to a REBOUNDx binary, and function points must be reset when a simulation is reloaded.");
    }
    return rebx;
}

FILE* rebx_input_inspect_binary(const char* const filename, enum rebx_input_binary_messages* warnings){
    FILE* inf = fopen(filename,"rb");
    if (!inf){
        *warnings |= REBX_INPUT_BINARY_ERROR_NOFILE;
        return NULL;
    }
    
    rebx_input_read_header(inf, warnings);
    return inf;
}

struct rebx_binary_field rebx_input_read_binary_field(FILE* inf){
    struct rebx_binary_field field;
    if (fread(&field, sizeof(field), 1, inf)){
        return field;
    }
    else{
        struct rebx_binary_field empty = {0};
        return empty;
    }
}
