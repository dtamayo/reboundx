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

/* Macro to write a single field to a binary file.
Each FIELD (e.g. the name of a param, or its type) has format in binary of :
- rebx_binary_field struct
- The value of the field
The field struct holds the type for the field value below so you know how to read it,
and the size is for the struct AND the value of the field which follows, so you know how to skip it.
Macro keeps track of size and writes it in at the end.
*/
 
#define REBX_WRITE_FIELD(typename, valueptr, typesize) {\
            long pos_field_rewrite = ftell(of);\
            struct rebx_binary_field field = {.type = REBX_BINARY_FIELD_TYPE_##typename, .size = 0};\
            fwrite(&field,sizeof(field),1,of);\
            long pos_start_field = ftell(of);\
            fwrite(valueptr,typesize,1,of);\
            long pos_end_field = ftell(of);\
            field.size = pos_end_field - pos_start_field;\
            fseek(of, pos_field_rewrite, SEEK_SET);\
            fwrite(&field, sizeof(field), 1, of);\
            fseek(of, 0, SEEK_END);\
        }

/*
 Each PARAM has format in binary of :
 - rebx_binary_field struct for param
 - param_type field
 - name field
 - value field
 - END field
 The binary_field struct type field says what follows is a param.
 Then there are binary_field/value pairs for each field in the param struct and an end field.
 The size in the top level binary_field for the param is for this whole set of entries so we can skip.
 Macro keeps track of size and writes it in at the end.
 */

// forward declaration
static void rebx_write_list(struct rebx_extras* rebx, enum rebx_binary_field_type list_type, struct rebx_node* list, FILE* of);

static void rebx_write_param(struct rebx_extras* rebx, struct rebx_param* param, FILE* of){
    long pos_param_rewrite = ftell(of);
    struct rebx_binary_field param_field = {.type = REBX_BINARY_FIELD_TYPE_PARAM, .size=0};
    fwrite(&param_field, sizeof(param_field), 1, of);
    long pos_start_param = ftell(of);
    
    REBX_WRITE_FIELD(PARAM_TYPE, &param->type,     sizeof(param->type));
    REBX_WRITE_FIELD(NAME,       param->name,      strlen(param->name) + 1); // +1 for \0 at end
    REBX_WRITE_FIELD(VALUE,      param->value,     rebx_sizeof(rebx, param->type));
    REBX_WRITE_FIELD(END,        NULL,             0);

    // Go back and write size of entire parameter and reset file position to end of file
    long pos_end_param = ftell(of);
    param_field.size = pos_end_param - pos_start_param;
    fseek(of, pos_param_rewrite, SEEK_SET);
    fwrite(&param_field, sizeof(param_field), 1, of);
    fseek(of, 0, SEEK_END);
}

static void rebx_write_registered_param(struct rebx_extras* rebx, struct rebx_param* param, FILE* of){
    long pos_param_rewrite = ftell(of);
    struct rebx_binary_field param_field = {.type = REBX_BINARY_FIELD_TYPE_REGISTERED_PARAM, .size=0};
    fwrite(&param_field, sizeof(param_field), 1, of);
    long pos_start_param = ftell(of);
    
    REBX_WRITE_FIELD(PARAM_TYPE, &param->type,     sizeof(param->type));
    REBX_WRITE_FIELD(NAME,       param->name,      strlen(param->name) + 1); // +1 for \0 at end
    REBX_WRITE_FIELD(END,        NULL,             0);
    
    // Go back and write size of entire parameter and reset file position to end of file
    long pos_end_param = ftell(of);
    param_field.size = pos_end_param - pos_start_param;
    fseek(of, pos_param_rewrite, SEEK_SET);
    fwrite(&param_field, sizeof(param_field), 1, of);
    fseek(of, 0, SEEK_END);
}

static void rebx_write_force(struct rebx_extras* rebx, struct rebx_force* force, FILE* of){
    long pos_effect_rewrite = ftell(of);
    struct rebx_binary_field force_field = {.type = REBX_BINARY_FIELD_TYPE_FORCE, .size=0};
    fwrite(&force_field, sizeof(force_field), 1, of);
    long pos_start_force = ftell(of);
    
    // must write name first so that force can be loaded on read
    REBX_WRITE_FIELD(NAME, force->name, strlen(force->name) + 1); // +1 for \0 at end
    
    REBX_WRITE_FIELD(PARAM_LIST, NULL,   0);
    rebx_write_list(rebx, REBX_BINARY_FIELD_TYPE_PARAM, force->ap, of); // Write all parameters
    REBX_WRITE_FIELD(END,           NULL,               0); // end of param list
    
    REBX_WRITE_FIELD(END,           NULL,               0); // end of force

    // Go back and write size of entire parameter and reset file position to end of file
    long pos_end_force = ftell(of);
    force_field.size = pos_end_force - pos_start_force;
    fseek(of, pos_effect_rewrite, SEEK_SET);
    fwrite(&force_field, sizeof(force_field), 1, of);
    fseek(of, 0, SEEK_END);
}

// Same as force, but only holds the name for later loading, rather than the whole parameter list
static void rebx_write_additional_force(struct rebx_extras* rebx, struct rebx_force* force, FILE* of){
    long pos_effect_rewrite = ftell(of);
    struct rebx_binary_field force_field = {.type = REBX_BINARY_FIELD_TYPE_ADDITIONAL_FORCE, .size=0};
    fwrite(&force_field, sizeof(force_field), 1, of);
    long pos_start_effect = ftell(of);
    
    REBX_WRITE_FIELD(NAME,          force->name,        strlen(force->name) + 1); // +1 for \0 at end
    REBX_WRITE_FIELD(END,           NULL,               0);

    // Go back and write size of entire parameter and reset file position to end of file
    long pos_end_effect = ftell(of);
    force_field.size = pos_end_effect - pos_start_effect;
    fseek(of, pos_effect_rewrite, SEEK_SET);
    fwrite(&force_field, sizeof(force_field), 1, of);
    fseek(of, 0, SEEK_END);
}

static void rebx_write_operator(struct rebx_extras* rebx, struct rebx_operator* operator, FILE* of){
    long pos_effect_rewrite = ftell(of);
    struct rebx_binary_field operator_field = {.type = REBX_BINARY_FIELD_TYPE_OPERATOR, .size=0};
    fwrite(&operator_field, sizeof(operator_field), 1, of);
    long pos_start_effect = ftell(of);
    
    REBX_WRITE_FIELD(NAME,          operator->name,             strlen(operator->name) + 1); // +1 for \0 at end
    
    REBX_WRITE_FIELD(PARAM_LIST, NULL,   0);
    rebx_write_list(rebx, REBX_BINARY_FIELD_TYPE_PARAM, operator->ap, of); // Write all parameters
    REBX_WRITE_FIELD(END,           NULL,               0); // end of param list
    
    REBX_WRITE_FIELD(END,           NULL,               0); // end of force

    // Go back and write size of entire parameter and reset file position to end of file
    long pos_end_effect = ftell(of);
    operator_field.size = pos_end_effect - pos_start_effect;
    fseek(of, pos_effect_rewrite, SEEK_SET);
    fwrite(&operator_field, sizeof(operator_field), 1, of);
    fseek(of, 0, SEEK_END);
}

static void rebx_write_step(struct rebx_extras* rebx, struct rebx_step* step, FILE* of){
    long pos_effect_rewrite = ftell(of);
    struct rebx_binary_field step_field = {.type = REBX_BINARY_FIELD_TYPE_STEP, .size=0};
    fwrite(&step_field, sizeof(step_field), 1, of);
    long pos_start_effect = ftell(of);
    
    // Need operator name to link it back when reading it back in
    REBX_WRITE_FIELD(NAME, step->operator->name,   strlen(step->operator->name) + 1);
    REBX_WRITE_FIELD(DT_FRACTION,   &step->dt_fraction,     sizeof(step->dt_fraction));
    REBX_WRITE_FIELD(END,           NULL,                   0);
    
    // Go back and write size of entire parameter and reset file position to end of file
    long pos_end_effect = ftell(of);
    step_field.size = pos_end_effect - pos_start_effect;
    fseek(of, pos_effect_rewrite, SEEK_SET);
    fwrite(&step_field, sizeof(step_field), 1, of);
    fseek(of, 0, SEEK_END);
}

static void rebx_write_particle(struct rebx_extras* rebx, struct reb_particle* particle, int index, FILE* of){
    long pos_particle_rewrite = ftell(of);
    struct rebx_binary_field particle_field = {.type = REBX_BINARY_FIELD_TYPE_PARTICLE, .size=0};
    fwrite(&particle_field, sizeof(particle_field), 1, of);
    long pos_start_particle = ftell(of);
    
    REBX_WRITE_FIELD(PARTICLE_INDEX,    &index, sizeof(index));
    
    REBX_WRITE_FIELD(PARAM_LIST, NULL,   0);
    rebx_write_list(rebx, REBX_BINARY_FIELD_TYPE_PARAM, particle->ap, of); // Write all parameters
    REBX_WRITE_FIELD(END,           NULL,               0); // end of param list
    
    // Write end marker for particle
    REBX_WRITE_FIELD(END,               NULL,   0);
    
    // Go back and write size of entire parameter and reset file position to end of file
    long pos_end_particle = ftell(of);
    particle_field.size = pos_end_particle - pos_start_particle;
    fseek(of, pos_particle_rewrite, SEEK_SET);
    fwrite(&particle_field, sizeof(particle_field), 1, of);
    fseek(of, 0, SEEK_END);
}

static void rebx_write_rebx(struct rebx_extras* rebx, FILE* of){
    long pos_rebx_rewrite = ftell(of);
    struct rebx_binary_field rebx_field = {.type = REBX_BINARY_FIELD_TYPE_REBX_STRUCTURE, .size=0};
    fwrite(&rebx_field, sizeof(rebx_field), 1, of);
    long pos_start_rebx = ftell(of);
    
    REBX_WRITE_FIELD(END,               NULL,               0);
    
    // Go back and write size of entire parameter and reset file position to end of file
    long pos_end_rebx = ftell(of);
    rebx_field.size = pos_end_rebx - pos_start_rebx;
    fseek(of, pos_rebx_rewrite, SEEK_SET);
    fwrite(&rebx_field, sizeof(rebx_field), 1, of);
    fseek(of, 0, SEEK_END);
}

static void rebx_write_list(struct rebx_extras* rebx, enum rebx_binary_field_type list_type, struct rebx_node* list, FILE* of){
    
    int N = rebx_len(list);
    while (N > 0){
        struct rebx_node* current = list;
        for(int i=0;i<N-1;i++){
            current=current->next;
        }
        switch(list_type){
            case REBX_BINARY_FIELD_TYPE_REGISTERED_PARAM:
            {
                rebx_write_registered_param(rebx, current->object, of);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_FORCE:
            {
                rebx_write_force(rebx, current->object, of);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_ADDITIONAL_FORCE:
            {
                rebx_write_additional_force(rebx, current->object, of);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_OPERATOR:
            {
                rebx_write_operator(rebx, current->object, of);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_PARAM:
            {
                rebx_write_param(rebx, current->object, of);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_STEP:
            {
                rebx_write_step(rebx, current->object, of);
                break;
            }
        }
        N--;
    }
}

void rebx_output_binary(struct rebx_extras* rebx, char* filename){
    struct reb_simulation* sim = rebx->sim;

    FILE* of = fopen(filename,"wb"); 
    if (of==NULL){
        reb_error(sim, "REBOUNDx error: Can not open file passed to rebx_output_binary.");
    }

    // Output header.
    const char str[] = "REBOUNDx Binary File. Version: ";
    char zero = '\0';
    size_t lenheader = strlen(str)+strlen(rebx_version_str);
    fwrite(str,sizeof(char),strlen(str),of);
    fwrite(rebx_version_str,sizeof(char), strlen(rebx_version_str),of);
    fwrite(&zero,sizeof(char),1,of);
    fwrite(rebx_githash_str,sizeof(char),62-lenheader,of);
    fwrite(&zero,sizeof(char),1,of);

    rebx_write_rebx(rebx, of);
    
    // Write registered parameters
    REBX_WRITE_FIELD(REGISTERED_PARAMETERS, NULL,   0);
    rebx_write_list(rebx, REBX_BINARY_FIELD_TYPE_REGISTERED_PARAM, rebx->registered_params, of);
    REBX_WRITE_FIELD(END,                   NULL,   0);
    
    REBX_WRITE_FIELD(ALLOCATED_FORCES,      NULL,   0);
    rebx_write_list(rebx, REBX_BINARY_FIELD_TYPE_FORCE, rebx->allocated_forces, of);
    REBX_WRITE_FIELD(END,                   NULL,   0);

    
    REBX_WRITE_FIELD(ALLOCATED_OPERATORS,   NULL,   0);
    rebx_write_list(rebx, REBX_BINARY_FIELD_TYPE_OPERATOR, rebx->allocated_operators, of);
    REBX_WRITE_FIELD(END,                   NULL,   0);

    REBX_WRITE_FIELD(ADDITIONAL_FORCES,     NULL,   0);
    rebx_write_list(rebx, REBX_BINARY_FIELD_TYPE_ADDITIONAL_FORCE, rebx->additional_forces, of);
    REBX_WRITE_FIELD(END,                   NULL,   0);

    REBX_WRITE_FIELD(PRE_TIMESTEP_MODIFICATIONS, NULL,   0);
    rebx_write_list(rebx, REBX_BINARY_FIELD_TYPE_STEP, rebx->pre_timestep_modifications, of);
    REBX_WRITE_FIELD(END,                       NULL,   0);
    
    REBX_WRITE_FIELD(POST_TIMESTEP_MODIFICATIONS, NULL,   0);
    rebx_write_list(rebx, REBX_BINARY_FIELD_TYPE_STEP, rebx->post_timestep_modifications, of);
    REBX_WRITE_FIELD(END,                       NULL,   0);

    REBX_WRITE_FIELD(PARTICLES,                 NULL,   0);
    for (int i=0; i<sim->N; i++){
        rebx_write_particle(rebx, &sim->particles[i], i, of);
    }
    REBX_WRITE_FIELD(END,                       NULL,   0);
    
    // Write end marker for binary
    REBX_WRITE_FIELD(END, NULL, 0);
    fclose(of);
}
