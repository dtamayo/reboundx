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

static void rebx_write_param(struct rebx_extras* rebx, struct rebx_param* param, FILE* of){
    long pos_param_rewrite = ftell(of);
    struct rebx_binary_field param_field = {.type = REBX_BINARY_FIELD_TYPE_PARAM, .size=0};
    fwrite(&param_field, sizeof(param_field), 1, of);
    long pos_start_param = ftell(of);
    
    REBX_WRITE_FIELD(PARAM_TYPE, &param->type,     sizeof(param->type));
    REBX_WRITE_FIELD(NAME,       param->name,      strlen(param->name) + 1); // +1 for \0 at end
    if(param->value){ // will be NULL for registered_params so don't write
        REBX_WRITE_FIELD(VALUE,      param->value,     rebx_sizeof(rebx, param->type));
    }
    REBX_WRITE_FIELD(END,        NULL,             0);

    // Go back and write size of entire parameter and reset file position to end of file
    long pos_end_param = ftell(of);
    param_field.size = pos_end_param - pos_start_param;
    fseek(of, pos_param_rewrite, SEEK_SET);
    fwrite(&param_field, sizeof(param_field), 1, of);
    fseek(of, 0, SEEK_END);
}
    
static void rebx_write_params(struct rebx_extras* rebx, struct rebx_node* ap, FILE* of){
    int nparams = rebx_len(ap);
    
    // Write params to binary in reverse order so we get same order back when read back in
    while (nparams > 0){
        struct rebx_node* current = ap;
        for(int i=0;i<nparams-1;i++){
            current=current->next;
        }
        rebx_write_param(rebx, current->object, of);
        nparams--;
    }
}

static void rebx_write_force(struct rebx_extras* rebx, struct rebx_force* force, FILE* of){
    long pos_effect_rewrite = ftell(of);
    struct rebx_binary_field force_field = {.type = REBX_BINARY_FIELD_TYPE_FORCE, .size=0};
    fwrite(&force_field, sizeof(force_field), 1, of);
    long pos_start_effect = ftell(of);
    
    REBX_WRITE_FIELD(NAME,          force->name,        strlen(force->name) + 1); // +1 for \0 at end
    REBX_WRITE_FIELD(FORCE_TYPE,    &force->force_type, sizeof(force->force_type));
    rebx_write_params(rebx, force->ap, of); // Write all parameters
    REBX_WRITE_FIELD(END,           NULL,               0);

    // Go back and write size of entire parameter and reset file position to end of file
    long pos_end_effect = ftell(of);
    force_field.size = pos_end_effect - pos_start_effect;
    fseek(of, pos_effect_rewrite, SEEK_SET);
    fwrite(&force_field, sizeof(force_field), 1, of);
    fseek(of, 0, SEEK_END);
}

// Same as force, but only holds the name for later loading, rather than the whole parameter list
static void rebx_write_additional_force(struct rebx_extras* rebx, struct rebx_force* force, FILE* of){
    long pos_effect_rewrite = ftell(of);
    struct rebx_binary_field force_field = {.type = REBX_BINARY_FIELD_TYPE_FORCE, .size=0};
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
    REBX_WRITE_FIELD(OPERATOR_TYPE, &operator->operator_type,   sizeof(operator->operator_type));
    rebx_write_params(rebx, operator->ap, of); // Write all parameters
    // Write end marker for effect
    REBX_WRITE_FIELD(END,           NULL,                       0);
    
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
    
    REBX_WRITE_FIELD(NAME,          step->name,             strlen(step->name) + 1);
    REBX_WRITE_FIELD(DT_FRACTION,   &step->dt_fraction,     sizeof(step->dt_fraction));
    // Need operator name to link it back when reading it back in
    REBX_WRITE_FIELD(OPERATOR_NAME, step->operator->name,   strlen(step->operator->name) + 1);
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
    rebx_write_params(rebx, particle->ap, of); // Write all parameters
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
    
    REBX_WRITE_FIELD(REBX_INTEGRATOR,   &rebx->integrator,  sizeof(rebx->integrator));
    REBX_WRITE_FIELD(END,               NULL,               0);
    
    // Go back and write size of entire parameter and reset file position to end of file
    long pos_end_rebx = ftell(of);
    rebx_field.size = pos_end_rebx - pos_start_rebx;
    fseek(of, pos_rebx_rewrite, SEEK_SET);
    fwrite(&rebx_field, sizeof(rebx_field), 1, of);
    fseek(of, 0, SEEK_END);
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
    rebx_write_params(rebx, rebx->registered_params, of);
    
    // Write all forces allocated
    int Nallocatedforces = rebx_len(rebx->allocated_forces);
    while (Nallocatedforces > 0){
        struct rebx_node* current = rebx->allocated_forces;
        for(int i=0;i<Nallocatedforces-1;i++){
            current=current->next;
        }
        rebx_write_force(rebx, current->object, of);
        Nallocatedforces--;
    }
    
    // Write all operators allocated
    int Nallocatedoperators = rebx_len(rebx->allocated_operators);
    while (Nallocatedoperators > 0){
        struct rebx_node* current = rebx->allocated_operators;
        for(int i=0;i<Nallocatedoperators-1;i++){
            current=current->next;
        }
        rebx_write_operator(rebx, current->object, of);
        Nallocatedoperators--;
    }
    
    // Write additional_forces
    int Nadditionalforces = rebx_len(rebx->additional_forces);
    while (Nadditionalforces > 0){
        struct rebx_node* current = rebx->additional_forces;
        for(int i=0;i<Nadditionalforces-1;i++){
            current=current->next;
        }
        rebx_write_additional_force(rebx, current->object, of);
        Nadditionalforces--;
    }
    
    // Write pre_timestep_modifications
    int Npretm = rebx_len(rebx->pre_timestep_modifications);
    while (Npretm > 0){
        struct rebx_node* current = rebx->pre_timestep_modifications;
        for(int i=0;i<Npretm-1;i++){
            current=current->next;
        }
        rebx_write_step(rebx, current->object, of);
        Npretm--;
    }
    
    // Write post_timestep_modifications
    int Nposttm = rebx_len(rebx->post_timestep_modifications);
    while (Nposttm > 0){
        struct rebx_node* current = rebx->post_timestep_modifications;
        for(int i=0;i<Nposttm-1;i++){
            current=current->next;
        }
        rebx_write_step(rebx, current->object, of);
        Nposttm--;
    }
    
    // Write particle parameters
    for (int i=0; i<sim->N; i++){
        rebx_write_particle(rebx, &sim->particles[i], i, of);
    }
    
    // Write end marker for binary
    REBX_WRITE_FIELD(END, NULL, 0);
    fclose(of);
}
