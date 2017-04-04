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

// Macro to write a single field to a binary file.
#define REBX_WRITE_FIELD(typename, valueptr, typesize, length) {\
            long pos_field_rewrite = ftell(of);\
            struct rebx_binary_field field = {.type = REBX_BINARY_FIELD_TYPE_##typename, .size = typesize};\
            fwrite(&field,sizeof(field),1,of);\
            long pos_start_field = ftell(of);\
            fwrite(valueptr,field.size,length,of);\
            long pos_end_field = ftell(of);\
            field.size = pos_end_field - pos_start_field;\
            fseek(of, pos_field_rewrite, SEEK_SET);\
            fwrite(&field, sizeof(field), 1, of);\
            fseek(of, 0, SEEK_END);\
        }

static void rebx_write_param(struct rebx_param* param, FILE* of){
    long pos_param_rewrite = ftell(of);
    struct rebx_binary_field param_field = {.type = REBX_BINARY_FIELD_TYPE_PARAM, .size=0};
    fwrite(&param_field, sizeof(param_field), 1, of);
    long pos_start_param = ftell(of);
    
    REBX_WRITE_FIELD(PARAM_TYPE,    &param->param_type,     sizeof(param->param_type),      1);
    REBX_WRITE_FIELD(PYTHON_TYPE,   &param->python_type,    sizeof(param->python_type),     1);
    REBX_WRITE_FIELD(NDIM,          &param->ndim,           sizeof(param->ndim),            1);
    REBX_WRITE_FIELD(NAME,          param->name,            strlen(param->name) + 1,        1); // +1 for \0 at end
    REBX_WRITE_FIELD(SHAPE,         param->shape,           sizeof(*param->shape),          param->ndim);
    REBX_WRITE_FIELD(CONTENTS,      param->contents,        rebx_sizeof(param->param_type), param->size);
    // Write end marker for param
    REBX_WRITE_FIELD(END,           NULL,                   0,                              0);

    // Go back and write size of entire parameter and reset file position to end of file
    long pos_end_param = ftell(of);
    param_field.size = pos_end_param - pos_start_param;
    fseek(of, pos_param_rewrite, SEEK_SET);
    fwrite(&param_field, sizeof(param_field), 1, of);
    fseek(of, 0, SEEK_END);
}
    
static void rebx_write_params(struct rebx_param* ap, FILE* of){
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

static void rebx_write_effect(struct rebx_effect* effect, FILE* of){
    long pos_effect_rewrite = ftell(of);
    struct rebx_binary_field effect_field = {.type = REBX_BINARY_FIELD_TYPE_EFFECT, .size=0};
    fwrite(&effect_field, sizeof(effect_field), 1, of);
    long pos_start_effect = ftell(of);
    
    REBX_WRITE_FIELD(NAME,          effect->name,           strlen(effect->name) + 1,       1); // +1 for \0 at end
    REBX_WRITE_FIELD(FORCE_AS_OPERATOR,          &effect->force_as_operator,           sizeof(effect->force_as_operator),            1);
    REBX_WRITE_FIELD(OPERATOR_ORDER,          &effect->operator_order,           sizeof(effect->operator_order),            1);
    rebx_write_params(effect->ap, of); // Write all parameters
    // Write end marker for effect
    REBX_WRITE_FIELD(END,           NULL,                   0,                              0);

    // Go back and write size of entire parameter and reset file position to end of file
    long pos_end_effect = ftell(of);
    effect_field.size = pos_end_effect - pos_start_effect;
    fseek(of, pos_effect_rewrite, SEEK_SET);
    fwrite(&effect_field, sizeof(effect_field), 1, of);
    fseek(of, 0, SEEK_END);
}

static void rebx_write_particle(struct reb_particle* particle, int index, FILE* of){
    long pos_particle_rewrite = ftell(of);
    struct rebx_binary_field particle_field = {.type = REBX_BINARY_FIELD_TYPE_PARTICLE, .size=0};
    fwrite(&particle_field, sizeof(particle_field), 1, of);
    long pos_start_particle = ftell(of);
    
    REBX_WRITE_FIELD(PARTICLE_INDEX,          &index,           sizeof(index),            1);
    rebx_write_params(particle->ap, of); // Write all parameters
    // Write end marker for particle
    REBX_WRITE_FIELD(END,           NULL,                   0,                              0);
    
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
    
    REBX_WRITE_FIELD(REBX_INTEGRATOR,          &rebx->integrator,           sizeof(rebx->integrator),            1);
    REBX_WRITE_FIELD(END,           NULL,                   0,                              0);
    
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
    
    for (int i=0; i<sim->N; i++){
        rebx_write_particle(&sim->particles[i], i, of);
    }
    
    // Write end marker for binary
    REBX_WRITE_FIELD(END,           NULL,                   0,                              0);
    fclose(of);
}
