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
#include "linkedlist.h"

/* Binary format is meant to allow for future structure modifications and backwards compatibility.
 
 It is a nested series of rebx_binary_field structs, with a binary_field_type enum that tells you how to handle what comes next, and a size so that you can skip this object if you don't recognize it (perhaps it's an older version of REBOUNDx reading a newer binary).
 
 At the outermost level you read in REBX_BINARY_FIELD_SNAPSHOT objects. Currently each binary only has one, but this allows us to later make it possible to append more snapshots if needed.
 
 Each snapshot currently holds the rebx structure and a list of particles, for which we store all  the params attached to them. So each of them would have a rebx_binary_field identifying them and telling you how large they are in case you need to skip it.
 
 In principle, the input.c file would have a different function for reading in each of these different types of objects. Each object can have its own set of objects in this nested fashion. Eventually you reach a basic type, whose data we want to read. In that case we use the size_to_skip as the size_to_read for fread, which is the same. These unambiguous blocks don't have an REBX_FIELD_TYPE_END field struct, only the abstract objects whose length is arbitrary (user could add different number of forces, or we could add fields to various structs with code updates).
 
 Clearest in an example,
 
 SNAPSHOT {type=SNAPSHOT, size=skip_to_next_snapshot}
    REBX {type=REBX_STRUCT, size=skip_to_particles}
        ...
    END (REBX)
    PARTICLES {type=PARTICLES, size=skip_to_END(SNAPSHOT)}
        PARTICLE {type=PARTICLE, size=skip_to_next_particle}
            PARTICLE_INDEX {type=PARTICLE_INDEX, size=size_to_read}
            INT
            PARAM_LIST {type=PARAM_LIST, size=skip_to_END(PARTICLE)}
                PARAM {type=PARAM, size=skip_to_next_param}
                    PARAM_TYPE {type=PARAM_TYPE, size=size_to_read}
                    ENUM
                    NAME {type=NAME, size=size_to_read}
                    STRING
                    PARAM_VALUE {type=PARAM_VALUE, size=size_to_read}
                    VALUE
                END (PARAM)
                ...
                PARAM {type=PARAM, size=skip_to_next_param}
                    PARAM_TYPE {type=PARAM_TYPE, size=size_to_read}
                    ENUM
                    NAME {type=NAME, size=size_to_read}
                    STRING
                    PARAM_VALUE {type=PARAM_VALUE, size=size_to_read}
                    VALUE
                END (PARAM)
            END (PARAM_LIST)
        END (PARTICLE)
        ...
        PARTICLE {type=PARTICLE, size=skip_to_next_particle}
        ...
        END (PARTICLE)
    END (PARTICLES)
 END (SNAPSHOT)
 
 // not implemented yet
 SNAPSHOT {type=SNAPSHOT, size=skip_to_next_snapshot}
 ...
 END (SNAPSHOT)
*/

/************************************************************
Macros to remove repetition in writing fields.
*************************************************************/

// Write a data field of binary_field_type typename with size typesize
// valueptr is a pointer to the memory to write
#define REBX_WRITE_DATA_FIELD(typename, valueptr, typesize) {\
struct rebx_binary_field field = {.type = REBX_BINARY_FIELD_TYPE_##typename, .size=typesize};\
fwrite(&field, sizeof(field), 1, of);\
fwrite(valueptr,typesize,1,of);\
}

/*  For the arbitrary objects, we write a preliminary field struct without a size (since we don't know it yet), and cache the file position to measure how large the object is later.*/
#define REBX_START_OBJECT_FIELD(name, typename)\
long pos_start_header_##name = ftell(of);\
struct rebx_binary_field header_##name = {.type = REBX_BINARY_FIELD_TYPE_##typename, .size=0};\
fwrite(&header_##name, sizeof(header_##name), 1, of);\
long pos_start_##name = ftell(of);\

/*  After we write all the data we need for the particular object, we calculate how long this segment is, and update the field struct with this size so we have option of skipping the whole object when reading.*/

#define REBX_END_OBJECT_FIELD(name) {\
REBX_WRITE_DATA_FIELD(END,        NULL,             0);\
long pos_end_##name = ftell(of);\
header_##name.size = pos_end_##name - pos_start_##name;\
fseek(of, pos_start_header_##name, SEEK_SET);\
fwrite(&header_##name, sizeof(header_##name), 1, of);\
fseek(of, 0, SEEK_END);\
}

/*  Write a list of listtype (e.g., ALLOCATED_FORCES) with nodes of type nodetype (e.g. ALLOCATED_FORCE), to the passed linkedlist (e.g. rebx->allocated_forces)*/

#define REBX_WRITE_LIST_FIELD(listtype, nodetype, linkedlist) {\
REBX_START_OBJECT_FIELD(list, listtype);\
rebx_write_list(rebx, REBX_BINARY_FIELD_TYPE_##nodetype, linkedlist, of);\
REBX_END_OBJECT_FIELD(list);\
}

static void rebx_write_list(struct rebx_extras* rebx, enum rebx_binary_field_type list_type, struct rebx_node* list, FILE* of);

static void rebx_write_force_param(struct rebx_extras* rebx, struct rebx_param* param, FILE* of){
    REBX_START_OBJECT_FIELD(force_param, PARAM);
    REBX_WRITE_DATA_FIELD(PARAM_TYPE, &param->type,     sizeof(param->type));
    REBX_WRITE_DATA_FIELD(NAME,       param->name,      strlen(param->name) + 1);
    struct rebx_force* force = param->value;
    REBX_WRITE_DATA_FIELD(PARAM_VALUE,      force->name,      strlen(force->name) + 1);
    REBX_END_OBJECT_FIELD(force_param);
}

static void rebx_write_param(struct rebx_extras* rebx, struct rebx_param* param, FILE* of){
    if (param->type == REBX_TYPE_POINTER){ // Don't write pointers because we won't know how to load them when we read binary. Need to add type to store in binaries.
        return;
    }
    
    if (param->type == REBX_TYPE_FORCE){ // Force already written to allocated_force list. For parce PARAMETERS we agree to store force name in param->value so that the reallocated force can be linked up when we read binary
        rebx_write_force_param(rebx, param, of);
        return;
    }
    REBX_START_OBJECT_FIELD(param, PARAM);
    REBX_WRITE_DATA_FIELD(PARAM_TYPE, &param->type,     sizeof(param->type));
    REBX_WRITE_DATA_FIELD(NAME,       param->name,      strlen(param->name) + 1);
    REBX_WRITE_DATA_FIELD(PARAM_VALUE,      param->value,     rebx_sizeof(rebx, param->type));
    REBX_END_OBJECT_FIELD(param);
}

static void rebx_write_registered_param(struct rebx_extras* rebx, struct rebx_param* param, FILE* of){
    REBX_START_OBJECT_FIELD(registered_param, REGISTERED_PARAM);
    REBX_WRITE_DATA_FIELD(PARAM_TYPE, &param->type,     sizeof(param->type));
    REBX_WRITE_DATA_FIELD(NAME,       param->name,      strlen(param->name) + 1);
    REBX_END_OBJECT_FIELD(registered_param);
}

static void rebx_write_force(struct rebx_extras* rebx, struct rebx_force* force, FILE* of){
    REBX_START_OBJECT_FIELD(force, FORCE);
    // must write name first so that force can be loaded on read
    REBX_WRITE_DATA_FIELD(NAME, force->name, strlen(force->name) + 1);
    REBX_WRITE_LIST_FIELD(PARAM_LIST, PARAM, force->ap);
    REBX_END_OBJECT_FIELD(force);
}

// Same as force, but only holds the name for later loading, rather than the whole parameter list
static void rebx_write_additional_force(struct rebx_extras* rebx, struct rebx_force* force, FILE* of){
    REBX_START_OBJECT_FIELD(additional_force, ADDITIONAL_FORCE);
    REBX_WRITE_DATA_FIELD(NAME, force->name, strlen(force->name) + 1);
    REBX_END_OBJECT_FIELD(additional_force);
}

static void rebx_write_operator(struct rebx_extras* rebx, struct rebx_operator* operator, FILE* of){
    REBX_START_OBJECT_FIELD(operator, OPERATOR);
    REBX_WRITE_DATA_FIELD(NAME, operator->name, strlen(operator->name) + 1);
    REBX_WRITE_LIST_FIELD(PARAM_LIST, PARAM, operator->ap);
    REBX_END_OBJECT_FIELD(operator);
}

static void rebx_write_step(struct rebx_extras* rebx, struct rebx_step* step, FILE* of){
    REBX_START_OBJECT_FIELD(step, STEP);
    // Need operator name to load it from source when reading it back in
    REBX_WRITE_DATA_FIELD(NAME, step->operator->name,   strlen(step->operator->name) + 1);
    REBX_WRITE_DATA_FIELD(STEP_DT_FRACTION,   &step->dt_fraction,     sizeof(step->dt_fraction));
    REBX_END_OBJECT_FIELD(step);
}

static void rebx_write_particle(struct rebx_extras* rebx, struct reb_particle* particle, int index, FILE* of){
    REBX_START_OBJECT_FIELD(particle, PARTICLE);
    REBX_WRITE_DATA_FIELD(PARTICLE_INDEX,    &index, sizeof(index));
    REBX_WRITE_LIST_FIELD(PARAM_LIST, PARAM, particle->ap);
    REBX_END_OBJECT_FIELD(particle);
}

static void rebx_write_rebx(struct rebx_extras* rebx, FILE* of){
    REBX_START_OBJECT_FIELD(rebx_structure, REBX_STRUCTURE);
    REBX_WRITE_LIST_FIELD(REGISTERED_PARAMETERS, REGISTERED_PARAM, rebx->registered_params);
    REBX_WRITE_LIST_FIELD(ALLOCATED_FORCES, FORCE, rebx->allocated_forces);
    REBX_WRITE_LIST_FIELD(ALLOCATED_OPERATORS, OPERATOR, rebx->allocated_operators);
    REBX_WRITE_LIST_FIELD(ADDITIONAL_FORCES, ADDITIONAL_FORCE, rebx->additional_forces);
    REBX_WRITE_LIST_FIELD(PRE_TIMESTEP_MODIFICATIONS, STEP, rebx->pre_timestep_modifications);
    REBX_WRITE_LIST_FIELD(POST_TIMESTEP_MODIFICATIONS, STEP, rebx->post_timestep_modifications);
    REBX_END_OBJECT_FIELD(rebx_structure);
}

// Write a particle field for each particle with a list of its parameters
static void rebx_write_particles(struct rebx_extras* rebx, FILE* of){
    struct reb_simulation* sim = rebx->sim; // checked sim valid in output_binray
    
    REBX_START_OBJECT_FIELD(particle_list, PARTICLES);
    for (int i=0; i<sim->N; i++){
        rebx_write_particle(rebx, &sim->particles[i], i, of);
    }
    REBX_END_OBJECT_FIELD(particle_list);
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

// Could be extended to include time or steps_done to make an archive
static void rebx_write_snapshot(struct rebx_extras* rebx, FILE* of){
    REBX_START_OBJECT_FIELD(snapshot, SNAPSHOT);
    rebx_write_rebx(rebx, of);
    rebx_write_particles(rebx, of);
    REBX_END_OBJECT_FIELD(snapshot);
}

void rebx_output_binary(struct rebx_extras* rebx, char* filename){
    FILE* of = fopen(filename,"wb");
    if (of==NULL){
        rebx_error(rebx, "REBOUNDx error: Can not open file passed to rebx_output_binary.");
    }
    if (rebx->sim == NULL){
        rebx_error(rebx, ""); // rebx_error gives meaningful err
        return;
    }
    // Write header.
    const char str[] = "REBOUNDx Binary File. Version: ";
    char zero = '\0';
    size_t lenheader = strlen(str)+strlen(rebx_version_str);
    fwrite(str,sizeof(char),strlen(str),of);
    fwrite(rebx_version_str,sizeof(char), strlen(rebx_version_str),of);
    fwrite(&zero,sizeof(char),1,of);
    fwrite(rebx_githash_str,sizeof(char),62-lenheader,of);
    fwrite(&zero,sizeof(char),1,of);

    rebx_write_snapshot(rebx, of);
    fclose(of);
}
