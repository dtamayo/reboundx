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

static void rebx_write_param(struct rebx_param* param, FILE* of){
    long pos_param_rewrite = ftell(of);
    struct rebx_binary_field param_field = {.type = REBX_BINARY_FIELD_TYPE_PARAM, .size=0};
    printf("before param_field %lu\n", ftell(of));
    fwrite(&param_field, sizeof(param_field), 1, of);
    printf("after param_field %lu\n", ftell(of));
    long pos_start_param = ftell(of);
    
    struct rebx_binary_field field;
    long pos_field_rewrite, pos_start_field, pos_end_field;
    
    // Write simple fields we know the length of ahead of time
    
    field.type = REBX_BINARY_FIELD_TYPE_PARAM_TYPE;
    field.size = sizeof(param->param_type);
    printf("before param_type_field %lu\n", ftell(of));
    fwrite(&field, sizeof(field), 1, of);
    printf("after param_type_field %lu\n", ftell(of));
    fwrite(&param->param_type, sizeof(param->param_type), 1, of);
    printf("after param_type value %lu\n", ftell(of));
    
    field.type = REBX_BINARY_FIELD_TYPE_NDIM;
    field.size = sizeof(param->ndim);
    printf("before ndim field %lu\n", ftell(of));
    fwrite(&field, sizeof(field), 1, of);
    printf("before ndim value %lu\n", ftell(of));
    fwrite(&param->ndim, sizeof(param->ndim), 1, of);
    printf("after ndim value %lu\n", ftell(of));
    
    field.type = REBX_BINARY_FIELD_TYPE_NAMELENGTH;
    size_t namelength = strlen(param->name) + 1; // +1 for \0 at end
    field.size = sizeof(namelength);
    printf("before namelength field %lu\n", ftell(of));
    fwrite(&field, sizeof(field), 1, of);
    printf("before namelength value %lu\n", ftell(of));
    fwrite(&namelength, sizeof(namelength), 1, of);
    printf("after namelength value %lu\n", ftell(of));
    
    // Write variable-sized fields by caching file position and then overwriting at end
    
    field.type = REBX_BINARY_FIELD_TYPE_NAME; // update size at the end
    pos_field_rewrite = ftell(of);
    fwrite(&field, sizeof(field), 1, of);
    pos_start_field = ftell(of);
    printf("before name %lu\n", ftell(of));
    fwrite(param->name, sizeof(*param->name), namelength, of);
    pos_end_field = ftell(of);
    field.size = pos_end_field - pos_start_field;
    //printf("%lu\t%lu\t%lu\n", pos_start_effect, pos_end_effect, field.size);
    fseek(of, pos_field_rewrite, SEEK_SET);
    fwrite(&field, sizeof(field), 1, of);
    fseek(of, 0, SEEK_END);
    
    if(param->ndim > 0){
        field.type = REBX_BINARY_FIELD_TYPE_SHAPE;
        pos_field_rewrite = ftell(of);
        fwrite(&field, sizeof(field), 1, of);
        pos_start_field = ftell(of);
        printf("before shape %lu\n", ftell(of));
        fwrite(param->shape, sizeof(*param->shape), param->ndim, of);
        printf("after shape %lu\n", ftell(of));
        pos_end_field = ftell(of);
        field.size = pos_end_field - pos_start_field;
        //printf("%lu\t%lu\t%lu\n", pos_start_effect, pos_end_effect, field.size);
        fseek(of, pos_field_rewrite, SEEK_SET);
        fwrite(&field, sizeof(field), 1, of);
        fseek(of, 0, SEEK_END);
    }
    
    // CONTENTS must come after PARAM_TYPE, NAME, NDIM AND SHAPE
    field.type = REBX_BINARY_FIELD_TYPE_CONTENTS;
    pos_field_rewrite = ftell(of);
    fwrite(&field, sizeof(field), 1, of);
    pos_start_field = ftell(of);
    fwrite(param->contents, rebx_sizeof(param->param_type), param->size, of);
    printf("after contents %lu\n", ftell(of));
    pos_end_field = ftell(of);
    field.size = pos_end_field - pos_start_field;
    //printf("%lu\t%lu\t%lu\n", pos_start_effect, pos_end_effect, field.size);
    fseek(of, pos_field_rewrite, SEEK_SET);
    fwrite(&field, sizeof(field), 1, of);
    fseek(of, 0, SEEK_END);

    // Write end marker
    field.type = REBX_BINARY_FIELD_TYPE_END;
    field.size = 0;
    fwrite(&field, sizeof(field), 1, of);
    
    // Go back and write size of entire parameter and reset file position to end of file
    long pos_end_param = ftell(of);
    param_field.size = pos_end_param - pos_start_param;
    //printf("%lu\t%lu\t%lu\n", pos_start_effect, pos_end_effect, field.size);
    fseek(of, pos_param_rewrite, SEEK_SET);
    fwrite(&param_field, sizeof(param_field), 1, of);
    fseek(of, 0, SEEK_END);
    //printf("%lu\n", ftell(of));
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
    printf("before effect_field %lu\n", ftell(of));
    fwrite(&effect_field, sizeof(effect_field), 1, of);
    printf("after effect_field %lu\n", ftell(of));
    long pos_start_effect = ftell(of);
    
    struct rebx_binary_field field;
    long pos_field_rewrite, pos_start_field, pos_end_field;

    // Write simple fields we know the length of ahead of time
    
    field.type = REBX_BINARY_FIELD_TYPE_NAMELENGTH;
    size_t namelength = strlen(effect->name) + 1; // +1 for \0 at end
    field.size = sizeof(namelength);
    printf("before namelength field %lu\n", ftell(of));
    fwrite(&field, sizeof(field), 1, of);
    printf("before namelength value %lu\n", ftell(of));
    fwrite(&namelength, sizeof(namelength), 1, of);
    printf("after namelength value %lu\n", ftell(of));
    
    // Write variable-sized fields by caching file position and then overwriting at end
    
    field.type = REBX_BINARY_FIELD_TYPE_NAME; // update size at the end
    pos_field_rewrite = ftell(of);
    fwrite(&field, sizeof(field), 1, of);
    pos_start_field = ftell(of);
    printf("before name %lu\n", ftell(of));
    fwrite(effect->name, sizeof(*effect->name), namelength, of);
    pos_end_field = ftell(of);
    field.size = pos_end_field - pos_start_field;
    //printf("%lu\t%lu\t%lu\n", pos_start_effect, pos_end_effect, field.size);
    fseek(of, pos_field_rewrite, SEEK_SET);
    fwrite(&field, sizeof(field), 1, of);
    fseek(of, 0, SEEK_END);

    // Write all parameters
    rebx_write_params(effect->ap, of);
    
    // Write end marker
    field.type = REBX_BINARY_FIELD_TYPE_END;
    field.size = 0;
    fwrite(&field, sizeof(field), 1, of);
    
    // Go back and write size of entire parameter and reset file position to end of file
    long pos_end_effect = ftell(of);
    effect_field.size = pos_end_effect - pos_start_effect;
    //printf("%lu\t%lu\t%lu\n", pos_start_effect, pos_end_effect, field.size);
    fseek(of, pos_effect_rewrite, SEEK_SET);
    fwrite(&effect_field, sizeof(effect_field), 1, of);
    fseek(of, 0, SEEK_END);
    //printf("%lu\n", ftell(of));
    
    /*//printf("bef effect %lu\n", ftell(of));
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
    //printf("%lu\n", ftell(of));*/
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
    
    // Write end marker for binary
    struct rebx_binary_field field = {.type = REBX_BINARY_FIELD_TYPE_END, .size=0};
    fwrite(&field, sizeof(field), 1, of);
    fclose(of);
}