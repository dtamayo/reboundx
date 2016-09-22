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
#include "rebound.h"
#include "reboundx.h"
#include "core.h"

static struct rebx_extras* rebx_create_extras_from_binary_with_messages(struct rebx_extras* rebx, const char* const filename, enum reb_input_binary_messages* warnings){
    FILE* inf = fopen(filename,"rb"); 
    if (inf){
        long objects = 0;
        // Input header.
        const char str[] = "REBOUNDx Binary File. Version: ";
        char readbuf[65], curvbuf[65];
        sprintf(curvbuf,"%s%s",str,rebx_version_str);
        for(size_t j=strlen(curvbuf);j<63;j++){
            curvbuf[j] = ' ';
        }
        curvbuf[63] = '\0';
        objects += fread(readbuf,sizeof(char),64,inf);
        if(strcmp(readbuf,curvbuf)!=0){
            *warnings |= REB_INPUT_BINARY_WARNING_VERSION;
        }
        printf("%s\n", readbuf);
        struct rebx_binary_field field;
        size_t elements_read = fread(&field, sizeof(field), 1, inf);
        switch (field.type){
            case REBX_BINARY_FIELD_TYPE_EFFECT:
            {
                break;
            }
            case REBX_BINARY_FIELD_TYPE_PARTICLE:
            {
                break;
            }
            case REBX_BINARY_FIELD_TYPE_PARAM:
            {
                break;
            }
            default:
                *warnings |= REB_INPUT_BINARY_WARNING_FIELD_UNKOWN;
                fseek(inf,field.size,SEEK_CUR); // type unrecognized (diff version?) try skipping
                break;
        }
        
        fclose(inf);
        
        return rebx;
    }
    *warnings |= REB_INPUT_BINARY_ERROR_NOFILE;
    return NULL;
}

struct rebx_extras* rebx_create_extras_from_binary(struct reb_simulation* sim, const char* const filename){
    enum reb_input_binary_messages warnings = REB_INPUT_BINARY_WARNING_NONE;
    struct rebx_extras* rebx = rebx_init(sim);
    rebx_create_extras_from_binary_with_messages(rebx, filename, &warnings);
    
    if (warnings & REB_INPUT_BINARY_WARNING_VERSION){
        reb_warning(sim,"REBOUNDx: Binary file was saved with a different version of REBOUNDx. Binary format might have changed and corrupted the loading. Check that effects and parameters are loaded as expected.");
    }
    if (warnings & REB_INPUT_BINARY_WARNING_FIELD_UNKOWN){
        reb_warning(sim,"REBOUNDx: Unknown field found in binary file. Any unknown fields not loaded.  This can happen if the binary was created with a later version of REBOUNDx than the one used to read it.");
    }
    if (warnings & REB_INPUT_BINARY_ERROR_NOFILE){
        reb_error(sim,"REBOUNDx: Cannot read binary file. Check filename and file contents.");
        rebx_free(rebx);
        rebx = NULL;
    }
    return rebx;
}
