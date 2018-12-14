/**
 * @file    core.c
 * @brief   Internal functions for manipulating linked lists
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

#include <string.h>
#include <stdlib.h>
#include "core.h"
#include "rebound.h"
#include "reboundx.h"

struct rebx_node* rebx_attach_node(struct rebx_extras* const rebx, struct rebx_node** head, enum rebx_node_type node_type){
    struct rebx_node* node = rebx_malloc(rebx, sizeof(*node));
    if (node == NULL){
        return NULL;
    }
    node->type = node_type;
    node->next = *head;
    *head = node;
    return node;
}

char* rebx_name_from_node(struct rebx_node* node){
    switch(node->type){
        case REBX_NODE_FORCE:
        {
            struct rebx_force* force = node->object;
            return force->name;
        }
        case REBX_NODE_STEP:
        {
            struct rebx_step* step = node->object;
            return step->operator->name;
        }
        case REBX_NODE_PARAM:
        {
            struct rebx_param* param = node->object;
            return param->name;
        }
    }
    return NULL;
}

struct rebx_node* rebx_get_node(struct rebx_node* head, const char* name){
    struct rebx_node* current = head;
    char* nodename;
    while(current != NULL){
        nodename = rebx_name_from_node(current);
        if(strcmp(nodename, name) == 0){
            return current;
        }
        current = current->next;
    }
    
    return NULL;   // name not found.
}

int rebx_remove_node(struct rebx_extras* const rebx, struct rebx_node** head, const char* name){
    // TODO free memory for deleted node (will need sim)
    struct rebx_node* current = *head;
    char* nodename = rebx_name_from_node(current);
    if(strcmp(nodename, name) == 0){
        *head = current->next;
        return 1;
    }
    
    while(current->next != NULL){
        nodename = rebx_name_from_node(current->next);
        if(strcmp(nodename, name) == 0){
            current->next = current->next->next;
            return 1;
        }
        current = current->next;
    }
    return 0;
}

int rebx_len(struct rebx_node* head){
    int len = 0;
    struct rebx_node* current = head;
    while(current != NULL){
        len++;
        current = current->next;
    }
    
    return len;
}
