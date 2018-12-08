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

static void rebx_push_node(struct rebx_node** head, struct rebx_node* node){
    node->next = *head;
    *head = node;
}

struct rebx_node* rebx_add_node(struct reb_simulation* const sim, struct rebx_node** head, void* object, const char* name){
    struct rebx_node* node = rebx_malloc(sim, sizeof(*node));
    if (node == NULL){
        return NULL;
    }
    node->name = rebx_malloc(sim, strlen(name) + 1); // +1 for \0 at end
    if (node->name == NULL){
        free(node);
        return NULL;
    }
    else{
        strcpy(node->name, name);
    }
    
    node->hash = reb_hash(name);
    node->object = object;
    rebx_push_node(head, node);
    return node;
}

struct rebx_node* rebx_get_node(struct rebx_node* head, uint32_t hash){
    struct rebx_node* current = head;
    while(current != NULL){
        if(current->hash == hash){
            return current;
        }
        current = current->next;
    }
    
    if (current == NULL){   // param_name not found.
        return NULL;
    }
    
    return current;
}

int rebx_remove_node(struct reb_simulation* const sim, struct rebx_node** head, uint32_t hash){
    // TODO free memory for deleted node (will need sim)
    struct rebx_node* current = *head;
    
    if(current->hash == hash){
        *head = current->next;
        return 1;
    }
    
    while(current->next != NULL){
        if(current->next->hash == hash){
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
