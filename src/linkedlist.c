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
#include "rebound.h"
#include "reboundx.h"

void rebx_push_node(struct rebx_node** head, struct rebx_node* node){
    node->next = *head;
    *head = node;
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

struct rebx_node* rebx_add_node(struct reb_simulation* const sim, struct rebx_node** apptr, const char* name){
    
    // Check it doesn't already exist in linked list
    void* ptr = rebx_get_node(*apptr, reb_hash(name));
    if (ptr != NULL){
        char str[300];
        sprintf(str, "REBOUNDx Error: Parameter '%s' passed to rebx_add_param already exists.\n", name);
        reb_error(sim, str);
        return NULL;
    }
    
    // Doesn't exist, allocate and add new one
    struct rebx_node* node = malloc(sizeof(*node));
    node->name = malloc(strlen(name) + 1); // +1 for \0 at end
    if (node->name != NULL){
        strcpy(node->name, name);
    }
    node->hash = reb_hash(name);
    rebx_push_node(apptr, node);
    
    return node;
}

int rebx_remove_node(struct rebx_node** head, uint32_t hash){
    // TODO free memory for deleted node
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
