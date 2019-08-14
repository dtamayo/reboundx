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

// Assumes caller has checked head and node != NULL, since theye're the ones that can give meaningful errors
void rebx_add_node(struct rebx_node** head, struct rebx_node* node){
    node->next = *head;
    *head = node;
    return;
}
/*
struct rebx_node* rebx_get_node(struct rebx_node* head, const char* name){
    struct rebx_node* current = head;
    while(current != NULL){
        if(strcmp(current->name, name) == 0){
            return current;
        }
        current = current->next;
    }
    
    return NULL;   // name not found.
}

int rebx_remove_node(struct rebx_node** head, const char* name){
    struct rebx_node* current = *head;
    // Treat edge case where head node is the one to be removed
    if(strcmp(current->name, name) == 0){
        *head = current->next;
        return 1;
    }
    
    while(current->next != NULL){
        if(strcmp(current->next->name, name) == 0){
            current->next = current->next->next;
            return 1;
        }
        current = current->next;
    }
    return 0;
}
*/

// Pass head of linked list and a pointer to the node->object to remove (e.g. force)
int rebx_remove_node(struct rebx_node** head, void* object){
    if (*head == NULL){
        return 0;
    }
    
    struct rebx_node* current = *head;
    if(current->object == object){ // edge case where force is first in list
        *head = current->next;
        free(current);
        return 1;
    }
    
    struct rebx_node* prev = current;
    current = current->next;
    while (current != NULL){
        if(current->object == object){
            prev->next = current->next;
            free(current);
            return 1;
        }
        prev = current;
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
