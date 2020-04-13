/**
 * @file    linkedlist.h
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

#ifndef _REBX_LINKEDLIST_H
#define _REBX_LINKEDLIST_H

/**
 * @brief Adds node to linked list
 * @details Caller is responsible for ensuring that the node passed is correctly initialized with name
 * @param head Pointer to the head of the linked list to be modified, e.g. &ps[1]->ap
 * @param node Pointer to an initialized rebx_node
 */
void rebx_add_node(struct rebx_node** head, struct rebx_node* node);
/**
 * @brief Gets node from linked list by name
 * @param head Head of linked list, e.g. ps[1]->ap
 * @param name Name of the node to retrieve
 * @return Pointer to the first node found with passed name, or NULL if none is found
 */
struct rebx_node* rebx_get_node(struct rebx_node* head, const char* name);

/**
 * @brief Removes node from linked list by name
 * @param head Pointer to the head of the linked list to be modified, e.g. &ps[1]->ap
 * @param name Name of the node to remove
 * @return 1 if node found and removed, 0 if not found
 */
int rebx_remove_node(struct rebx_node** head, void* object);

/**
 * @brief Get length of linked list
 * @param head Pointer to the head of the linked list, e.g. &ps[1]->ap
 * @return Length of linked list.
 */
int rebx_len(struct rebx_node* head);

#endif
