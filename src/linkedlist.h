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

void rebx_push_node(struct rebx_node** head, struct rebx_node* node);
struct rebx_node* rebx_get_node(struct rebx_node* head, uint32_t hash);
struct rebx_node* rebx_add_node(struct reb_simulation* const sim, struct rebx_node** apptr, const char* name);
int rebx_remove_node(struct reb_simulation* const sim, struct rebx_node** head, uint32_t hash);
#endif
