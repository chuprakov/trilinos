/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2007) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 */

#ifndef util_taskpool_h
#define util_taskpool_h

#if defined( __cplusplus )
extern "C" {
#endif

int phdmesh_taskpool_resize( unsigned p_size );

typedef int (*phdmesh_taskpool_routine)( void * routine_data ,
                                         unsigned p_size ,
                                         unsigned p_rank );

/* Execute
 *   (*routine)( routine_data , p_size , p_rank )
 * for each p_rank in [ 0 .. p_size )
 * Enable the requested number of locks.
 * Nested calls to this routine are illegal and immediately returns -1.
 * otherwise the union (bitwise-or) of each routine's return value is returned.
 */
int phdmesh_taskpool_run(
  phdmesh_taskpool_routine routine ,
  void * routine_data ,
  unsigned number_locks );

/** Issue a lock, provide a string to print for an error */
void phdmesh_taskpool_lock( unsigned , const char * const );

/** Unlock, provide a string to print for an error */
void phdmesh_taskpool_unlock( unsigned , const char * const );

#if defined( __cplusplus )
}
#endif

#endif

