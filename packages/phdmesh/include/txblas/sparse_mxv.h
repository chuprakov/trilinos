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
 * @author  H. Carter Edwards
 */

#ifndef txblas_sparse_mxv_h
#define txblas_sparse_mxv_h

#if defined( __cplusplus )
extern "C" {
#endif

/*--------------------------------------------------------------------*/

typedef struct txblas_SparseMatrixEntry_struct {
  double   val ;
  unsigned row ;
  unsigned col ;
} txblas_SparseMatrixEntry ;

/*--------------------------------------------------------------------*/
/* Preprocess the matrix entries and generate blocking data
 * for subsequent matrix-vector multiplication.
 * The 'blocking' array must be dimensioned to at least:
 *   blocking[ 1 + ceil( number_row / block_size ) ]
 */
int txblas_rbcr_prep(
  const unsigned number_row ,
  const unsigned block_size ,
  const unsigned number_coef ,
        txblas_SparseMatrixEntry a[] ,
        unsigned blocking[] );

void txblas_rbcr_mxv(
  const unsigned number_row ,
  const unsigned block_size ,
  const txblas_SparseMatrixEntry a[] ,
  const unsigned blocking[] ,
  const double   x[] ,
        double   y[] );

/*--------------------------------------------------------------------*/
/* Preprocess the matrix entries and generate blocking data
 * for subsequent matrix-vector multiplication.
 * The 'blocking' array must be dimensioned to at least:
 *   blocking[ number_row ]
 */
int txblas_cr_prep(
  const unsigned number_row ,
  const unsigned number_coef ,
        txblas_SparseMatrixEntry a[] ,
        unsigned blocking[] );

void txblas_cr_mxv(
  const unsigned number_row ,
  const txblas_SparseMatrixEntry a[] ,
  const unsigned blocking[] ,
  const double   x[] ,
        double   y[] );

/*--------------------------------------------------------------------*/

#if defined( __cplusplus )
} /* extern "C" */
#endif

#endif

