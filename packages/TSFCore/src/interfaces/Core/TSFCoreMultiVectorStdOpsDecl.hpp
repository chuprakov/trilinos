// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// ////////////////////////////////////////////////////////////////////
// TSFCoreMultiVectorStdOpsDecl.hpp

#ifndef TSFCORE_MULTI_VECTOR_STD_OPS_DECL_HPP
#define TSFCORE_MULTI_VECTOR_STD_OPS_DECL_HPP

#include "TSFCoreTypes.hpp"

namespace TSFCore {

/** \defgroup TSFCORE_MultiVectorStdOps_grp Collection of standard multi-vector operations.
 */
//@{

///
/** Multi-vector dot product.
 *
 * @param  V1   [in]
 * @param  V2   [in]
 * @param  dot  [out] Array (size <tt>V1->domain()->dim()</tt>) of the dot products
 *              <tt>dot[j-1] = dot(*V1.col(j),*V2.col(j))</tt> computed using
 *              a single reduction.
 */
template<class Scalar>
void dot( const MultiVector<Scalar>& V1, const MultiVector<Scalar>& V2, Scalar dot[] );

///
/** Take the one norm of a multi-vector.
 *
 * @param V  [in]
 */
template<class Scalar>
Scalar norm_1( const MultiVector<Scalar>& V );

///
/** V = alpha*V.
 *
 * Note, if alpha==0.0, then V=alpha is performed and if alpha==1.0,
 * then nothing is done.
 */
template<class Scalar>
void scale( Scalar alpha, MultiVector<Scalar>* V );

///
/** V = alpha.
 */
template<class Scalar>
void assign( MultiVector<Scalar>* V, Scalar alpha );

///
/** V = U
 */
template<class Scalar>
void assign( MultiVector<Scalar>* V, const MultiVector<Scalar>& U );

///
/** alpha*U + V -> V
 */
template<class Scalar>
void update( Scalar alpha, const MultiVector<Scalar>& U, MultiVector<Scalar>* V );

///
/** alpha[j-1]*beta*U(j) + V(j) - > V(j), for j = 1 ... U.domain()->dim()
 */
template<class Scalar>
void update( Scalar alpha[], Scalar beta, const MultiVector<Scalar>& U, MultiVector<Scalar>* V );

///
/** U(j) + alpha[j-1]*beta*V(j) - > V(j), for j = 1 ... U.domain()->dim()
 */
template<class Scalar>
void update( const MultiVector<Scalar>& U, Scalar alpha[], Scalar beta, MultiVector<Scalar>* V );

///
/** Generate a random multi-vector with elements uniformly distrubuted
 * elements.
 * 
 * The elements <tt>V->col(j)-getEle(i)</tt> are randomly generated between
 * <tt>[l,u]</tt>.
 *
 * The seed is set using <tt>seed_randomize()</tt>
 */
template<class Scalar>
void randomize( Scalar l, Scalar u, MultiVector<Scalar>* V );

//@}

} // end namespace TSFCore

#endif // TSFCORE_MULTI_VECTOR_STD_OPS_DECL_HPP
