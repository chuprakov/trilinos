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

// ///////////////////////////////////////////////////////////////
// TSFCoreMPIVectorSpaceFactoryStdDecl.hpp

#ifndef TSFCORE_MPI_VECTOR_SPACE_FACTORY_STD_DECL_HPP
#define TSFCORE_MPI_VECTOR_SPACE_FACTORY_STD_DECL_HPP

#include "TSFCoreVectorSpaceFactory.hpp"

namespace TSFCore {

///
/** Implementation of a vector-space factory for a <tt>MPIVectorSpaceStd</tt> objects.
 *
 * This will create either serial (<tt>mpiComm==MPI_COMM_NULL</tt>) or locally
 * replicated (<tt>mpiComm!=MPI_COMM_NULL</tt>) vector space objects.
 */
template<class Scalar>
class MPIVectorSpaceFactoryStd : public VectorSpaceFactory<Scalar> {
public:

	///
	MPIVectorSpaceFactoryStd( MPI_Comm  mpiComm );

	/** @name Overridden from VectorSpaceFactory */
	//@{
	///
	Teuchos::RefCountPtr<const VectorSpace<Scalar> > createVecSpc(int dim) const;
	//@}

private:

	MPI_Comm  mpiComm_;

  MPIVectorSpaceFactoryStd(); // Not defined and not to be called!
  	
}; // end class MPIVectorSpaceFactoryStd

} // end namespace TSFCore

#endif  // TSFCORE_MPI_VECTOR_SPACE_FACTORY_STD_DECL_HPP
