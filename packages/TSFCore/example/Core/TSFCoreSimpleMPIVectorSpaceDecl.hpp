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

// ////////////////////////////////////////////////////////////////////////////
// TSFCoreSimpleMPIVectorSpaceDecl.hpp

#ifndef TSFCORE_SIMPLE_MPI_VECTOR_SPACE_DECL_HPP
#define TSFCORE_SIMPLE_MPI_VECTOR_SPACE_DECL_HPP

#include "TSFCoreMPIVectorSpaceBase.hpp"

namespace TSFCore {

///
/** Simple MPI-based vector class.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class SimpleMPIVectorSpace : virtual public MPIVectorSpaceBase<Scalar> {
public:

	///
	/** Constructs to a parallel vector with a fixed number of elements per process.
	 *
	 * @param  mpiComm     [in] The MPI communicator.  This object must be maintained
	 *                     by the client the entire time that <tt>this</tt> is in use.
	 * @pram  localSubDim  [in] The number of elements per processor.  This number must
	 *                     be the same on every processor in this simple implementation.
	 */
	SimpleMPIVectorSpace( MPI_Comm mpiComm, const Index localSubDim );

	/** @name Overridden from VectorSpece */
	//@{

	///
	Teuchos::RefCountPtr<Vector<Scalar> > createMember() const;
	///
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > clone() const;

	//@}

	/** @name Overridden from MPIVectorSpaceBase */
	//@{

	///
	MPI_Comm mpiComm() const;
	///
 	Index localSubDim() const;

	//@}

private:

	// //////////////////////////////////////
	// Private data members

	MPI_Comm           mpiComm_;
	Index              localSubDim_;

	// Not defined and not to be called
	SimpleMPIVectorSpace();
	
}; // end class SimpleMPIVectorSpace

} // end namespace TSFCore

#endif // TSFCORE_SIMPLE_MPI_VECTOR_SPACE_DECL_HPP
