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
// TSFCoreEpetraVectorSpaceFactory.hpp

#ifndef TSFCORE_EPETRA_VECTOR_SPACE_FACTORY_HPP
#define TSFCORE_EPETRA_VECTOR_SPACE_FACTORY_HPP

#include "TSFCoreEpetraTypes.hpp"
#include "TSFCoreVectorSpaceFactory.hpp"

namespace TSFCore {

///
/** \brief Concrete <tt>VectorSpaceFactory</tt> adapter subclass for
 * <tt>Epetra_Comm</tt> that creates locally-replicated
 * <tt>EpetraVectorSpace</tt> objects given a dimension.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup TSFCore_adapters_Epetra_grp
 */
class EpetraVectorSpaceFactory : public VectorSpaceFactory<double> {
public:

	///
	typedef double Scalar;

	/** @name Constructors / initializers */
	//@{

	///
	/** Constructs to an uninitialized state.
	 *
	 * See the postconditions for <tt>setUninitialized()</tt>.
	 */
	EpetraVectorSpaceFactory();

	///
	/** Calls <tt>initialize()</tt>.
	 */
	EpetraVectorSpaceFactory(
		const Teuchos::RefCountPtr<const Epetra_Comm>  &epetra_comm
		);

	///
	/** Initialize given a communicator.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>epetra_comm.get()!=NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->epetra_com.get()==epetra_comm.get()</tt>
	 * </ul>
	 */
	void initialize(
		const Teuchos::RefCountPtr<const Epetra_Comm>  &epetra_comm
		);

	///
	/** Set uninitialized and return the underlying <tt>Epetra_BlockMap</tt>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->epetra_map.get()==NULL</tt>
	 * </ul>
	 */
	void setUninitialized(
		Teuchos::RefCountPtr<const Epetra_Comm> *epetra_comm = NULL
		);

	///
	/** 
   * Return the underlying Epetra communicator.
   * (This had been declared but unimplemented. Fixed by KL)
	 */
	Teuchos::RefCountPtr<const Epetra_Comm> epetra_comm() const
  {return epetra_comm_;}

	//@}

	/** @name Overridden from VectorSpaceFactory */
	//@{

	///
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > createVecSpc(int dim) const;

	//@}

private:

#ifdef DOXYGEN_COMPILE
	Epetra_Comm                                     *epetra_comm;
#else	
	Teuchos::RefCountPtr<const Epetra_Comm>         epetra_comm_;
#endif

}; // end class EpetraVectorSpaceFactory

} // end namespace TSFCore

#endif  // TSFCORE_EPETRA_VECTOR_SPACE_FACTORY_HPP
