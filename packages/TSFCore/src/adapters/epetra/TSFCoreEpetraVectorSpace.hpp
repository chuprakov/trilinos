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
// TSFCoreEpetraVectorSpace.hpp

#ifndef TSFCORE_EPETRA_VECTOR_SPACE_HPP
#define TSFCORE_EPETRA_VECTOR_SPACE_HPP

#include "TSFCoreEpetraTypes.hpp"
#include "TSFCoreMPIVectorSpaceBase.hpp"

namespace TSFCore {

/** \brief Concrete <tt>VectorSpace</tt> adapter subclass for
 * <tt>Epetra_Map</tt> that creates <tt>Epetra</tt> adapted vectors
 * and multi-vectors.
 *
 * This uses an <tt>Epetra_Map</tt> object to implement the
 * <tt>MPIVectorSpaceBase</tt> interface.  By implementing the
 * <tt>MPIVectorSpaceBase</tt> interface, this implementation allows
 * the seemless collaboration of different vectors and multi-vectors
 * through the interface classes <tt>MPIVectorBase</tt> and
 * <tt>MPIMultiVectorBase</tt>.
 *
 * The fact that this class embeds a <tt>Epetra_Map</tt> object means
 * that only maps that have elements of size one can be used to define
 * a vector space.  General <tt>Epetra_BlockMap</tt>s can not be used.
 * This is not a serious limitation since <tt>Epetra_Operator</tt>'s
 * domain and range maps are of type <tt>Epetra_Map</tt>.
 *
 * This class works properly even if Epetra is not compiled with
 * support for MPI (i.e. <tt>HAVE_MPI</tt> is not defined when
 * compiling and linking).  If MPI support is not compiled into
 * Epetra, then the dummy implementation defined in
 * <tt>RTOp_mpi.h</tt> is used instead.
 *
 * The default copy constructor and assignment operators are allowed
 * since they have the correct semantics.
 *
 * \ingroup TSFCore_adapters_Epetra_grp
 */
class EpetraVectorSpace : public MPIVectorSpaceBase<RTOp_value_type> {
public:

	///
	typedef RTOp_value_type Scalar;

	/** @name Constructors / initializers */
	//@{

	/** Constructs to an uninitialized state.
	 *
	 * See the postconditions for <tt>setUninitialized()</tt>.
	 */
	EpetraVectorSpace();

	/** Calls <tt>initialize()</tt>.
	 */
	EpetraVectorSpace(
		const Teuchos::RefCountPtr<const Epetra_Map>  &epetra_map
		);

	/** Initialize given an <tt>Epetra_Map</tt>.
	 *
	 * @param  epetra_map  [in] Smart pointer to <tt>const Epetra_Map</tt> object.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>epetra_map.get()!=NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->epetra_map().get()==epetra_map.get()</tt>
	 * <li> <tt>this->dim()==epetra_map->NumGlobalElements()</tt>
	 * <li> if (<tt>PETRA_COMM_MPI</tt> is defined and <tt>dynamic_cast<const Epetra_MpiComm*>(epetra_map.get()) != NULL</tt>)
	 *      <ul><li><tt>this->mpiComm() == dynamic_cast<const Epetra_MpiComm&>(epetra_map)->Comm()</tt></ul>
	 *      else
	 *      <ul><li><tt>this->mpiComm() == MPI_COMM_NULL</tt></ul>
	 * <li> <tt>this->localOffset() == epetra_map->MinMyGID() - epetra_map->IndexBase()</tt>
	 * <li> <tt>this->localSubdim() == epetra_map->NumMyElements()</tt>
	 * </ul>
	 */
	void initialize(
		const Teuchos::RefCountPtr<const Epetra_Map>  &epetra_map
		);

	/** Set uninitialized and return the underlying <tt>Epetra_Map</tt>
	 *
	 * @param  epetra_map  [in/out]  If <tt>epetra_map!=NULL</tt> on input then
	 *                      <tt>*epetra_map</tt> will be set to <tt>this->epetra_map()</tt>.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->epetra_map().get()==NULL</tt>
	 * <li> <tt>this->dim()==0</tt>
	 * <li> <tt>this->mpiComm() == MPI_COMM_NULL</tt>
	 * <li> <tt>this->localOffset() == -1</tt>
	 * <li> <tt>this->localSubdim() == -1)</tt>
	 * </ul>
	 */
	void setUninitialized(
		Teuchos::RefCountPtr<const Epetra_Map> *epetra_map = NULL
		);

	/** Return a smart pointer to the underlying <tt>Epetra_Map</tt> object.
	 */
	Teuchos::RefCountPtr<const Epetra_Map> epetra_map() const;

	//@}

	/** @name Overridden from VectorSpece */
	//@{

	/// Returns an allocated <tt>EpetraVector</tt> object. 
	Teuchos::RefCountPtr<Vector<Scalar> > createMember() const;
#ifdef TSFCORE_EPETRA_USE_EPETRA_DOMAIN_VECTOR_SPACE
	/// Returns a <tt>EpetraVectorSpaceFactory</tt> object.
	Teuchos::RefCountPtr< const VectorSpaceFactory<Scalar> > smallVecSpcFcty() const;
#endif
	/// Returns an allocated <tt>EpetraMultiVector</tt> object. 
	Teuchos::RefCountPtr< MultiVector<Scalar> > createMembers(int numMembers) const;
	/// Creates a deep copy
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > clone() const;

	//@}

	/** @name Overriddend from MPIVectorSpaceBase */
	//@{
	
	/// Overridden
	MPI_Comm mpiComm() const;
	/// Overridden
	Index localSubDim() const;

	//@}

private:

	Teuchos::RefCountPtr<const Epetra_Map>                 epetra_map_;
	MPI_Comm                                               mpiComm_;
	Index                                                  localSubDim_;
#ifdef TSFCORE_EPETRA_USE_EPETRA_DOMAIN_VECTOR_SPACE
	Teuchos::RefCountPtr<const EpetraVectorSpaceFactory>   smallVecSpcFcty_;
#endif

}; // end class EpetraVectorSpace

// //////////////////////////////////////////////
// Inline members

inline
Teuchos::RefCountPtr<const Epetra_Map>
EpetraVectorSpace::epetra_map() const
{
	return epetra_map_;
}

} // end namespace TSFCore

#endif // TSFCORE_EPETRA_VECTOR_SPACE_HPP
