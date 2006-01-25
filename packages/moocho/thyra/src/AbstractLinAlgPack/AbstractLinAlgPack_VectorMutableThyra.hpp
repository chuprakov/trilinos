// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef ALAP_VECTOR_MUTABLE_Thyra_HPP
#define ALAP_VECTOR_MUTABLE_Thyra_HPP

#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorApplyOpSerialBase.hpp"
#include "AbstractLinAlgPack_VectorSpaceThyra.hpp"
#include "Thyra_VectorBase.hpp"

namespace AbstractLinAlgPack {

///
/** <tt>VectorMutable</tt> adapter subclass for <tt>Thyra::VectorBase</tt>.
 */
class VectorMutableThyra : public VectorMutable, private VectorApplyOpSerialBase {
public:

	/** @name Constructors / Initializers */
	//@{

	///
	/** Construct to uninitialized.
	 *
	 * Postconditioins:<ul>
	 * <li><tt>this->thyra_vec().get() == NULL</tt>
	 * </ul>
	 */
	VectorMutableThyra();
	///
	/** Calls <tt>this->initialize()</tt>.
	 */
	VectorMutableThyra( const Teuchos::RefCountPtr<Thyra::VectorBase<value_type> >& thyra_vec );
	///
	/** Initalize given a smart pointer to a <tt>Thyra::Vetor</tt> object.
	 *
	 * @param  thyra_vec  [in] Smart pointer to Thyra vector <tt>this</tt> will adapt.
	 *
	 * Preconditioins:<ul>
	 * <li><tt>thyra_vec.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditioins:<ul>
	 * <li><tt>this->thyra_vec().get() == thyra_vec.get()</tt>
	 * </ul>
	 */
	void initialize( const Teuchos::RefCountPtr<Thyra::VectorBase<value_type> >& thyra_vec );
	///
	/** Set to uninitialized and return smart pointer to the internal <tt>Thyra::VectorBase</tt> object.
	 *
	 * Postconditioins:<ul>
	 * <li><tt>this->thyra_vec().get() == NULL</tt>
	 * </ul>
	 */
	Teuchos::RefCountPtr<Thyra::VectorBase<value_type> > set_uninitialized();
	///
	/** Return a (converted) smart pointer to the internal smart pointer to the <tt>Thyra::VectorBase</tt> object.
	 *
	 * If <tt>this->thyra_vec().count() == 2</tt>, then <tt>this</tt>
	 * has so ownership of the <tt>*this->thyra_vec()</tt> object.
	 */
	Teuchos::RefCountPtr<const Thyra::VectorBase<value_type> > thyra_vec() const;

	//@}

	/** @name Methods overridden from Vector */
	//@{

	///
	const VectorSpace& space() const;
	///
	void apply_op(
		const RTOpPack::RTOp       &op
		,const size_t              num_vecs
		,const Vector*             vecs[]
		,const size_t              num_targ_vecs
		,VectorMutable*            targ_vecs[]
		,RTOpPack::ReductTarget    *reduct_obj
		,const index_type          first_ele
		,const index_type          sub_dim
		,const index_type          global_offset
		) const;
	///
	index_type dim() const;
	///
	void get_sub_vector( const Range1D& rng, RTOpPack::SubVector* sub_vec ) const;
	///
	void free_sub_vector( RTOpPack::SubVector* sub_vec ) const;

	//@}

	/** @name Methods overridden from VectorMutable */
	//@{

	///
	void get_sub_vector( const Range1D& rng, RTOpPack::MutableSubVector* sub_vec );
	///
	void commit_sub_vector( RTOpPack::MutableSubVector* sub_vec );
	///
	void set_sub_vector( const RTOpPack::SparseSubVector& sub_vec );

	//@}

private:

#ifdef DOXYGEN_COMPILE
	Thyra::VectorBase<value_type>                          *thyra_vector;
#else
	Teuchos::RefCountPtr<Thyra::VectorBase<value_type> >   thyra_vec_;
#endif
	VectorSpaceThyra                                       space_;

	// Not defined and not to be called!
	//VectorMutableThyra(const VectorMutableThyra&);
	VectorMutableThyra& operator=(const VectorMutableThyra&);

}; // class VectorMutableThyra 

// ////////////////////////////////
// Inline functions

inline
Teuchos::RefCountPtr<const Thyra::VectorBase<value_type> >
VectorMutableThyra::thyra_vec() const
{
	return thyra_vec_;
}

} // end namespace AbstractLinAlgPack

#endif  // ALAP_VECTOR_MUTABLE_Thyra_HPP
