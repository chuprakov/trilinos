// ////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_VectorMutableThyra.hpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

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
	VectorMutableThyra(const VectorMutableThyra&);
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
