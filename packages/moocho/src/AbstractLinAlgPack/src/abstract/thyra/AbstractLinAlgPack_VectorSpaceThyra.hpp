// ///////////////////////////////////////////////////////////////
// VectorSpaceThyra.hpp
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

#ifndef ALAP_VECTOR_SPACE_Thyra_HPP
#define ALAP_VECTOR_SPACE_Thyra_HPP

#include "AbstractLinAlgPack/src/abstract/interfaces/VectorSpace.hpp"
#include "Thyra_VectorSpaceBase.hpp"

namespace AbstractLinAlgPack {

///
/** <tt>VectorSpace</tt> adapter subclass for <tt>Thyra::VectorSpaceBase<value_type> </tt>.
 *
 * Note that the default copy constructor and assignment operators are
 * allowed which yield in shallow copy, not deep copy.
 */
class VectorSpaceThyra : public VectorSpace {
public:

	/** @name Constructors / initializers */
	//@{

	///
	/** Construct to uninitialized.
	 *
	 * Postconditioins:<ul>
	 * <li><tt>this->thyra_vec().get() == NULL</tt>
	 * </ul>
	 */
	VectorSpaceThyra();
	///
	/** Calls <tt>this->initialize()</tt>.
	 */
	VectorSpaceThyra(
		const Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<value_type> >    &thyra_vec_spc
		,const inner_prod_ptr_t                                                  &inner_prod    = Teuchos::null
		);
	///
	/** Initalize given a smart pointer to a <tt>Thyra::VetorSpace</tt> object.
	 *
	 * @param  thyra_vec_spc  [in] Smart pointer to Thyra vector <tt>this</tt> will adapt.
	 * @param  inner_prod    [in] Smart pointer to an inner product.  If <tt>inner_prod.get()==NULL</tt>
	 *                       then a <tt>InnerProductThyra</tt> object will be used which will
	 *                       point to this.
	 *
	 * Preconditioins:<ul>
	 * <li><tt>thyra_vec_spc.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditioins:<ul>
	 * <li><tt>this->thyra_vec_spc().get() == thyra_vec_spc.get()</tt>
	 * <li>[<tt>inner_prod.get()!=NULL</tt>]
	 *     <tt>this->inner_prod().get()==inner_prod.get()</tt>
	 * <li>[<tt>inner_prod.get()==NULL</tt>]
	 *     <tt>dynamic_cast<const InnerProductThyra*>(this->inner_prod().get()).thyra_vec_spc().get()==thyra_vec_spc.get()</tt>
	 * </ul>
	 */
	void initialize(
		const Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<value_type> >    &thyra_vec_spc
		,const inner_prod_ptr_t                                                  &inner_prod    = Teuchos::null
		);
	///
	/** Set to uninitialized and return smart pointer to the internal <tt>Thyra::VectorSpaceBase<value_type> </tt> object.
	 *
	 * Postconditioins:<ul>
	 * <li><tt>this->thyra_vec_spc().get() == NULL</tt>
	 * </ul>
	 */
	Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<value_type> > set_uninitialized();
	///
	/** Return a (converted) smart pointer to the internal smart pointer to the <tt>Thyra::VectorSpaceBase<value_type> </tt> object.
	 *
	 * If <tt>this->thyra_vec_spc().count() == 1</tt>, then <tt>this</tt>
	 * has sole ownership of the <tt>*this->thyra_vec_spc()</tt> object.
	 */
	const Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<value_type> >& thyra_vec_spc() const;

	//@}

	/** @name Overridden from VectorSpace */
	//@{

	///
	space_ptr_t clone() const;
	///
	bool is_compatible(const VectorSpace& vec_spc ) const;
	///
	bool is_in_core() const;
	///
	index_type dim() const;
	///
	vec_mut_ptr_t create_member() const;
	///
	space_fcty_ptr_t small_vec_spc_fcty() const;
	///
	multi_vec_mut_ptr_t create_members(size_type num_vecs) const;

	//@}

private:

#ifdef DOXYGEN_COMPILE
	const Thyra::VectorSpaceBase<value_type>                              *thyra_vector_space;
#else
	Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<value_type> >  thyra_vec_spc_;
#endif

}; // end class VectorSpaceThyra

// ///////////////////////////////
// Inline functions

inline
const Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<value_type> >&
VectorSpaceThyra::thyra_vec_spc() const
{
	return thyra_vec_spc_;
}

} // end namespace AbstractLinAlgPack

#endif  // ALAP_VECTOR_SPACE_Thyra_HPP
