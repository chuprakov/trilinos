// ///////////////////////////////////////////////////////////////
// InnerProductThyra.hpp
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

#ifndef ALAP_INNER_PRODUCT_Thyra_H
#define ALAP_INNER_PRODUCT_Thyra_H

#include "AbstractLinAlgPack/src/abstract/interfaces/InnerProduct.hpp"
#include "Thyra_VectorSpaceBase.hpp"

namespace AbstractLinAlgPack {

///
/** Implements the inner product using <tt>Thyra::VectorSpaceBase::scalarProd()</tt>.
 */
class InnerProductThyra : public InnerProduct {
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
	InnerProductThyra();
	///
	/** Calls <tt>this->initialize()</tt>.
	 */
	InnerProductThyra( const Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<value_type> > &thyra_vec_spc );
	///
	/** Initalize given a smart pointer to a <tt>Thyra::VetorSpace</tt> object.
	 *
	 * @param  thyra_vec_spc  [in] Smart pointer to Thyra vector
	 *
	 * Preconditioins:<ul>
	 * <li><tt>thyra_vec_spc.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditioins:<ul>
	 * <li><tt>this->thyra_vec_spc().get() == thyra_vec_spc.get()</tt>
	 * </ul>
	 */
	void initialize( const Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<value_type> > &thyra_vec_spc );
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

	/** @name Overridden from InnerProduct */
	//@{
	///
	value_type inner_prod(const Vector& v1, const Vector& v2) const;
	//@}

private:

#ifdef DOXYGEN_COMPILE
	const Thyra::VectorSpaceBase<value_type>                              *thyra_vector_space;
#else
	Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<value_type> >  thyra_vec_spc_;
#endif

}; // end class InnerProductThyra

// ///////////////////////////////
// Inline functions

inline
const Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<value_type> >&
InnerProductThyra::thyra_vec_spc() const
{
	return thyra_vec_spc_;
}

} // end namespace AbstractLinAlgPack

#endif  // ALAP_INNER_PRODUCT_Thyra_H
