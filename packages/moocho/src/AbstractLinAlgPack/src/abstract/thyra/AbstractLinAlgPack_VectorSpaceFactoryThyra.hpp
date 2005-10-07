// ///////////////////////////////////////////////////////////////
// VectorSpaceFactoryThyra.hpp
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

#ifndef ALAP_VECTOR_SPACE_FACTORY_Thyra_HPP
#define ALAP_VECTOR_SPACE_FACTORY_Thyra_HPP

#include "AbstractLinAlgPack/src/abstract/interfaces/VectorSpaceFactory.hpp"
#include "Thyra_VectorSpaceFactoryBase.hpp"

namespace AbstractLinAlgPack {

///
/** <tt>VectorSpaceFactory</tt> adapter subclass for <tt>Thyra::VectorSpaceBase</tt>.
 */
class VectorSpaceFactoryThyra : public VectorSpaceFactory {
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
	VectorSpaceFactoryThyra();
	///
	/** Calls <tt>this->initialize()</tt>.
	 */
	VectorSpaceFactoryThyra( const Teuchos::RefCountPtr<const Thyra::VectorSpaceFactoryBase<value_type> > &thyra_vec_spc_fcty );
	///
	/** Initalize given a smart pointer to a <tt>Thyra::VetorSpaceFactory</tt> object.
	 *
	 * @param  thyra_vec_spc_fcty  [in] Smart pointer to Thyra vector
	 *
	 * Preconditioins:<ul>
	 * <li><tt>thyra_vec_spc_fcty.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditioins:<ul>
	 * <li><tt>this->thyra_vec_spc_fcty().get() == thyra_vec_spc_fcty.get()</tt>
	 * </ul>
	 */
	void initialize( const Teuchos::RefCountPtr<const Thyra::VectorSpaceFactoryBase<value_type> > &thyra_vec_spc_fcty );
	///
	/** Set to uninitialized and return smart pointer to the internal <tt>Thyra::VectorSpaceBase</tt> object.
	 *
	 * Postconditioins:<ul>
	 * <li><tt>this->thyra_vec_spc_fcty().get() == NULL</tt>
	 * </ul>
	 */
	Teuchos::RefCountPtr<const Thyra::VectorSpaceFactoryBase<value_type> > set_uninitialized();
	///
	/** Return a (converted) smart pointer to the internal smart pointer to the <tt>Thyra::VectorSpaceBase</tt> object.
	 *
	 * If <tt>this->thyra_vec_spc_fcty().count() == 1</tt>, then <tt>this</tt>
	 * has sole ownership of the <tt>*this->thyra_vec_spc_fcty()</tt> object.
	 */
	const Teuchos::RefCountPtr<const Thyra::VectorSpaceFactoryBase<value_type> >& thyra_vec_spc_fcty() const;

	//@}

	/** @name Overridden from VectorSpaceFactory */
	//@{

	///
	space_ptr_t create_vec_spc(index_type dim) const;

	//@}
	
private:

#ifdef DOXYGEN_COMPILE
	const Thyra::VectorSpaceFactoryBase<value_type>                             *thyra_vector_space_factory;
#else
	Teuchos::RefCountPtr<const Thyra::VectorSpaceFactoryBase<value_type> >  thyra_vec_spc_fcty_;
#endif

}; // end class VectorSpaceFactoryThyra

// ///////////////////////////////
// Inline functions

inline
const Teuchos::RefCountPtr<const Thyra::VectorSpaceFactoryBase<value_type> >&
VectorSpaceFactoryThyra::thyra_vec_spc_fcty() const
{
	return thyra_vec_spc_fcty_;
}

} // end namespace AbstractLinAlgPack

#endif  // ALAP_VECTOR_SPACE_FACTORY_Thyra_HPP
