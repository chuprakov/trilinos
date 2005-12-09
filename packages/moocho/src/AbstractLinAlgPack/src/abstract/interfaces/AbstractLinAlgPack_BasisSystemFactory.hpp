// ////////////////////////////////////////////////////////////
// AbstractLinAlgPack_BasisSystemFactory.hpp
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

#ifndef ABSTRACT_LIN_ALG_PACK_BASIS_SYSTEM_FACTORY_H
#define ABSTRACT_LIN_ALG_PACK_BASIS_SYSTEM_FACTORY_H

#include "AbstractLinAlgPack_Types.hpp"
#include "Teuchos_AbstractFactory.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace OptionsFromStreamPack {
	class OptionsFromStream;
}

namespace AbstractLinAlgPack {

///
/** Interface for a factory object that will create <tt>BasisSystem</tt> objects.
 *
 * 
 */
class BasisSystemFactory : public Teuchos::AbstractFactory<BasisSystem>
{
public:

	/** @name Public types */
	//@{

	///
	typedef Teuchos::RefCountPtr<
		const OptionsFromStreamPack::OptionsFromStream>             options_ptr_t;

	//@}

	///
	virtual ~BasisSystemFactory() {}

	///
	/** Set the options that will be used to determine what basis system will be returned
	 * from <tt>this->create()</tt>.
	 *
	 * Note that it is allowed for the client to alter <tt>*options.get()</tt> after
	 * this method is called so <tt>this</tt> had better read the options inside of the
	 * <tt>this->create()</tt> method.
	 */
	virtual void set_options( const options_ptr_t& options ) = 0;

	///
	/** Get the <tt>OptionsFromStream</tt> object being used to extract the options from.
	 */
	virtual const options_ptr_t& get_options() const = 0;

}; // end class BasisSystemFactory

}  // end namespace AbstractLinAlgPack

#endif // ABSTRACT_LIN_ALG_PACK_BASIS_SYSTEM_FACTORY_H
