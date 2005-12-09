// ////////////////////////////////////////////////////////////
// AbstractLinAlgPack_BasisSystemFactoryStd.hpp
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

#ifndef SPARSE_SOLVER_PACK_BASIS_SYSTEM_FACTORY_H
#define SPARSE_SOLVER_PACK_BASIS_SYSTEM_FACTORY_H

#include "AbstractLinAlgPack_Types.hpp"
#include "AbstractLinAlgPack_BasisSystemFactory.hpp"
#include "Teuchos_AbstractFactory.hpp"

namespace AbstractLinAlgPack {

///
/** Default implementation for <tt>BasisSystemPermDirectSparse</tt> obejcts
 * using <tt>DirectSparseSolver</tt> object.
 *
 * Several direct sparse solvers are supported by default.  These include:
 * <ul>
 * <li> DENSE (using LAPACK xGETRF())
 * <li> MA28
 * <li> MA48 (using MA28 for BasisSystemPerm::select_basis()) (not yet)
 * <li> SuperLU
 * </ul>
 *
 * These solvers are supported only if the proper macros are defined.
 * 
 * ToDo: Create a DirectSparseSolverFactory interface and use this
 * to allow clients to add new DirectSparseSolvers ...
 *
 */
class BasisSystemFactoryStd
	: public AbstractLinAlgPack::BasisSystemFactory
{
public:

	///
	BasisSystemFactoryStd(); // ToDo: Add arguments!

	/** @name Overridden from BasisSystemFactory */
	//@{

	///
	void set_options( const options_ptr_t& options );
	///
	const options_ptr_t& get_options() const;

	//@}

	/** @name Overridden from AbstractFactory */
	//@{

	///
	obj_ptr_t create() const;

	//@}

private:

	// ////////////////////////
	// Private types

	enum EDirectLinearSolverType { LA_DENSE, LA_MA28, LA_MA48, LA_SUPERLU };

	// ////////////////////////
	// Private data members

	mutable EDirectLinearSolverType  direct_linear_solver_type_;
	options_ptr_t                    options_;

	// ////////////////////////
	// Private member functions
	
	void read_options() const;

}; // end class BasisSystemFactoryStd

}  // end namespace AbstractLinAlgPack

#endif // SPARSE_SOLVER_PACK_BASIS_SYSTEM_FACTORY_H
