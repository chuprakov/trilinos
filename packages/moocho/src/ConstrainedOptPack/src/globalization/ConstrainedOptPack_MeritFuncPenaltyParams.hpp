// //////////////////////////////////////////////////////////////////////////////////
// ConstrainedOptPack_MeritFuncPenaltyParams.hpp
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

#ifndef MERIT_FUNC_PENALTY_PARAMS_H
#define MERIT_FUNC_PENALTY_PARAMS_H

#include "ConstrainedOptPack_Types.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"

namespace ConstrainedOptPack {

///
/** This class provides interface for setting and retrieving a penalty parameter
  * that many merit functions use {abstract}.
  */
class MeritFuncPenaltyParams {
public:

	///
	class CanNotResize : public std::logic_error
	{public: CanNotResize(const std::string& what_arg) : std::logic_error(what_arg) {}};

	///
	virtual ~MeritFuncPenaltyParams() {}

	/** @name To be overridden by subclasses */
	//@{

	///
	/** Set the vector space for \c to use for the penalty parameters.
	  */
	virtual void set_space_c( const VectorSpace::space_ptr_t& space_c ) = 0;

	/// Get the vector of penalty parameters for setting them
	virtual VectorMutable& set_mu() = 0;

	/// Get the vector of penalty parameters for viewing them
	virtual const Vector& get_mu() const = 0;

	//@}

};	// end class MeritFuncPenaltyParams

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_PENALTY_PARAMS_H
