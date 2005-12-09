// //////////////////////////////////////////////////////////////////////////////////
// ConstrainedOptPack_MeritFuncPenaltyParam.hpp
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

#ifndef MERIT_FUNC_PENALTY_PARAM_H
#define MERIT_FUNC_PENALTY_PARAM_H

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

///
/** This class provides interface for setting and retrieving a penalty parameter
  * that many merit functions use {abstract}.
  */
class MeritFuncPenaltyParam {
public:

	///
	virtual ~MeritFuncPenaltyParam() {}

	/// Set the penalty parameter mu
	virtual void mu(value_type mu) = 0;

	/// Get the value of mu
	virtual value_type mu() const = 0;

};	// end class MeritFuncPenaltyParam

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_PENALTY_PARAM_H
