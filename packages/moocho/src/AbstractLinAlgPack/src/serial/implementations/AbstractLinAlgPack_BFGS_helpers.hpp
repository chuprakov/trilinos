// ////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_BFGS_helpers.hpp
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

#ifndef BFGS_HELPERS_H
#define BFGS_HELPERS_H

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

///
/** @name Functions to be used in BFGS updating.
 */
//@{

///
/** Check that s'*y is sufficiently positive and print the result if it is not.
 *
 * @param  s           [in] DVector (size n): Secant update vector B*s=y.
 * @param  y           [in] DVector (size n): Secant update vector for B*s=y.
 * @param  sTy         [in] If sTy != NULL then *sTy must contain the value computed
 *                     from dot(s,y).  If sTy == NULL, then this value will be computed
 *                     internally.  This argument is included so that the same computation
 *                     does not have to be performed more than once by the client and
 *                     this function.
 * @param  out         [in/out] If out==NULL then no output will be printed if the
 *                     condition fails.  If out!=NULL and the function returns false
 *                     then a small amount of output will be sent to *out.
 * @param  func_name   [in] The name of the function this is being called from.
 *                     If the condition is not met then this name in included
 *                     in the printout to *out.  If func_name == NULL then this
 *                     string will obviously not be printed.
 *
 * @return If s'*y >= sqrt(mach_epsilon) * ||s||2 * ||y||2 then this function will return true.
 * Otherwise it will return false.
 */
bool BFGS_sTy_suff_p_d(
	const Vector    &s
	,const Vector   &y
	,const value_type     *sTy        = NULL
	,std::ostream         *out        = NULL
	,const char           func_name[] = NULL
	);

//@}

} // end namespace AbstractLinAlgPack

#endif // BFGS_HELPERS_H
