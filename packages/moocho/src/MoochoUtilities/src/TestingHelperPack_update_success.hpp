// //////////////////////////////////////////////////////////////////////
// update_success.
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

#ifndef UPDATE_SUCCESS_H
#define UPDATE_SUCCESS_H

#include "Moocho_ConfigDefs.hpp"

namespace TestingHelperPack {

/** @name namespace TestingHelperPack
  *
  * @memo Helper functions for writting general testing code.
  */
//@{

/// Helper function for updating a flag for if an operation returned false
bool update_success(bool result_check, bool* success);

///
/** Global flag for setting if and exception is thrown if a result_check
  * == false causes an exception to be thrown.
  */
extern bool throw_except_on_fail;

//	end namespace TestingHelperPack
//@}

}	// end namespace TestingHelperPack

#endif	// UPDATE_SUCCESS_H
