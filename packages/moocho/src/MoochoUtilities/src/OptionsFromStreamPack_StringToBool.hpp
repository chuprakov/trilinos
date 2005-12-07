// //////////////////////////////////////////////////////////////////////
// OptionsFromStreamPack_StringToBool.hpp
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

#ifndef STRING_TO_BOOL_H
#define STRING_TO_BOOL_H

#include "Moocho_ConfigDefs.hpp"

namespace OptionsFromStreamPack {

///
/** Convert a string "true" or "false" into bool #true# or #false#.
  *
  * If the input string is not "true" or "false" then the exception
  * "InputException" will be thrown and the message will include
  * the name of the option this value is for.
  */
bool StringToBool( const char* opt_name, const char* str );

}	// end namespace OptionsFromStreamPack 

#endif	// STRING_TO_BOOL_H
