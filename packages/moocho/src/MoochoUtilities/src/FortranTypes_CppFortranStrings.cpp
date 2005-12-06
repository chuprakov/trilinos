// //////////////////////////////////////////////////////////
// CppFortranStrings.pp
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
//
// ToDo: 8/9/99: Make sure that these are valid ASCII conversions
// on each target platform.  For example, UNIX does not use
// all 8 bits for characters.

#include "FortranTypes_CppFortranStrings.hpp"

int FortranTypes::convert_to_f_int_string( const char string[], f_int i_string[]
	, f_int* string_len )
{
	*string_len = 0;
	for( ; *string != 0; ++(*string_len) )
		*i_string++ = *string++;	// Sse standard C conversions by default
	return 0;	// success
}

int FortranTypes::convert_from_f_int_string( const f_int i_string[], f_int string_len
	, char string[] )
{
	for( f_int i = 0; i < string_len; )
		*string++ = *i_string++;	// Sse standard C conversions by default
	*string = 0; // Null terminate the target string.
	return 0;	// success
}
