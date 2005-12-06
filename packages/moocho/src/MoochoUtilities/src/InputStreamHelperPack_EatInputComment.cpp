// //////////////////////////////////////////////////////////////////////////////
// InputStreamHelperPack_EatInputComment.cpp
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

#include <ctype.h>	// Change to <cctype> then it works properly
#include <string>
#include <istream>
#include <iomanip>

#include "InputStreamHelperPack_EatInputComment.hpp"

std::istream& InputStreamHelperPack::eat_comment_lines(std::istream& is
	, char comment_identifier)
{
	const int many = 1000;
	is >> std::ws;
	char c;
	while( ( c = is.peek() ) == comment_identifier || c == '\n' ) {
		is.ignore(many,'\n');
		is >> std::ws; // Eat whitespace up to the next line
	}
	return is;
}
