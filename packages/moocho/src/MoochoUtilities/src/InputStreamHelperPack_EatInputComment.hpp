// //////////////////////////////////////////////////////////////////////////////
// InputStreamHelperPack_EatInputComment.hpp
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

#ifndef EAT_INPUT_COMMENT_H
#define EAT_INPUT_COMMENT_H

#include "Moocho_ConfigDefs.hpp"

namespace InputStreamHelperPack {

/** @name namespace InputStreamHelperPack
 *
 * @memo Basic input stream helper functions.
 */
//@{

///
/** Discards comment lines from an input stream.
 *
 * Call this function to discard text that includes comment lines and
 * lines with only newline chars.  Here comment lines are in the tradition
 * of Fortran in that the comment line must begin with the comment
 * identification character (user supplied) and ends with a newline char
 * ('\n').
 *
 * In particular this function will discard any white
 * space before a set of consecutive comment lines and
 * the comment lines and will leave the input steam
 * at the first non-comment line that does not start with '\n'.
 */
std::istream& eat_comment_lines(std::istream& is, char comment_identifier);

//@}

}	// end namespace InputStreamHelperPack

#endif	// EAT_INPUT_COMMENT_H
