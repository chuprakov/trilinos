// //////////////////////////////////////////////////////////////////////
// ConstrainedOptPack_QPKWIK_Output.hpp
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
// These are subroutines that are called by QPKWIK to output data for
// debugging purposes.

#include <ostream>

namespace QPKWIK_Output {

/// 
/** Class for setting the output stream (destructor unsets it).
  *
  * Create an object of this type to set a stream for outputing
  * before you call QPKWIK.
  * This is not tread safe.
  */
class set_output {
public:
	///
	set_output(std::ostream* out);
	///
	~set_output();
private:
	// not defined and not to be called
	set_output();
	set_output(const set_output&);
	set_output& operator=(const set_output&);
};	// end class set_output

// Output stream to use (default == 0, no output).
extern std::ostream* out;

}	// end namespace QPKWIK_Output
