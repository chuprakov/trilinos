// /////////////////////////////////////////////////////////////////////////
// OptionsFromStreamPack_SetOptionsFromStream.hpp
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

#ifndef SET_OPTIONS_FROM_STREAM_H
#define SET_OPTIONS_FROM_STREAM_H

#include "Moocho_ConfigDefs.hpp"

namespace OptionsFromStreamPack {

class OptionsFromStream;

///
/** Abstact interface for objects that have options to be set
  * that are contained in an OptionsFromStreamObject.
  *
  * \ingroup OptionsFromStreamPack_grp
  */
class SetOptionsFromStream {
public:
	
	///
	virtual ~SetOptionsFromStream() {}
	
	/// Call to set options from a stream
	virtual void set_options( const OptionsFromStream& options ) = 0;

};	// end class SetOptionsFromStream

}	// end namespace OptionsFromStreamPack

#endif	// SET_OPTIONS_FROM_STREAM_H
