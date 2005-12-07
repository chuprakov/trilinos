// //////////////////////////////////////////////////////////////////
// OptionsFromStreamPack_OptionsFromStreamExceptions.hpp
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

#ifndef OPTIONS_FROM_STREAM_EXCEPTIONS_H
#define OPTIONS_FROM_STREAM_EXCEPTIONS_H

#include "Moocho_ConfigDefs.hpp"

namespace OptionsFromStreamPack {

/** \defgroup OptionsFromStreamExceptions Exception classes.
 * \ingroup OptionsFromStreamPack_grp
 */
//@{

/// Input from stream error
class InputException : public std::logic_error
{public: InputException(const std::string& what_arg) : std::logic_error(what_arg) {}};

/// Access error
class AccessException : public std::logic_error
{public: AccessException(const std::string& what_arg) : std::logic_error(what_arg) {}};

//@}

}	// end namespace OptionsFromStreamPack 

#endif	// OPTIONS_FROM_STREAM_EXCEPTIONS_H
