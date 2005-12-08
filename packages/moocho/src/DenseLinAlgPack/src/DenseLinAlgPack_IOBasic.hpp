// ///////////////////////////////////////////////////////////////////////////////////
// DenseLinAlgPack_IOBasic.hpp
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
// Basic declarations for I/O

#ifndef LINALGPACK_IO_BASIC_H
#define LINALGPACK_IO_BASIC_H

#include <ios>
#include <stdexcept>

#include "DenseLinAlgPack_Types.hpp"

namespace DenseLinAlgPack {
namespace LinAlgPackIO {

/// Exception throw on input error
class InputException : public std::logic_error
{public: InputException(const std::string& what_arg) : std::logic_error(what_arg) {}};

///
typedef std::ios_base::fmtflags fmtflags;	
	
/// Format flags
enum { ignore_dim_bit = 0x0001, no_insert_newlines_bit = 0x0002 };

}	// end namespace LinAlgPackIO
}	// end namespace DenseLinAlgPack



#endif	// LINALGPACK_IO_BASIC_H
