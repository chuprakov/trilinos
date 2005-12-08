// ///////////////////////////////////////////////////////
// DenseLinAlgLAPack_Types.hpp
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

#ifndef LIN_ALG_LA_PACK_TYPES_H
#define LIN_ALG_LA_PACK_TYPES_H

#include <stdexcept>

#include "DenseLinAlgPack_Types.hpp"

namespace DenseLinAlgLAPack {

// Include public types from DenseLinAlgPack
#include "DenseLinAlgPack_PublicTypes.ud"

/// Exception for factorization error
class FactorizationException : public std::logic_error
{public: FactorizationException(const std::string& what_arg) : std::logic_error(what_arg) {}};

}	// end namespace DenseLinAlgLAPack

#endif	// LIN_ALG_LA_PACK_TYPES_H
