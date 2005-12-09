// //////////////////////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_SpVectorOut.hpp
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
// Output stream operator for SpVector

#ifndef SP_VECTOR_OUT_H
#define SP_VECTOR_OUT_H

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

///
/** SpVectorSlice output stream operator.
  */
std::ostream& operator<<(std::ostream& os, const SpVectorSlice& svs);

}	// end namespace AbstractLinAlgPack

#endif // SP_VECTOR_OUT_H
