// ///////////////////////////////////////////
// AbstractLinAlgPack_GenPermMatrixSliceOut.hpp
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

#ifndef GEN_PERM_MATRIX_SLICE_OUT_H
#define GEN_PERM_MATRIX_SLICE_OUT_H

#include <iosfwd>

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

///
/** Print the contents of a GenPermMatrixSlice object in machine readable format.
 *
 * The format is:
 * 
 \begin{varbatim}
  rows cols nz
   row_i(1):col_j(1) ... row_i(nz):col_j(nz)
  \end{verbatim}
 */
std::ostream& operator<<( std::ostream& out, const GenPermMatrixSlice& P );

} // end namespace AbstractLinAlgPack

#endif  // GEN_PERM_MATRIX_SLICE_OUT_H
