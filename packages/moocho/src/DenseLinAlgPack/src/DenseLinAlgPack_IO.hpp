// //////////////////////////////////////////////////////////////////////////////
// DenseLinAlgPackIO.hpp
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
// Includes all the headers for the DenseLinAlgPack input/output functions and
// stream operators for DVector, DVectorSlice, DMatrix, and DMatrixSlice.
// The utility function for "eating" comment lines is also included.

#ifndef LINALGPACK_IO_H
#define LINALGPACK_IO_H

#include "EatInputComment.hpp"
#include "DVectorIn.hpp"
#include "DVectorOut.hpp"
#include "DMatrixIn.hpp"
#include "DMatrixOut.hpp"
#include "DenseLinAlgPackInFormatDef.hpp"
#include "DenseLinAlgPackOutFormatDef.hpp"

// Include namelookups for templated operator functions.  MS VS++ 5.0 standard
// nonconformance problem.

#include "DenseLinAlgPackIO_NameLookups.hpp"

#endif // LINALGPACK_IO_H
