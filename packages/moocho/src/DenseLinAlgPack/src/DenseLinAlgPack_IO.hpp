// //////////////////////////////////////////////////////////////////////////////
// LinAlgPackIO.h
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
// Includes all the headers for the LinAlgPack input/output functions and
// stream operators for Vector, VectorSlice, GenMatrix, and GenMatrixSlice.
// The utility function for "eating" comment lines is also included.

#ifndef LINALGPACK_IO_H
#define LINALGPACK_IO_H

#include "Misc/include/EatInputComment.h"
#include "VectorIn.h"
#include "VectorOut.h"
#include "GenMatrixIn.h"
#include "GenMatrixOut.h"
#include "LinAlgPackInFormatDef.h"
#include "LinAlgPackOutFormatDef.h"

// Include namelookups for templated operator functions.  MS VS++ 5.0 standard
// nonconformance problem.

#include "LinAlgPackIO_NameLookups.h"

#endif // LINALGPACK_IO_H
