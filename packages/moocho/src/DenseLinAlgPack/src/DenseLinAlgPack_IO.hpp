// //////////////////////////////////////////////////////////////////////////////
// LinAlgPackIO.h
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