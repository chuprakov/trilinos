/*/////////////////////////////////////////////////////
// RTOp_config.h
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
// If the macro RTOp_USE_MPI is defined, then these
// declarations will be MPI compatible.  If not then
// dummy MPI declarations will be used.
//
*/

#ifndef REDUCT_TRANS_VECTOR_OPERATORS_CONFIG_H
#define REDUCT_TRANS_VECTOR_OPERATORS_CONFIG_H

/* KL: Read the configuration header
// RAB: 2003/3/5
*/
#ifndef MOOCHO_NO_TRILINOS
#include "TSFCoreUtils_ConfigDefs.hpp"
#endif

/* KL: Test value of HAVE_MPI, conforming to autotools conventions for
// conditional compilation.
*/
#ifdef HAVE_MPI
#include "mpi.h"       /* Use real MPI declarations */
#else
#include "RTOp_mpi.h"  /* Use dummy MPI declarations */
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_config.h Platform dependent configuration options for RTOp.
 *
 * These typedefs and macros can be adjusted to the specific requirements of the platform.
 * For example, <tt>long double</tt> could be used instead of \c double if greater
 * precsion is needed.
 *
 * Also included are a few macros for MPI interoperability.  Also included in this default
 * header file is is \c RTOp_mpi.h which contains dummy MPI declarations (which are defined
 * in \c RTOp_mpi.c) for a subset of the MPI functions that are correct for the serial
 * case.
 */
/*@{*/

typedef MPI_Datatype            RTOp_Datatype;   /*!< Compatible with MPI_Datatype? */
typedef double                  RTOp_value_type; /*!< Compatible with fortran DOUBLE PRECISION? */
typedef int                     RTOp_index_type; /*!< Must be a signed integer!  Compatible with fortran INTEGER? */
typedef char                    RTOp_char_type;  /*!< Compatible with fortran CHARACTER? */
#define RTOpMPI_VALUE_TYPE      MPI_DOUBLE       /*!< (MPI only) Compatible with fortran DOUBLE PRECISION? */
#define RTOpMPI_INDEX_TYPE      MPI_INT          /*!< (MPI only) Compatible with fortran INTEGER? */
#define RTOpMPI_CHAR_TYPE       MPI_CHAR         /*!< (MPI only) Compatible with fortran CHARACTER? */

/* The maxinum number of characters in a name of a reduction/transformation operator class */
#define RTOp_MAX_REDUCT_TRANS_OP_CLASS_NAME  50

/*@}*/

#ifdef __cplusplus
}
#endif

#endif /* REDUCT_TRANS_VECTOR_OPERATORS_CONFIG_H */
