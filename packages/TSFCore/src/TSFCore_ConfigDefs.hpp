/*
// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#ifndef _TSFCORE_CONFIGDEFS_H_
#define _TSFCORE_CONFIGDEFS_H_

#include <TSFCoreUtils_ConfigDefs.hpp>

#ifndef __cplusplus
#define __cplusplus
#endif

/*
 * The macros PACKAGE, PACKAGE_NAME, etc, get defined for each package and need to
 * be undef'd here to avoid warnings when this file is included from another package.
 * KL 11/25/02
 */
#ifdef PACKAGE
#undef PACKAGE
#endif

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif

#ifdef VERSION
#undef VERSION
#endif

#include <TSFCoreConfig.h>

#ifdef HAVE_OSTREAM
#include <ostream>
#elif defined(HAVE_IOSTREAM)
#include <iostream>
#elif defined(HAVE_IOSTREAM_H)
#include <iostream.h>
#endif

#ifdef HAVE_IOMANIP
#include <iomanip>
#else
#include <iomanip.h>
#endif

#ifdef HAVE_VALARRAY
#include <valarray>
#endif

#ifdef HAVE_STDEXCEPT
#include <stdexcept>
#endif

#ifdef HAVE_TYPEINFO
#include <typeinfo>
#endif

#ifdef HAVE_CASSERT
#include <cassert>
#else
#include <assert.h>
#endif

#ifdef HAVE_VECTOR
#include <vector>
#endif

#ifdef HAVE_ALGORITHM
#include <algorithm>
#elif defined(HAVE_ALGO_H)
#include <algo.h>
#else
#include <algorithm.h>
#endif

#ifdef HAVE_MEMORY
#include <memory>
#else
#include <memory.h>
#endif

#ifdef HAVE_STRING
#include <string>
#else
#include <string.h>
#endif

#ifdef HAVE_IOSFWD
#include <iosfwd>
#endif


/* Every line that begins with 'using' should eventually be dependent
   on some check within the configure script */


#ifndef TFLOP
#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif
#else /* TFLOP defined */
#ifdef HAVE_IOMANIP
#include <iomanip>
#else
#include <iomanip.h>
#endif
#ifdef HAVE_STRING
using std::string;
#endif
#ifdef HAVE_IOSTREAM
using std::istream;
using std::ostream;
using std::cerr;
using std::cout;
using std::endl;
#endif
#endif

#endif /*_TSFCORE_CONFIGDEFS_H_*/
