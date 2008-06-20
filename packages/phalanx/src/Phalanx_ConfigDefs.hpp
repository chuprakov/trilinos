
//@HEADER
/*
************************************************************************

                   New_Package Example Package
               Copyright (2004) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
//@HEADER

#ifndef PHALANX_CONFIGDEFS_H
#define PHALANX_CONFIGDEFS_H

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

#include <Phalanx_config.hpp>


#ifdef HAVE_ALGORITHM
#include <algorithm>
#elif HAVE_ALGORITHM_H
#include <algorithm.h>
#else
#include <algo.h>
#endif

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#ifdef HAVE_CSTDDEF
#include <cstddef>
#else
#include <stddef.h>
#endif

#ifdef HAVE_IOMANIP
#include <iomanip>
#else
#include <iomanip.h>
#endif

#ifdef HAVE_IOSTREAM
#include <iostream>
#else
#include <iostream.h>
#endif

#ifdef HAVE_ITERATOR
#include <iterator>
#else
#include <iterator.h>
#endif

#ifdef HAVE_MAP
#include <map>
#else
#include <map.h>
#endif

#ifdef HAVE_SSTREAM
#include <sstream>
#else
#include <sstream.h>
#endif

#ifdef HAVE_STRING
#include <string>
#else
#include <string.h>
#endif

#ifdef HAVE_TYPEINFO
#include <typeinfo>
#else
#include <typeinfo.h>
#endif

#ifdef HAVE_VECTOR
#include <vector>
#else
#include <vector.h>
#endif

#endif /* PHALANX_CONFIGDEFS_H */
