#ifndef TSFDEFS_H
#define TSFDEFS_H

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

#include "TSFConfig.h"

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

/* define HAVE_MPI=1 if we are linking with MPI */

#define MPICH_SKIP_MPICXX


/* define HAVE_PETRA=1 to link with the Petra vector/matrix classes */

#define HAVE_PETRA 1


/* figure out if we need to set the PETRA_MPI flag to tell petra to use MPI */

#ifdef HAVE_MPI
#define HAVE_PETRA_MPI 1
#define EPETRA_MPI
#else
#define HAVE_PETRA_MPI 0
#endif

typedef double TSFReal;

/* define HAVE_EXPAT=1 to link with the expat XML parser */

#define HAVE_EXPAT 0


/* define flags indicating debug and opt level */

#define TSF_DEBUG_LEVEL 1

#define TSF_OPT_LEVEL 0



#define HAVE_DENSE_VECTOR_BOUNDSCHECK 0

#define HAVE_ARRAY_BOUNDSCHECK 0


/* ------------ specification of header file flavor --------------- */


#if HAVE_CMATH
#include <cmath>
#elif HAVE_MATH_H
#include <math.h>
#else
#error "Found neither cmath nor math.h"
#endif

#if HAVE_CSTDIO
#include <cstdio>
#elif HAVE_STDIO_H
#include <stdio.h>
#else
#error "Found neither cstdio nor stdio.h"
#endif


#if HAVE_CSTDLIB
#include <cstdlib>
#elif HAVE_STDLIB_H
#include <stdlib.h>
#else
#error "Found neither cstdlib nor stdlib.h"
#endif

#if HAVE_STRING
#include <string>
#elif HAVE_STRING_H
#include <string.h>
#else
#error "Found neither string nor string.h"
#endif

#if HAVE_IOSTREAM
#include <iostream>
#elif HAVE_IOSTREAM_H
#include <iostream.h>
#else
#error "Found neither iostream nor iostream.h"
#endif

#if HAVE_SSTREAM
#include <sstream>
typedef std::ostringstream TSFOStringStream;
#elif HAVE_SSTREAM_H
#include <sstream.h>
typedef ostringstream TSFOStringStream;
#elif HAVE_STRSTREAM
#include <strstream>
typedef std::ostrstream TSFOStringStream;
#elif HAVE_STRSTREAM_H
#include <strstream.h>
typedef ostrstream TSFOStringStream;
#else
#error "Found neither sstream, sstream.h, strstream.h, nor strstream"
#endif




#endif
