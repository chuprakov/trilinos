# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------
#
# Adapted from Tammy Kolda's configure system for NOX.
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------

# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.

AC_INIT(TSF, 1.04)

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.

AC_CONFIG_SRCDIR(src/Support/TSFBlas.h)

# --------- Create configuration header 

AC_CONFIG_HEADER(include/TSFConfig.h:Build/TSFConfig.h.in)

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)


# Specify required version of autoconf.

AC_PREREQ(2.52)

# Specify revision number of configure.ac that was used to generate
# configure.

dnl Commented out this line for now because it causes headaches in
dnl committing.
dnl AC_REVISION($Revision$)

# specify architechture

AC_ARG_WITH(arch,
[  --with-arch            Choose arch (options are LINUX|SunWS5.0)],
[ ARCH=$withval],
[AC_MSG_ERROR(must set --with-arch to one of LINUX|SunWS5.0)]
)

# ------------------------------------------------------------------------
# Check for user overrides (to be used later)
# ------------------------------------------------------------------------

if test -n "${CXXFLAGS}"; then
  USERDEF_CXXFLAGS=yes
fi

if test "$srcdir" = "."; then
	PROJECT_ROOT=`pwd`
else
	PROJECT_ROOT=$srcdir
fi

# ----------- choose whether to compile with debugging -----------------

AC_ARG_ENABLE(opt,
[  --enable-opt            enable optimized compilation],
[  
AC_MSG_CHECKING(whether we are compiling optimized)
AC_MSG_RESULT([${enableval}])
 OPT=${enableval}
],
[
AC_MSG_CHECKING(whether we are compiling optimized)
AC_MSG_RESULT([no])
OPT=no]
)

AC_ARG_ENABLE(debug,
[  --enable-debug            produce debugging information],
[  
AC_MSG_CHECKING(whether we are producing debugging information)
AC_MSG_RESULT([${enableval}])
 DEBUG=${enableval}
],
[
AC_MSG_CHECKING(whether we are producing debugging information)
AC_MSG_RESULT([yes])
DEBUG=yes]
)



AC_ARG_ENABLE(boundscheck,
[  --enable-boundscheck            enable dense vector bounds checking],
[
	AC_MSG_CHECKING(whether to do dense vector bounds checking)
	AC_MSG_RESULT([${enableval}])
	case $enableval in
  yes)
 	AC_DEFINE(HAVE_DENSE_VECTOR_BOUNDSCHECK, 1,[enable dense vector bounds checks])
  ;;
	no)
 	AC_DEFINE(HAVE_DENSE_VECTOR_BOUNDSCHECK, 0,[disable dense vector bounds checks])
  ;;
esac
],
[
	AC_MSG_CHECKING(whether to do dense vector bounds checking)
	AC_MSG_RESULT([yes])
	AC_DEFINE(HAVE_DENSE_VECTOR_BOUNDSCHECK, 1,[enable dense vector bounds checks])],
)


# ------------------------------------------------------------------------
# Math libraries
# ------------------------------------------------------------------------


AC_ARG_WITH(math,
[  --with-math              specify math libraries (default -llapack -lblas -lm)],
[
MATH_LIB=${withval}
AC_MSG_CHECKING(using alternate math libraries)
AC_MSG_RESULT([${MATH_LIB}])
],
[
MATH_LIB="-llapack -lblas"
AC_MSG_CHECKING(using math libraries)
AC_MSG_RESULT([${MATH_LIB}])
]
)

AC_ARG_WITH(mathdir,
[  --with-mathdir              specify path to math libraries],
[
MATH_LIBDIR=${withval}
AC_MSG_CHECKING(path to math libraries)
AC_MSG_RESULT([${MATH_LIBDIR}])
]
)



# ------------------------------------------------------------------------
# MPI configure command line option
# ------------------------------------------------------------------------

AC_ARG_ENABLE(mpi,
[  --enable-mpi            enable MPI using mpiCC C++ compiler],
[
HAVE_PKG_MPI=yes 
MPI_CXX=mpiCC
],
[HAVE_PKG_MPI=no]
)

AC_ARG_WITH(mpi,
[  --with-mpi              specify root directory of MPI installation],
[
HAVE_PKG_MPI=yes
MPI_DIR=${withval}
AC_MSG_CHECKING(MPI directory)
AC_MSG_RESULT([${MPI_DIR}])
]
)

AC_ARG_WITH(mpi-cxx,
[  --with-mpi-cxx          specify MPI C++ compiler],
[
HAVE_PKG_MPI=yes 
MPI_CXX=${withval}
AC_MSG_CHECKING(user-defined MPI C++ compiler)
AC_MSG_RESULT([${MPI_CXX}])
]
)

AC_ARG_WITH(mpi-include,
[  --with-mpi-include      specify include directory for MPI],
[
HAVE_PKG_MPI=yes
MPI_INC=${withval}
AC_MSG_CHECKING(user-defined MPI includes)
AC_MSG_RESULT([${MPI_INC}])
]
)

AC_ARG_WITH(mpi-libs,
[  --with-mpi-libs         specify MPI libraries],
[
HAVE_PKG_MPI=yes
MPI_LIBS=${withval}
AC_MSG_CHECKING(user-defined MPI libraries)
AC_MSG_RESULT([${MPI_LIBS}])
]
)

AC_ARG_WITH(mpi-libdir,
[  --with-mpi-libdir       specify location of MPI libraries],
[
HAVE_PKG_MPI=yes
MPI_LIBDIR=${withval}
AC_MSG_CHECKING(user-defined MPI libraries)
AC_MSG_RESULT([${MPI_LIBS}])
]
)

AC_MSG_CHECKING(whether we are using MPI)
AC_MSG_RESULT([${HAVE_PKG_MPI}])

# ------------------------------------------------------------------------
# Check for MPI compiler (must be done *before* AC_PROG_CXX)
# ------------------------------------------------------------------------



if test -n "${MPI_CXX}"; then
  AC_CHECK_PROG(MPI_CXX_EXISTS, ${MPI_CXX}, yes, no)
  if test "X${MPI_CXX_EXISTS}" = "Xyes"; then
    CXX=${MPI_CXX}
  else
    AC_MSG_ERROR([MPI C++ compiler (${MPI_CXX}) not found.])
  fi
fi

# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CXX(CC c++ g++ gcc cxx)
AC_PROG_CC(cc gcc)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Checks for extra flags, etc
# ------------------------------------------------------------------------

AC_ARG_WITH(flags,
[  --with-flags            add extra preprocessor flags ("-Dflag -Idir")],
[XTRAFLAGS=${withval}]
)

AC_MSG_CHECKING([for extra preprocessor flags])
if test -z "${XTRAFLAGS}"; then
  AC_MSG_RESULT(none)
else
  AC_MSG_RESULT([${XTRAFLAGS}])
  CPPFLAGS="${CPPFLAGS} ${XTRAFLAGS}"
fi


AC_MSG_CHECKING([for special C++ flags])
if test "X${USERDEF_CXXFLAGS}" = "Xyes"; then

  AC_MSG_RESULT([${CXXFLAGS} (user-specified via CXXFLAGS)])

else

  if test -n "${SPECIAL_CXXFLAGS}"; then
    AC_MSG_RESULT([${SPECIAL_CXXFLAGS}])
    CXXFLAGS="${CXXFLAGS} ${CXXOPTIONS}"
  else
    AC_MSG_RESULT([${CXXFLAGS} (default)])
  fi

fi




# ------------------------------------------------------------------------
# MPI link check
# ------------------------------------------------------------------------

if test "X${HAVE_PKG_MPI}" = "Xyes"; then

	AC_DEFINE(HAVE_MPI, 1)

  if test -n "${MPI_DIR}" && test -z "${MPI_INC}"; then
    MPI_INC=${MPI_DIR}/include
  fi

  if test -n "${MPI_DIR}" && test -z "${MPI_LIBDIR}"; then
    MPI_LIBDIR=${MPI_DIR}/lib
  fi

  if test -n "${MPI_INC}"; then
    CPPFLAGS="${CPPFLAGS} ${MPI_INC}"
  fi

  if test -n "${MPI_LIBDIR}"; then
    LDFLAGS="${LDFLAGS} ${MPI_LIBDIR}"
  fi

  if test -z "${MPI_LIBS}" && test -n "${MPI_LIBDIR}"; then
    MPI_LIBS="-lmpi"
  fi

  if test -n "${MPI_LIBS}"; then
    LIBS = "${MPI_LIBS} ${LIBS}"
  fi

  AC_LANG_CPLUSPLUS 
  AC_MSG_CHECKING(whether MPI will link using C++ compiler)
  AC_TRY_LINK([#include <mpi.h>],
  [int c; char** v; MPI_Init(&c,&v);],
  [AC_MSG_RESULT(yes)], 
  [AC_MSG_RESULT(no)  
   AC_MSG_ERROR(MPI cannot link)]
  )

else

	AC_DEFINE(HAVE_MPI, 0)

fi



# ------------------------------------------------------------------------
# Checks for header file style: <header> or <header.h>
# ------------------------------------------------------------------------

#iostream
AC_LANG_CPLUSPLUS 
AC_MSG_CHECKING([for new header style])
AC_TRY_LINK([#include <iostream>],[],
[
  AC_MSG_RESULT(yes)
  AC_DEFINE(HAVE_NEW_HEADERS,1)
])


#stringstream flavor: strstream or sstream. If strstream exists, use it,
# because gnu 2.95's sstream is botched. 
AC_LANG_CPLUSPLUS 
AC_MSG_CHECKING([for strstream])
AC_TRY_LINK([#include <strstream>],[],
[
  AC_MSG_RESULT(yes)
  AC_DEFINE(HAVE_STRSTREAM,1)
])

AC_LANG_CPLUSPLUS 
AC_MSG_CHECKING([for strstream.h])
AC_TRY_LINK([#include <strstream.h>],[],
[
  AC_MSG_RESULT(yes)
  AC_DEFINE(HAVE_STRSTREAM,1)
])


# can we say using std::ostream?
AC_LANG_CPLUSPLUS 
AC_MSG_CHECKING(whether "using std::ostream" is okay)
AC_TRY_LINK(
[
#ifdef HAVE_IOSTREAM
#include <iostream>
#elif defined(HAVE_IOSTREAM_H)
#include <iostream.h>
#endif
],
[using std::ostream;],
[AC_MSG_RESULT(yes)], 
[
  AC_MSG_RESULT(no)  
  AC_DEFINE(NO_USING_STD_OSTREAM,1)
]
)


# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

# ------------------------------------------------------------------------
# Specify package directory
# ------------------------------------------------------------------------

AC_ARG_WITH(pkgdir,
[  --with-pkgdir           specify root directory of add'l packages],
[
USE_PKGDIR=yes
PKGDIR=${withval}
]
)

if test "X${USE_PKGDIR}" = "Xyes"; then

   AC_MSG_CHECKING(whether package directory is valid)
   if test -d ${PKGDIR}; then
      AC_MSG_RESULT(yes)
   else
      AC_MSG_NOTICE(no)
      AC_MSG_ERROR([Specified package directory does not exist])
   fi

   LDFLAGS="${LDFLAGS} -L${PKGDIR}/lib"
   CPPFLAGS="${CPPFLAGS} -I${PKGDIR}/include"

fi

# --------- Specify additional include paths -----------

AC_ARG_WITH(include,
[  --with-include          Specifies include directories (\"-Idir1 -Idir2 ...\")],
[
XTRAINC=$withval
echo Using extra include directories: $withval
]
)

# --------- Specify additional libraries -----------

AC_ARG_WITH(libs,
[  --with-libs             Specifies libraries (\"-Ldir -llib ...\")],
[
XTRALIBS=$withval
echo Using extra libraries: $withval
]
)

# ------------------------------------------------------------------------
# Checks for types
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for structures
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

# THIS TEST ISN'T WORKING - WILL COME BACK TO IT LATER.
# Check that all libraries will link with C++
dnl if test X${LIBS} != X; then
dnl   AC_MSG_CHECKING([[whether libraries (${LIBS}) will link]])
dnl   AC_LANG_CPLUSPLUS 
dnl   AC_TRY_LINK([],[],
dnl   [AC_MSG_RESULT(yes)],
dnl   [
dnl   AC_MSG_RESULT(no)
dnl   AC_MSG_ERROR(cannot link)
dnl   ]
dnl   )
dnl fi

# ------------------------------------------------------------------------
# Checks for compiler characteristics
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

# -- DEC CXX --
dnl if test X${target_cpu} = Xalpha && test X${CXX} = Xcxx; then
dnl    XTRALDADD="cxx_repository/*.o"
dnl fi

dnl AC_MSG_CHECKING([for extra object files to be linked])
dnl if test X$XTRALDADD = X; then
dnl AC_MSG_RESULT(none)
dnl else
dnl AC_MSG_RESULT($XTRALDADD)
dnl fi

AC_SUBST(XTRALDADD)

# ------------------------------------------------------------------------
# Checks for library functions
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for system services
# ------------------------------------------------------------------------


AC_F77_LIBRARY_LDFLAGS
AC_F77_DUMMY_MAIN


# ------------------------------------------------------------------------
# do we use Trilinos?
# ------------------------------------------------------------------------

AC_ARG_ENABLE(trilinos,
[  --disable-trilinos     disable hooks to trilinos],
[HAVE_PETRA=${enableval}],
[HAVE_PETRA=yes]
)

if test "X${HAVE_PETRA}" = "Xyes"; then

  # Set TRILINOS_COMM (used several places)
  if test "X${HAVE_PKG_MPI}" = "Xyes"; then
    TRILINOS_COMM="MPI"
  else
    TRILINOS_COMM="SERIAL"
  fi

	AC_MSG_CHECKING(whether TRILINOS_HOME is defined)
  
	if test -n "${TRILINOS_HOME}"; then
     AC_MSG_RESULT([yes (${TRILINOS_HOME})])
  else
       AC_MSG_RESULT([no (using ${HOME}/Trilinos)])
       TRILINOS_HOME=${HOME}/Trilinos
  fi
 
    AC_MSG_CHECKING(whether TRILINOS_HOME is valid)
    if test -d ${TRILINOS_HOME}; then
       AC_MSG_RESULT(yes)
    else
       AC_MSG_RESULT(no)
       AC_MSG_ERROR([TRILINOS_HOME is invalid])
    fi
 
    AC_MSG_CHECKING(whether TRILINOS_ARCH is defined)
    if test -n "${TRILINOS_ARCH}"; then
       AC_MSG_RESULT([yes (${TRILINOS_ARCH})])
    else
       AC_MSG_RESULT(no)
       AC_MSG_ERROR([TRILINOS_ARCH is invalid])
    fi

    AC_MSG_CHECKING(whether TRILINOS_ID is defined)
    if test -n "${TRILINOS_ID}"; then
       AC_MSG_RESULT([yes (${TRILINOS_ID})])
    else
       AC_MSG_RESULT(no)
    fi
 
    AC_MSG_CHECKING(whether TRILINOS_TARGET is defined)
    if test -n "${TRILINOS_TARGET}"; then
       AC_MSG_RESULT([yes (${TRILINOS_TARGET})])
    else
       TRILINOS_TARGET="${TRILINOS_ARCH}.${TRILINOS_COMM}${TRILINOS_ID}"
       AC_MSG_RESULT([no (using ${TRILINOS_TARGET})])
    fi

    EPETRA_INC=${TRILINOS_HOME}/packages/epetra/src
    AZTECOO_INC=${TRILINOS_HOME}/packages/aztecoo/src
    IFPACK_INC=${TRILINOS_HOME}/packages/ifpack/src
    TRILINOS_LIBDIR=-L${TRILINOS_HOME}/lib/${TRILINOS_TARGET}
    TRILINOS_LIBS="-laztecoo -lifpack -lepetra -ly12m -lspblas"

  # If we make it through all the tests, define the correct libs and includes!!
  TRILINOS_INC="-I${EPETRA_INC} -I${IFPACK_INC} -I${AZTECOO_INC}"

	AC_DEFINE(HAVE_PETRA, 1)
else
	AC_DEFINE(HAVE_PETRA, 0)	
fi


AC_SUBST(TRILINOS_INC)
AC_SUBST(TRILINOS_LIBDIR)
AC_SUBST(TRILINOS_LIBS)
AC_SUBST(HAVE_PETRA)

# ------------------------------------------------------------------------
# Perform substitutions in output files
# ------------------------------------------------------------------------

AC_SUBST(ARCH)
AC_SUBST(ac_aux_dir)
AC_SUBST(PROJECT_ROOT)
AC_SUBST(HAVE_PKG_EXPAT)
AC_SUBST(MPI_LIBS)
AC_SUBST(HAVE_PKG_MPI)
AC_SUBST(MPI_LIBDIR)
AC_SUBST(MPI_INC)
AC_SUBST(OPT)
AC_SUBST(DEBUG)
AC_SUBST(FLIBS)
AC_SUBST(MATH_LIB)
AC_SUBST(MATH_LIBDIR)



# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------


AC_OUTPUT(Build/Makefile.site:Build/Makefile.site.in Build/Makefile.export:Build/Makefile.export.in doc/Doxygen/Doc.cfg:doc/Doxygen/Doc.cfg.in)


