# @HEADER
# @HEADER


# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside meros. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#
# There is, at present, no distinction between most of the categories, 
# other than those which start with EXTRA_.  All categories which do not
# start with EXTRA_ are included in all builds.  Categories which start
# with EXTRA_ are included in some builds (typically controlled by 
# --with-xyz or --without-xyz)
#
CORE_H = \
          $(srcdir)/Aztec2TSF.h  $(srcdir)/ReadPetra.h  \
          $(srcdir)/TSFLinearOperator2EpetraRowMatrix.h \
          $(srcdir)/KayLoghinRightOperatorSource.h      \
          $(srcdir)/RightBlockNSOperatorSource.h        \
          $(srcdir)/Vbr2Petra.h                         \
          $(srcdir)/KayLoghinSchurFactory.h             \
          $(srcdir)/SchurFactoryBase.h                  \
          $(srcdir)/NSBlockPreconditionerFactory.h  $(srcdir)/SchurFactory.h

CORE = \
          $(srcdir)/Aztec2TSF.cpp  $(srcdir)/ReadPetra.cpp  \
          $(srcdir)/TSFLinearOperator2EpetraRowMatrix.cpp \
          $(srcdir)/KayLoghinRightOperatorSource.cpp      \
          $(srcdir)/RightBlockNSOperatorSource.cpp        \
          $(srcdir)/Vbr2Petra.cpp                         \
          $(srcdir)/KayLoghinSchurFactory.cpp             \
          $(srcdir)/SchurFactoryBase.cpp                  \
          $(srcdir)/NSBlockPreconditionerFactory.cpp      \
          $(srcdir)/SchurFactory.cpp

UTIL_H = 

UTIL = 

ENUM_H = 

UNIPROC_H = 

UNIPROC = 

SERIAL_H = 

SERIAL = 

#
#  Files which are included in select builds are typically have
#  names that start with EXTRA
#
#  EXTRA_MPI is control by --with-mpi and --without-mpi which is
#  set up by TAC_ARC_CONFIG_MPI 
#
EXTRA_MPI_H =

EXTRA_MPI = 

if HAVE_MPI 
MPI_H = \
	$(EXTRA_MPI_H)
MPI = \
	$(EXTRA_MPI)

else

MPI_H =
MPI =

endif


USERBASE_H = 


LANG_WRAPS_H = 

LANG_WRAPS = 

DIST_H = 

DIST = 

COMP_H =

COMP = 


AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src \
              -I$(top_srcdir)/../TSF/src/LinearOperators  \
              -I$(top_srcdir)/../TSF/src/Implementations  \
              -I$(top_srcdir)/../TSF/src/LinearProblem    \
              -I$(top_srcdir)/../TSF/src/Vectors          \
              -I$(top_srcdir)/../TSF/src/Support          \
              -I$(top_srcdir)/../TSF/src/Preconditioners  \
              -I$(top_srcdir)/../TSF/src/VectorSpaces     \
              -I$(top_srcdir)/../TSF/src/LinearSolvers    \
              -I$(top_srcdir)/../aztecoo/src              \
              -I$(top_builddir)/../TSF/src/LinearOperators  \
              -I$(top_builddir)/../TSF/src/Implementations  \
              -I$(top_builddir)/../TSF/src/LinearProblem    \
              -I$(top_builddir)/../TSF/src/Vectors          \
              -I$(top_builddir)/../TSF/src/Support          \
              -I$(top_builddir)/../TSF/src/Preconditioners  \
              -I$(top_builddir)/../TSF/src/VectorSpaces     \
              -I$(top_builddir)/../TSF/src/LinearSolvers    \
              -I$(top_builddir)/../aztecoo/src              

# ------------------------------------------------------------------------
# meros  library specifications
# ------------------------------------------------------------------------

MEROS_LIB = libmeros.a


MEROS_H = \
	$(CORE_H) \
	$(UTIL_H) \
	$(ENUM_H) \
	$(UNIPROC_H) \
	$(SERIAL_H) \
	$(USERBASE_H) \
	$(DIST_H) \
	$(COMP_H) \
	$(MPI_H)

libmeros_a_SOURCES = \
	$(CORE) \
	$(UTIL) \
	$(UNIPROC) \
	$(SERIAL) \
	$(DIST) \
	$(COMP) \
	$(MPI)


#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libbpack_a_AR = $(ALTERNATE_AR)
else

libmeros_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libmeros_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(MEROS_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(pkgincludedir)
# ------------------------------------------------------------------------

pkginclude_HEADERS = $(MEROS_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(pkgexecincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: meros_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

pkgexecincludedir = $(exec_prefix)/include/meros
nodist_pkgexecinclude_HEADERS = meros_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






