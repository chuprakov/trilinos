# @HEADER
# ************************************************************************
# 
# Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
#
#   Copyright (2003) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
# @HEADER


# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#
# There is, at present, no distinction between most of the categories, 
# other than those which start with EXTRA_.  All categories which do not
# start with EXTRA_ are included in all builds.  Categories which start
# with EXTRA_ are included in some builds (typically controlled by 
# --with-xyz or --without-xyz)
#

MOOCHO_UTILITIES_HPP = \
  BLAS_Cpp_Types.hpp\
  FortranTypes_CppFortranStrings.hpp\
  FortranTypes_f_open_file.hpp\
  InputStreamHelperPack_EatInputComment.hpp\
  OptionsFromStreamPack_OptionsFromStreamExceptions.hpp\
  OptionsFromStreamPack_OptionsFromStream.hpp\
  OptionsFromStreamPack_SetOptionsFromStream.hpp\
  OptionsFromStreamPack_SetOptionsFromStreamNode.hpp\
  OptionsFromStreamPack_SetOptionsToTargetBase.hpp\
  OptionsFromStreamPack_StringToBool.hpp\
  OptionsFromStreamPack_StringToIntMap.hpp\
  ProfileHackPack_profile_hack.hpp\
  SerializationPack_Serializable.hpp\
  SpecialASCII.hpp\
  StopWatchPack_stopwatch.hpp\
  StrideIterPack_StrideIter.hpp\
  TestingHelperPack_update_success.hpp

MOOCHO_UTILITIES_CPP = \
  FortranTypes_CppFortranStrings.cpp\
  FortranTypes_f_open_file.cpp\
  InputStreamHelperPack_EatInputComment.cpp\
  OptionsFromStreamPack_OptionsFromStream.cpp\
  OptionsFromStreamPack_SetOptionsFromStreamNode.cpp\
  OptionsFromStreamPack_StringToBool.cpp\
  OptionsFromStreamPack_StringToIntMap.cpp\
  ProfileHackPack_profile_hack.cpp\
  StopWatchPack_stopwatch.cpp\
  TestingHelperPack_update_success.cpp

MOOCHO_UTILITIES_F = \
  ff_open_file.f\
  FortranCppStrings.f\
  FSeconds.f

#
# Specify locations of include files within Trilinos
#

include $(top_builddir)/Makefile.export.moocho

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(MOOCHO_INCLUDES))
else
EXPORT_INCLUDES = $(MOOCHO_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# --------------------------------------------------------------------
# Moocho library specifications
# --------------------------------------------------------------------

MOOCHO_LIB = libmoocho.a

MOOCHO_H = \
  $(srcdir)/Moocho_ConfigDefs.hpp \
	$(MOOCHO_UTILITIES_HPP)

libmoocho_a_SOURCES = \
	$(MOOCHO_UTILITIES_CPP) \
  $(MOOCHO_UTILITIES_F)

# --------------------------------------------------------------------
# For using a special archiver
# --------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libmoocho_a_AR = $(ALTERNATE_AR)

else

libmoocho_a_AR = $(AR) cru

endif

# --------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# --------------------------------------------------------------------

# I have not been able to find any kind of "or" for Automake
# conditionals. JW

 libmoocho_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# --------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# --------------------------------------------------------------------

lib_LIBRARIES = $(MOOCHO_LIB)

# --------------------------------------------------------------------
# List of all headers to install in $(includedir)
# --------------------------------------------------------------------

include_HEADERS = $(MOOCHO_H) 

# --------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# --------------------------------------------------------------------

# SPECIAL NOTE: MoochoConfig.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = Moocho_Config.h $(MOOCHO_H)

# --------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# --------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Library Target for Circular Dependencies
# ------------------------------------------------------------------------
all-libs:
	$(MAKE) libmoocho.a
	cd $(top_builddir)/../teuchos ; $(MAKE) all-libs
