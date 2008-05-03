
EXEEXT = .exe

#
# Actual programs
#

noinst_PROGRAMS = simpleEpetraExternCInterface

simpleEpetraExternCInterface_SOURCES = \
  simpleEpetraExternCInterface.c \
  dummy.cpp
simpleEpetraExternCInterface_DEPENDENCIES = $(common_dependencies)
simpleEpetraExternCInterface_LDADD = $(common_ldadd)

#
# Common stuff
#

common_dependencies = \
  $(top_builddir)/src/libctrilinos.a

include $(top_builddir)/Makefile.export.CTrilinos

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(CTRILINOS_INCLUDES))
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(CTRILINOS_LIBS))
else
EXPORT_INCLUDES = $(CTRILINOS_INCLUDES)
EXPORT_LIBS = $(CTRILINOS_LIBS)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

common_ldadd = $(EXPORT_LIBS)

# This is already added part of CTRILINOS_LIBS and therefore automake does not need to add this again!
#LIBS =


EXTRA_DIST = definition
