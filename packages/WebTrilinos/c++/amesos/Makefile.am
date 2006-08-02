include $(top_builddir)/Makefile.export

EXAMPLE = ex_simple
EXEEXT = .exe

noinst_PROGRAMS = $(EXAMPLE)

ex_simple_SOURCES = $(srcdir)/ex_simple.cpp
ex_simple_CXXFLAGS = $(WEBTRILINOS_INCLUDES) 
ex_simple_LDADD = $(WEBTRILINOS_LIBS)
