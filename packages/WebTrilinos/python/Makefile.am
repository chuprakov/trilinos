EXTRA_DIST = index.html help.html \
	     setpath.py \
             ex_Amesos.py  ex_Comm.py \
	     ex_Map.py     ex_ML.py \
	     ex_AztecOO.py ex_IFPACK.py \
	     ex_Matrix.py  ex_MultiVector.py  

EXAMPLE_SUITE = ex_Amesos  ex_Comm\
                ex_Map     ex_ML \
                ex_AztecOO ex_IFPACK \
                ex_Matrix  ex_MultiVector  

SETPATH = setpath.py

all: $(EXAMPLE_SUITE)
	        if [ $(top_srcdir) != @abs_top_builddir@ ]; then cp $(srcdir)/$(SETPATH) . ; fi

ex%: $(srcdir)/ex%.py
	        if [ $(top_srcdir) != @abs_top_builddir@ ]; then cp $(srcdir)/$@.py . ; fi

CLEANFILES = $(EXAMPLE_SUITE:=.py) $(SETPATH) *.pyc
