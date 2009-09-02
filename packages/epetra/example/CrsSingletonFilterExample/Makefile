TRILINOS_INSTALLATION=..
PUC=IFPACK
PLC=ifpack
ROOT=cxx_main
 
include $(TRILINOS_INSTALLATION)/include/Makefile.export.$(PLC).macros
include $(TRILINOS_INSTALLATION)/include/Makefile.export.$(PLC)
include $(TRILINOS_INSTALLATION)/include/Makefile.export.epetraext.macros
include $(TRILINOS_INSTALLATION)/include/Makefile.export.epetraext
 
INCDIR=$($(PUC)_INCLUDES) $(EPETRAEXT_INCLUDES)
LIBDIR=$($(PUC)_LIBS) $(EPETRAEXT_LIBS)
 
CXXFLAGS=$(INCDIR) $($(PLC)_CXXFLAGS) $($(PLC)_DEFS) $($(PUC)_CXXFLAGS) $($(PUC)_DEFS)
CXX=$($(PLC)_CXX) $($(PUC)_CXX)
 
TARGET=$(ROOT).exe
SRC=$(ROOT).cpp
 
OBJ=$(SRC:.cpp=.o)
 
$(TARGET): $(OBJ)
	$(CXX) -o $(TARGET) $(OBJ) $(LIBDIR)
 
clean:
	rm -rf $(OBJ) $(TARGET) *~
 
again: clean $(TARGET)
