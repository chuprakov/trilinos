



# list of source subdirectories

SUBDIRS = Support Vectors VectorSpaces LinearOperators NonlinearOperators \
Implementations LinearSolvers Preconditioners NonlinearSolvers \
LinearProblem TestProblems Readers XML

# We build lib by going into each subdirectory 
# and invoking make there. 
#
# We make a phony .lib target for each subdirectory name; i.e., if
# the subdirs are A, B, and C, we create targets A.lib, B.lib, and C.lib.
# The rule for making a .lib target is to go into the associated subdirectory
# and running make lib.
# We use this system instead of the simpler foreach because foreach will
# glob all the commands together on one line, with the consequence that 
# make keeps running even after an error is detected.

lib:
	gmake $(LIBS)

clean: 
	$(RM) */*.o */*.d */*~ */#*#

spotless: 
	$(RM) */*.o */*.d */*~ */#*#



# Macro to create a list of subdir names with .sub suffixes
SUBS = $(foreach SUBDIR, $(SUBDIRS), $(SUBDIR).sub)

# Transforms list of .sub names to .lib names
LIBS = $(SUBS:%.sub=%.lib)

# building subdir.lib goes into subdir and runs make lib. Since no real file
# called subdir.lib is ever created, this will always detect an outstanding
# dependency and thus will always re-run. 

%.lib :
	(cd $* ; gmake lib)

include ../Build/Makefile.site