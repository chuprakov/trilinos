############################################################################
#
# Makefile for doc subdirectory
#
############################################################################

include ../Build/Makefile.site

manual:
	gmake developersGuide
	gmake usersGuide

developersGuide:
	(cd Doxygen; gmake manual)

installGuide:
	(cd Installation; gmake installGuide)

usersGuide:
	(cd UsersGuide; gmake usersGuide)

clean:
	$(RM) *~ core #*
	$(RM) */*~ */core */#*
	$(RM) -r html latex man


