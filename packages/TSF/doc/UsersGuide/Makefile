# Makefile for documentation
#

include ../../Build/Makefile.site

usersGuide:
	gmake ps
	gmake html

ps:
	latex UsersGuide
	latex UsersGuide
	dvips -o UsersGuide.ps UsersGuide

html: 
	latex2html -init_file latex2html.init -split 5 -no_math UsersGuide