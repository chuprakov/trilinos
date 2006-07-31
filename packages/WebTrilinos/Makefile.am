#np# Typically the only line that needs to be edited in this file is the last 
#np# line.  Set 'SUBDIRS= (list all subdirectories that are autotool'ed here)'
#np# List all subdirectories even if one or more are configured/built only
#np# conditionally.

# @HEADER
# ************************************************************************
# 
#                     New_Package Example Package
#                 Copyright (2004) Sandia Corporation
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

## #######################################################################
## Options to automake (rarely used - don't worry about it)
## #######################################################################
AUTOMAKE_OPTIONS = foreign

## #######################################################################
## Aclocal command (rarely used - don't worry about it)
## #######################################################################
ACLOCAL	= aclocal -I config

## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
EXTRA_DIST = \
config/ac_cxx_bool.m4                     config/depcomp \
config/ac_cxx_have_std.m4                 config/install-sh \
config/ac_cxx_have_stl.m4                 config/missing \
config/ac_cxx_mutable.m4                  config/mkinstalldirs \
config/ac_cxx_namespaces.m4               config/tac_arg_check_mpi.m4 \
config/ac_cxx_new_for_scoping.m4          config/tac_arg_config_mpi.m4 \
config/ac_cxx_reinterpret_cast.m4         config/generate-makeoptions.pl \
config/replace-install-prefix.pl          config/tac_arg_enable_can_use_package.m4 \
config/tac_arg_enable_compiler_options.m4 config/tac_arg_enable_export-makefiles.m4 \
config/ac_cxx_std_ios_base_fmtflags.m4    config/tac_arg_enable_feature.m4 \
config/ac_cxx_std_sprintf.m4              config/tac_arg_enable_option.m4 \
config/ac_prompt_user.m4                  config/tac_arg_with_ar.m4 \
config/ac_prompt_user_no_define.m4        config/tac_arg_with_blaslib.m4 \
config/acx_blas.m4                        config/tac_arg_with_compiler.m4 \
config/acx_check_dos_filesys.m4           config/tac_arg_with_flags.m4 \
config/acx_lapack.m4                      config/tac_arg_with_incdirs.m4 \
config/acx_pthread.m4                     config/tac_arg_with_lapacklib.m4 \
config/config.guess                       config/tac_arg_with_libdirs.m4 \
config/config.sub                         config/tac_arg_with_package.m4 \
bootstrap \
index.html \
ChangeLog-WebTrilinos     README-WebTrilinos \
common/index.html         common/basic.html \
common/setup.html         common/testing.html \
common/favicon.ico        common/footer.html \
common/util.php           common/config.php \
common/load_images.js     common/header.html \
common/header.php         common/styles.css \
common/WebTrilinos_config.h.in common/WebTrilinos_config.h.in~ \
data/index.html \
data/HB/index.html \
data/HB/fidap005.rua \
data/XML/index.html \
data/XML/sample.xml \
tmp/index.html \
img/beginner.png             img/rd100.png \
img/trilinos_normal.png      img/recirc1.png \
img/trilinos_over.png        img/eth_logo.gif \
img/recirc2.png              img/trilinos_title_normal.png \
img/expert.png               img/sandia_normal.png \
img/trilinos_title_over.png  img/intermediate.png \
img/sandia_over.png          img/webtrilinos_logo.jpg \
img/laplace2d_stencil.png    img/spacer.png \
img/matrix_portal_logo.png   img/star2d_stencil.png \
MatrixPortal/calculator.html MatrixPortal/index.html \
MatrixPortal/step1.html      MatrixPortal/step1C.html \
MatrixPortal/step2.html      MatrixPortal/step2.html \
MatrixPortal/step3.html      MatrixPortal/step4.html \
MatrixPortal/MATLAB_to_XML.m \
MatrixPortal/step_process.py \
MatrixPortal/help/footer.html \
MatrixPortal/help/header.html \
MatrixPortal/help/step1_collections.html \
MatrixPortal/help/step1_generator.html \
MatrixPortal/help/step1_upload.html \
MatrixPortal/help/step2_multilevel.html \
MatrixPortal/help/step2_incomplete_factorization.html \
MatrixPortal/help/step2_performance.html \
MatrixPortal/help/step2_general.html \
MatrixPortal/help/step2_relaxation.html \
MatrixPortal/help/step2_chebyshev.html \
MatrixPortal/help/step4_plot.html \
MatrixPortal/help/step_workflow.html \
MatrixPortal/help/template.html

## #######################################################################
## Tools in the auxillary directory 
## #######################################################################
AUX_DIST = config/install-sh config/missing config/mkinstalldirs 
#
#  Again, I hope that AUX_DIR(../../config) eliminates these
#  config/install-sh config/missing config/mkinstalldirs 

## #######################################################################
## Files to be deleted by 'make maintainer-clean'
## #######################################################################
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 autom4te.cache/* \
	configure config.status config.log \
	$(AUX_DIST) 

all-local:
	@echo "Trilinos package WebTrilinos built successfully"
	@echo 
	@echo "        +----------------------------------------------------+"
	@echo "        | Trilinos package WebTrilinos built successfully    |"
	@echo "        |                                                    |"
	@echo "        | The WebTrilinos Users' Guide is located in         |"
	@echo "        | Trilinos/packages/WebTrilinos/doc/UG.pdf. You can  |"
	@echo "        | also consult the WebTrilinos web site for more     |"
	@echo "        | details. ** Note that WebTrilinos must be manually |"
	@echo "        | installed to be operational. **                    |"
	@echo "        |                                                    |"
	@echo "        | You can type \"make run-tests\" to verify that the   |"
	@echo "        | c++ examples work, and \"make run-pytests\" to verify|"
	@echo "        | that the python scripts are ok.                    |"
	@echo "        |                                                    |"
	@echo "        | Please send questions/comments to:                 |"
	@echo "        | - WebTrilinos-users@software.sandia.gov            |"
	@echo "        | - WebTrilinos-developers@software.sandia.gov       |"
	@echo "        |                                                    |"
	@echo "        | Have fun!                                          |"
	@echo "        +----------------------------------------------------+"
	@echo 

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################

SUBDIRS = c++ python

## #######################################################################
## Other make targets
## #######################################################################

TRILINOS_HOME_DIR=@abs_top_srcdir@/../..
TRILINOS_BUILD_DIR=@abs_top_builddir@/../..

if HAVE_MPI
WEBTRILINOS_CHECK_COMM=mpi
else
WEBTRILINOS_CHECK_COMM=serial
endif

run-tests:
	@echo
	@echo "================================"
	@echo "Running the WebTrilinos tests..."
	@echo "================================"
	@echo
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(WEBTRILINOS_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=WebTrilinosTests \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=WebTrilinos

run-pytests:
	@echo
	@echo "======================================="
	@echo "Running the Python/WebTrilinos tests..."
	@echo "======================================="
	@echo
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(WEBTRILINOS_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=WebTrilinosPyTests \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=WebTrilinos
