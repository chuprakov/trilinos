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
  $(srcdir)/MoochoUtilities/src/BLAS_Cpp_Types.hpp\
  $(srcdir)/MoochoUtilities/src/FortranTypes_CppFortranStrings.hpp\
  $(srcdir)/MoochoUtilities/src/FortranTypes_f_open_file.hpp\
  $(srcdir)/MoochoUtilities/src/InputStreamHelperPack_EatInputComment.hpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_OptionsFromStreamExceptions.hpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_OptionsFromStream.hpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_SetOptionsFromStream.hpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_SetOptionsFromStreamNode.hpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_SetOptionsToTargetBase.hpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_StringToBool.hpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_StringToIntMap.hpp\
  $(srcdir)/MoochoUtilities/src/ProfileHackPack_profile_hack.hpp\
  $(srcdir)/MoochoUtilities/src/SerializationPack_Serializable.hpp\
  $(srcdir)/MoochoUtilities/src/SpecialASCII.hpp\
  $(srcdir)/MoochoUtilities/src/StopWatchPack_stopwatch.hpp\
  $(srcdir)/MoochoUtilities/src/StrideIterPack_StrideIter.hpp\
  $(srcdir)/MoochoUtilities/src/TestingHelperPack_update_success.hpp

MOOCHO_UTILITIES_CPP = \
  $(srcdir)/MoochoUtilities/src/FortranTypes_CppFortranStrings.cpp\
  $(srcdir)/MoochoUtilities/src/FortranTypes_f_open_file.cpp\
  $(srcdir)/MoochoUtilities/src/InputStreamHelperPack_EatInputComment.cpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_OptionsFromStream.cpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_SetOptionsFromStreamNode.cpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_StringToBool.cpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_StringToIntMap.cpp\
  $(srcdir)/MoochoUtilities/src/ProfileHackPack_profile_hack.cpp\
  $(srcdir)/MoochoUtilities/src/StopWatchPack_stopwatch.cpp\
  $(srcdir)/MoochoUtilities/src/TestingHelperPack_update_success.cpp

MOOCHO_UTILITIES_F = \
  $(srcdir)/MoochoUtilities/src/ff_open_file.f\
  $(srcdir)/MoochoUtilities/src/FortranCppStrings.f\
  $(srcdir)/MoochoUtilities/src/FSeconds.f

ITERATION_PACK_HPP = \
  $(srcdir)/IterationPack/src/IterationPack_Algorithm.hpp \
  $(srcdir)/IterationPack/src/IterationPack_AlgorithmSetOptions.hpp \
  $(srcdir)/IterationPack/src/IterationPack_AlgorithmState.hpp \
  $(srcdir)/IterationPack/src/IterationPack_AlgorithmStep.hpp \
  $(srcdir)/IterationPack/src/IterationPack_AlgorithmTracker.hpp \
  $(srcdir)/IterationPack/src/IterationPack_AlgorithmTrackerComposite.hpp \
  $(srcdir)/IterationPack/src/IterationPack_CastIQMember.hpp \
  $(srcdir)/IterationPack/src/IterationPack_IterQuantity.hpp \
  $(srcdir)/IterationPack/src/IterationPack_IterQuantityAccess.hpp \
  $(srcdir)/IterationPack/src/IterationPack_IterQuantityAccessContiguous.hpp \
  $(srcdir)/IterationPack/src/IterationPack_IterQuantityAccessContiguousDecl.hpp \
  $(srcdir)/IterationPack/src/IterationPack_IterQuantityAccessContiguousDef.hpp \
  $(srcdir)/IterationPack/src/IterationPack_Types.hpp \
  $(srcdir)/IterationPack/src/IterationPack_cast_iq.hpp \
  $(srcdir)/IterationPack/src/IterationPack_print_algorithm_step.hpp

ITERATION_PACK_CPP = \
  $(srcdir)/IterationPack/src/IterationPack_Algorithm.cpp \
  $(srcdir)/IterationPack/src/IterationPack_AlgorithmSetOptions.cpp \
  $(srcdir)/IterationPack/src/IterationPack_AlgorithmState.cpp \
  $(srcdir)/IterationPack/src/IterationPack_AlgorithmTracker.cpp \
  $(srcdir)/IterationPack/src/IterationPack_AlgorithmTrackerComposite.cpp \
  $(srcdir)/IterationPack/src/IterationPack_CastIQMember.cpp \
  $(srcdir)/IterationPack/src/IterationPack_IterQuantity.cpp \
  $(srcdir)/IterationPack/src/IterationPack_cast_iq.cpp \
  $(srcdir)/IterationPack/src/IterationPack_print_algorithm_step.cpp

RTOP_PACK_BASE_H = \
  $(srcdir)/RTOpPack/src/base/RTOp.h\
  $(srcdir)/RTOpPack/src/base/RTOpToMPI.h\
  $(srcdir)/RTOpPack/src/base/RTOp_SparseSubVector.h\
  $(srcdir)/RTOpPack/src/base/RTOp_apply_op_mpi.h\
  $(srcdir)/RTOpPack/src/base/RTOp_apply_op_serial.h\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_free_free.h\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_index_vtbl.h\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_null_vtbl.h\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_value_index_vtbl.h\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_value_value_vtbl.h\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_value_vtbl.h\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_values_vtbl.h\
  $(srcdir)/RTOpPack/src/base/RTOp_reduct_max_value.h\
  $(srcdir)/RTOpPack/src/base/RTOp_reduct_min_value.h\
  $(srcdir)/RTOpPack/src/base/RTOp_reduct_sum_value.h\
  $(srcdir)/RTOpPack/src/base/RTOp_reduct_sum_values.h\
  $(srcdir)/RTOpPack/src/base/check_nan_inf.h

RTOP_PACK_BASE_C = \
  $(srcdir)/RTOpPack/src/base/RTOp.c\
  $(srcdir)/RTOpPack/src/base/RTOpToMPI.c\
  $(srcdir)/RTOpPack/src/base/RTOp_SparseSubVector.c\
  $(srcdir)/RTOpPack/src/base/RTOp_apply_op_mpi.c\
  $(srcdir)/RTOpPack/src/base/RTOp_apply_op_serial.c\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_free_free.c\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_index_vtbl.c\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_null_vtbl.c\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_value_index_vtbl.c\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_value_value_vtbl.c\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_value_vtbl.c\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_values_vtbl.c\
  $(srcdir)/RTOpPack/src/base/RTOp_reduct_max_value.c\
  $(srcdir)/RTOpPack/src/base/RTOp_reduct_min_value.c\
  $(srcdir)/RTOpPack/src/base/RTOp_reduct_sum_value.c\
  $(srcdir)/RTOpPack/src/base/RTOp_reduct_sum_values.c

RTOP_PACK_BASE_HPP = \
  $(srcdir)/RTOpPack/src/base/RTOpPack_RTOpC.hpp\
  $(srcdir)/RTOpPack/src/base/RTOpPack_RTOpCPostMod.hpp\
  $(srcdir)/RTOpPack/src/base/RTOp_get_reduct_op.hpp

RTOP_PACK_BASE_CPP = \
  $(srcdir)/RTOpPack/src/base/RTOpPack_RTOpC.cpp\
  $(srcdir)/RTOpPack/src/base/check_nan_inf.cpp

RTOp_PACK_STD_OPS_H = \
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_combined_nu_comp_err.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_comp_err_with_mu.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_dot_prod.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_find_nan_inf.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_fraction_to_boundary.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_fraction_to_zero_boundary.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_get_ele.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_get_sub_vector.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_log_bound_barrier.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max_abs_ele.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max_inequ_viol.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max_near_feas_step.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max_rel_step.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max_step.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_norms.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_num_bounded.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_num_nonzeros.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_sum.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_sum_abs.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_add_scalar.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_assign_scalar.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_assign_vectors.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_axpy.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_correct_multipliers.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_ele_wise_divide.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_ele_wise_prod.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_ele_wise_sqrt.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_force_in_bounds.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_inv_of_difference.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_max_abs_vec_scalar.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_max_vec_scalar.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_multiplier_step.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_random_vector.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_scale_vector.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_set_ele.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_set_sub_vector.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_sign.h

RTOp_PACK_STD_OPS_C = \
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_combined_nu_comp_err.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_comp_err_with_mu.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_dot_prod.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_find_nan_inf.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_fraction_to_boundary.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_fraction_to_zero_boundary.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_get_ele.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_get_sub_vector.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_log_bound_barrier.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max_abs_ele.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max_inequ_viol.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max_near_feas_step.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max_rel_step.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max_step.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_norms.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_num_bounded.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_num_nonzeros.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_sum.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_sum_abs.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_add_scalar.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_assign_scalar.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_assign_vectors.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_axpy.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_correct_multipliers.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_ele_wise_divide.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_ele_wise_prod.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_ele_wise_sqrt.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_force_in_bounds.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_inv_of_difference.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_max_abs_vec_scalar.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_max_vec_scalar.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_multiplier_step.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_random_vector.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_scale_vector.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_set_ele.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_set_sub_vector.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_sign.c

RTOp_PACK_STD_OPS_HPP = \
  $(srcdir)/RTOpPack/src/std_ops/RTOpPack_initRTOpStdOpsLibCpp.hpp

RTOp_PACK_STD_OPS_CPP = \
  $(srcdir)/RTOpPack/src/std_ops/RTOpPack_initRTOpStdOpsLibCpp.cpp

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
	$(MOOCHO_UTILITIES_HPP) \
  $(ITERATION_PACK_HPP) \
  $(RTOP_PACK_BASE_H) \
  $(RTOP_PACK_BASE_HPP) \
  $(RTOP_PACK_STD_OPS_H) \
  $(RTOP_PACK_STD_OPS_HPP)

libmoocho_a_SOURCES = \
	$(MOOCHO_UTILITIES_CPP) \
  $(MOOCHO_UTILITIES_F) \
  $(ITERATION_PACK_CPP) \
  $(RTOP_PACK_BASE_C) \
  $(RTOP_PACK_BASE_CPP) \
  $(RTOP_PACK_STD_OPS_C) \
  $(RTOP_PACK_STD_OPS_CPP)

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
