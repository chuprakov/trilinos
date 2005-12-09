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
  $(srcdir)/RTOpPack/src/base/RTOp_get_reduct_op.hpp\
  $(srcdir)/RTOpPack/src/base/RTOpPack_print_sub_vector.hpp

RTOP_PACK_BASE_CPP = \
  $(srcdir)/RTOpPack/src/base/RTOpPack_RTOpC.cpp\
  $(srcdir)/RTOpPack/src/base/check_nan_inf.cpp\
  $(srcdir)/RTOpPack/src/base/RTOpPack_print_sub_vector.cpp

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

DENSE_LIN_ALG_PACK_HPP = \
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgLAPack.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgLAPack_Types.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_AssertOp.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_BLAS_Cpp.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixAsTriSym.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixAssign.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixClass.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixIn.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixInFunc.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixOp.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixOpBLAS.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixOut.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixOutFunc.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorAssign.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorClass.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorClassExt.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorClassTmpl.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorIn.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgLAPack_Types.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_AssertOp.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_BLAS_Cpp.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixAsTriSym.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixAssign.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixClass.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixIn.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixInFunc.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixOp.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixOpBLAS.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixOut.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixOutFunc.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorAssign.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorClass.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorClassExt.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorClassTmpl.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorIn.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorInFunc.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorOp.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorOut.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorOutFunc.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_ElementaryMatVec.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_IO.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_IOBasic.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_IOFormat.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_IO_NameLookups.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_IVector.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_InFormat.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_InFormatDecl.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_InFormatDef.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_InvCholUpdate.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_LAPACK_C_Decl.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_LAPACK_Cpp.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_LinAlgOpPack.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_MatVecCompare.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_MatlabPack.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_Options.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_OutFormat.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_PermIn.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_PermOut.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_PermVecMat.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_Types.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_assert_print_nan_inf.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_delete_row_col.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_extended_value_type.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_random_vector.hpp

DENSE_LIN_ALG_PACK_CPP = \
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_BLAS_Cpp.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixClass.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixInFunc.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixOpBLAS.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixOutFunc.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_InvCholUpdate.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_AssertOp.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_IOFormat.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_MatVecCompare.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_MatlabPack.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_PermIn.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_PermOut.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_PermVecMat.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorClassTmpl.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorInFunc.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorOpBLAS.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorOutFunc.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_random_vector.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_delete_row_col.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_assert_print_nan_inf.cpp \
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_LAPACK_C_Decl.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgLAPack.cpp

ABSTRACT_LIN_ALG_PACK_HPP = \
  $(srcdir)/AbstractLinAlgPack/src/AbstractLinAlgPack_Types.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_AssertOp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_BasisSystem.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_BasisSystemFactory.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_BasisSystemPerm.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_EtaVector.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_GenPermMatrixSlice.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_GenPermMatrixSliceIterator.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_GenPermMatrixSliceOut.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_InnerProduct.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_InnerProductDot.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_LinAlgOpPack.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixBase.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixNonsing.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixOp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixOpNonsing.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixOpNonsingAggr.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixOpOut.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixOpSubView.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixPermAggr.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixSymDiag.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixSymInitDiag.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixSymNonsing.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixSymOp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixSymOpNonsing.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixSymSecant.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MultiVector.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MultiVectorMutable.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_Permutation.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_PermutationOut.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SpVecIndexLookupClass.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SpVecIndexLookupClassDecl.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SpVecIndexLookupClassDef.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SpVectorClass.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SpVectorOut.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SpVectorView.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SparseElement.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SparseVectorClass.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SparseVectorClassDecl.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SparseVectorClassDef.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_Vector.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorApplyOpSerialBase.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorMutable.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorMutableSubView.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorOut.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorSpace.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorSpaceFactory.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorSpaceSubSpace.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorStdOps.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorSubView.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_apply_op_helper.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_compare_element_indexes.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/thyra/AbstractLinAlgPack_InnerProductThyra.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/thyra/AbstractLinAlgPack_MatrixOpNonsingThyra.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/thyra/AbstractLinAlgPack_MatrixOpThyra.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/thyra/AbstractLinAlgPack_MultiVectorMutableThyra.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/thyra/AbstractLinAlgPack_ThyraAccessors.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/thyra/AbstractLinAlgPack_VectorMutableThyra.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/thyra/AbstractLinAlgPack_VectorSpaceFactoryThyra.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/thyra/AbstractLinAlgPack_VectorSpaceThyra.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_BasisSystemComposite.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_BasisSystemTester.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_BasisSystemTesterSetOptions.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_MatrixComposite.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_MatrixOpNonsingTester.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_MatrixSymDiagStd.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_MatrixSymIdent.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_MatrixZero.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_TestMatrixSymSecant.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_VectorAuxiliaryOps.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_VectorMutableBlocked.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_VectorSpaceBlocked.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_VectorSpaceTester.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_VectorSpaceTesterSetOptions.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_assert_print_nan_inf.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_BFGS_helpers.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMPartitionOp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMPartitionOut.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixClass.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixIn.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixOut.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixOutFunc.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixPartitionViewSubclass.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixPartitionedViewClass.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixPartitionedViewClassDecl.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixPartitionedViewClassDef.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixTmplConvertToSparseCompressedColumn.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixTmplConvertToSparseCompressedColumnDecl.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixTmplConvertToSparseCompressedColumnDef.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixTmplOp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixTmplOpDecl.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixTmplOpDef.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixTmplOutFunc.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixTmplOutFuncDecl.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixTmplOutFuncDef.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixWithPartitionedView.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixWithPartitionedViewSubclass.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_MatrixBasisNonbasisStd.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_MatrixCOORTmplItfc.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_MatrixSparseCOORSerial.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_MatrixSymDiagSparseStd.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_MatrixSymPosDefCholFactor.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_MatrixVectorTemplateOp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_MatrixVectorTemplateOpDef.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_MatrixWithOpConcreteEncap.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_PermutationSerial.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SortByDescendingAbsValue.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SpVectorOp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparseCOOElementTransView.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparseCOOPtrElement.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparseCOOReadMatrix.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparseElementOut.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparsePtrElement.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparseVectorOp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparseVectorOpDecl.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparseVectorOpDef.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparseVectorOutput.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparseVectorOutputDef.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparseVectorSliceOp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_TransSparseCOOElementViewIter.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_VectorDenseEncap.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_dense_Vp_StPtMtV.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_rank_2_chol_update.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_sparse_bounds.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_sparse_bounds_diff.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_ConvertToCSC.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_GenPermMatrixSliceOp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_LinAlgOpPackHack.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixBasisNonbasis.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixConvertToSparse.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixConvertToSparseEncap.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixExtractInvCholFactor.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixExtractSparseElements.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixLoadSparseElements.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixNonsingSerial.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixOpGetGMS.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixOpGetGMSMutable.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixOpGetGMSTri.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixOpNonsingSerial.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixOpSerial.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymAddDelUpdateable.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymDenseInitialize.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymDiagSparse.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymNonsingSerial.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymOpGetGMSSym.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymOpGetGMSSymMutable.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymOpNonsingSerial.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymOpSerial.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MultiVectorMutableCols.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MultiVectorMutableDense.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_VectorMutableDense.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_VectorSpaceFactorySerial.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_VectorSpaceSerial.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_BasisSystemFactoryStd.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_BasisSystemPermDirectSparse.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolver.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolverDense.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolverImp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolverMA28.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolverMA28SetOptions.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolverSuperLU.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_MA28CommonBlockEncap.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_MA28Solver.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_MA28_CppDecl.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_MatrixScaling_Strategy.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_SuperLUSolver.hpp

ABSTRACT_LIN_ALG_PACK_CPP = \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_AssertOp.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_AbstractLinAlgOpPack.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_apply_op_helper.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_BasisSystem.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_GenPermMatrixSlice.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_GenPermMatrixSliceIterator.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_GenPermMatrixSliceOut.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_InnerProductDot.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixBase.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixOp.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixNonsing.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixOpNonsing.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixOpNonsingAggr.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixOpSubView.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixSymOp.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixSymNonsing.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixSymOpNonsing.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixPermAggr.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MultiVector.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MultiVectorMutable.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_Vector.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorApplyOpSerialBase.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorSpace.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorSpaceSubSpace.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorSpaceFactory.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorSubView.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorMutable.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorMutableSubView.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorStdOps.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SparseVectorClass.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SpVectorOut.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SpVectorView.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_assert_print_nan_inf.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_BasisSystemComposite.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_BasisSystemTester.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_BasisSystemTesterSetOptions.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_MatrixOpNonsingTester.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_MatrixSymIdent.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_MatrixSymDiagStd.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_MatrixComposite.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_MatrixZero.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_VectorSpaceTester.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_VectorSpaceTesterSetOptions.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_VectorSpaceBlocked.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_VectorMutableBlocked.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_VectorAuxiliaryOps.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_TestMatrixSymSecant.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/thyra/AbstractLinAlgPack_InnerProductThyra.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/thyra/AbstractLinAlgPack_VectorSpaceFactoryThyra.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/thyra/AbstractLinAlgPack_VectorSpaceThyra.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/thyra/AbstractLinAlgPack_VectorMutableThyra.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/thyra/AbstractLinAlgPack_MatrixOpThyra.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/thyra/AbstractLinAlgPack_MatrixOpNonsingThyra.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/thyra/AbstractLinAlgPack_MultiVectorMutableThyra.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/thyra/AbstractLinAlgPack_ThyraAccessors.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_LinAlgOpPackHack.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_GenPermMatrixSliceOp.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_VectorSpaceSerial.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_VectorSpaceFactorySerial.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_VectorMutableDense.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MultiVectorMutableDense.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MultiVectorMutableCols.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixConvertToSparseEncap.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixOpSerial.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixNonsingSerial.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymOpSerial.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymNonsingSerial.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymDiagSparse.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixExtractSparseElements.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_BFGS_helpers.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_MatrixSymDiagSparseStd.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_MatrixSymPosDefCholFactor.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_MatrixSparseCOORSerial.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_PermutationSerial.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SpVectorOp.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_rank_2_chol_update.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_BasisSystemPermDirectSparse.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_BasisSystemFactoryStd.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolverImp.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolverDense.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolverMA28.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolverMA28SetOptions.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolverSuperLU.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_MA28_CppDecl.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_MA28CommonBlockEncap.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_MA28Solver.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_SuperLUSolver.cpp

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
	$(RTOP_PACK_STD_OPS_HPP) \
	$(DENSE_LIN_ALG_PACK_HPP) \
	$(ABSTRACT_LIN_ALG_PACK_HPP)

libmoocho_a_SOURCES = \
	$(MOOCHO_UTILITIES_CPP) \
  $(MOOCHO_UTILITIES_F) \
  $(ITERATION_PACK_CPP) \
  $(RTOP_PACK_BASE_C) \
  $(RTOP_PACK_BASE_CPP) \
  $(RTOP_PACK_STD_OPS_C) \
  $(RTOP_PACK_STD_OPS_CPP) \
  $(DENSE_LIN_ALG_PACK_CPP) \
	$(ABSTRACT_LIN_ALG_PACK_CPP)

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
