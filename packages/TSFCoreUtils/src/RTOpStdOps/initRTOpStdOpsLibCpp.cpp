// ////////////////////////////////////////////////
// initRTOpStdOpsLibCpp.cpp

#include <assert.h>

#include "initRTOpStdOpsLibCpp.hpp"
#include "RTOpCPostMod.hpp"
#include "AbstractFactoryStd.hpp"
#include "RTOp_ROp_combined_nu_comp_err.h"
#include "RTOp_ROp_max_rel_step.h"
#include "RTOp_TOp_ele_wise_sqrt.h"
#include "RTOp_ROp_comp_err_with_mu.h"
#include "RTOp_ROp_max_step.h"
#include "RTOp_TOp_force_in_bounds.h"
#include "RTOp_ROp_dot_prod.h"
#include "RTOp_ROp_norms.h"
#include "RTOp_TOp_inv_of_difference.h"
#include "RTOp_ROp_find_nan_inf.h"
#include "RTOp_ROp_num_bounded.h"
#include "RTOp_TOp_max_abs_vec_scalar.h"
#include "RTOp_ROp_fraction_to_boundary.h"
#include "RTOp_ROp_num_nonzeros.h"
#include "RTOp_TOp_max_vec_scalar.h"
#include "RTOp_ROp_fraction_to_zero_boundary.h"
#include "RTOp_ROp_sum.h"
#include "RTOp_ROp_sum_abs.h"
#include "RTOp_TOp_multiplier_step.h"
#include "RTOp_ROp_get_ele.h"
#include "RTOp_TOp_add_scalar.h"
#include "RTOp_TOp_random_vector.h"
#include "RTOp_ROp_get_sub_vector.h"
#include "RTOp_TOp_assign_scalar.h"
#include "RTOp_TOp_scale_vector.h"
#include "RTOp_ROp_log_bound_barrier.h"
#include "RTOp_TOp_assign_vectors.h"
#include "RTOp_TOp_set_ele.h"
#include "RTOp_ROp_max_abs_ele.h"
#include "RTOp_TOp_axpy.h"
#include "RTOp_TOp_set_sub_vector.h"
#include "RTOp_ROp_max.h"
#include "RTOp_TOp_correct_multipliers.h"
#include "RTOp_TOp_sign.h"
#include "RTOp_ROp_max_inequ_viol.h"
#include "RTOp_TOp_ele_wise_divide.h"
#include "RTOp_ROp_max_near_feas_step.h"
#include "RTOp_TOp_ele_wise_prod.h"

namespace {

void add_op_factory(RTOpPack::RTOpServer<RTOp_value_type> *op_server, const RTOp_RTOp_vtbl_t &op_vtbl )
{
	namespace mmp = MemMngPack;
	namespace rtop = RTOpPack;
	typedef mmp::AbstractFactoryStd<rtop::RTOp,rtop::RTOpC,rtop::RTOpCPostMod>  OpFactory_t;
	op_server->add_op_factory(
		mmp::rcp(
			new OpFactory_t(
				rtop::RTOpCPostMod(&op_vtbl)
				,mmp::AllocatorNew<rtop::RTOpC>() // RAB: 2003/10/03: MipsPro needs this *$&%^#&!!!
				)
			)
		);
}

} // namespace

void RTOpPack::initRTOpStdOpsLibCpp( RTOpPack::RTOpServer<RTOp_value_type> *op_server )
{
	namespace mmp = MemMngPack;
	typedef mmp::AbstractFactoryStd<RTOp,RTOpC,RTOpCPostMod>  OpFactory_t;
	
	add_op_factory(op_server,RTOp_ROp_combined_nu_comp_err_vtbl);
	add_op_factory(op_server,RTOp_ROp_combined_nu_comp_err_one_only_vtbl);
	add_op_factory(op_server,RTOp_ROp_comp_err_with_mu_vtbl);
	add_op_factory(op_server,RTOp_ROp_dot_prod_vtbl);
	add_op_factory(op_server,RTOp_ROp_find_nan_inf_vtbl);
	add_op_factory(op_server,RTOp_ROp_fraction_to_boundary_vtbl);
	add_op_factory(op_server,RTOp_ROp_fraction_to_zero_boundary_vtbl);
	add_op_factory(op_server,RTOp_ROp_get_ele_vtbl);
	add_op_factory(op_server,RTOp_ROp_get_sub_vector_vtbl);
	add_op_factory(op_server,RTOp_ROp_log_bound_barrier_vtbl);
	add_op_factory(op_server,RTOp_ROp_max_abs_ele_vtbl);
	add_op_factory(op_server,RTOp_ROp_max_vtbl);
	add_op_factory(op_server,RTOp_ROp_max_inequ_viol_vtbl);
	add_op_factory(op_server,RTOp_ROp_max_near_feas_step_vtbl);
	add_op_factory(op_server,RTOp_ROp_max_rel_step_vtbl);
	add_op_factory(op_server,RTOp_ROp_max_step_vtbl);
	add_op_factory(op_server,RTOp_ROp_norm_1_vtbl);
	add_op_factory(op_server,RTOp_ROp_norm_2_vtbl);
	add_op_factory(op_server,RTOp_ROp_norm_inf_vtbl);
	add_op_factory(op_server,RTOp_ROp_num_bounded_vtbl);
	add_op_factory(op_server,RTOp_ROp_num_nonzeros_vtbl);
	add_op_factory(op_server,RTOp_ROp_sum_vtbl);
	add_op_factory(op_server,RTOp_ROp_sum_abs_vtbl);
	add_op_factory(op_server,RTOp_TOp_add_scalar_vtbl);
	add_op_factory(op_server,RTOp_TOp_assign_scalar_vtbl);
	add_op_factory(op_server,RTOp_TOp_assign_vectors_vtbl);
	add_op_factory(op_server,RTOp_TOp_axpy_vtbl);
	add_op_factory(op_server,RTOp_TOp_Correct_Multipliers_vtbl);
	add_op_factory(op_server,RTOp_TOp_ele_wise_divide_vtbl);
	add_op_factory(op_server,RTOp_TOp_ele_wise_prod_vtbl);
	add_op_factory(op_server,RTOp_TOp_ele_wise_sqrt_vtbl);
	add_op_factory(op_server,RTOp_TOp_force_in_bounds_vtbl);
	add_op_factory(op_server,RTOp_TOp_force_in_bounds_buffer_vtbl);
	add_op_factory(op_server,RTOp_TOp_inv_of_difference_vtbl);
	add_op_factory(op_server,RTOp_TOp_max_abs_vec_scalar_vtbl);
	add_op_factory(op_server,RTOp_TOp_max_vec_scalar_vtbl);
	add_op_factory(op_server,RTOp_TOp_multiplier_step_vtbl);
	add_op_factory(op_server,RTOp_TOp_random_vector_vtbl);
	add_op_factory(op_server,RTOp_TOp_scale_vector_vtbl);
	add_op_factory(op_server,RTOp_TOp_set_ele_vtbl);
	add_op_factory(op_server,RTOp_TOp_set_sub_vector_vtbl);
	add_op_factory(op_server,RTOp_TOp_sign_vtbl);
}
