// /////////////////////////////////////////////////////
// NLPFirstOrderDirect.cpp

#include <assert.h>

#include "NLPInterfacePack/include/NLPFirstOrderDirect.h"
#include "AbstractLinAlgPack/include/MatrixWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "Range1D.h"
#include "ThrowException.h"

namespace NLPInterfacePack {

// NLPFirstOrderDirect

size_type NLPFirstOrderDirect::r() const
{
	return this->con_decomp().size();
}

Range1D NLPFirstOrderDirect::var_dep() const
{
	return Range1D(1,m());
}
Range1D NLPFirstOrderDirect::var_indep() const
{
	return Range1D(m()+1,n());
}
Range1D NLPFirstOrderDirect::con_decomp() const
{
	return Range1D(1,m());
}

Range1D NLPFirstOrderDirect::con_undecomp() const
{
	return Range1D::Invalid;
}

const NLPFirstOrderDirect::mat_space_ptr_t&
NLPFirstOrderDirect::space_GcU() const
{
	return NULL;
}

const NLPFirstOrderDirect::mat_space_ptr_t&
NLPFirstOrderDirect::space_Gh() const
{
	return NULL;
}

const NLPFirstOrderDirect::mat_space_ptr_t&
NLPFirstOrderDirect::space_V() const
{
	return NULL;
}

const NLPFirstOrderDirect::mat_space_ptr_t&
NLPFirstOrderDirect::space_P() const
{
	return NULL;
}

const NLPFirstOrderDirect::mat_mut_space_ptr_t&
NLPFirstOrderDirect::space_GcUD() const
{
	return NULL;
}

const NLPFirstOrderDirect::mat_mut_space_ptr_t&
NLPFirstOrderDirect::space_GhD() const
{
	return NULL;
}

size_type NLPFirstOrderDirect::mI() const
{
	return 0;
}

NLP::vec_space_ptr_t NLPFirstOrderDirect::space_h() const
{
	return NULL;
}

const VectorWithOp& NLPFirstOrderDirect::hl() const
{
	THROW_EXCEPTION( true, NoBounds, "NLPFirstOrderDirect::hl(), Error, default is for mI() == 0" );
	return xl(); // will never execute.
}

const VectorWithOp& NLPFirstOrderDirect::hu() const
{
	THROW_EXCEPTION( true, NoBounds, "NLPFirstOrderDirect::hl(), Error, default is for mI() == 0" );
	return xu(); // will never execute.
}

void NLPFirstOrderDirect::imp_calc_h(
	const VectorWithOp& x, bool newx, const ZeroOrderInfo& zero_order_info) const
{
	assert(0); // Should never be called!
}

} // end namespace NLPIntefacePack
