// ////////////////////////////////////////////////////////////////////////////
// VectorSpaceSerial.cpp

#include <assert.h>

#include "SparseLinAlgPack/include/VectorSpaceSerial.h"
#include "SparseLinAlgPack/include/VectorWithOpMutableDense.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "LinAlgPack/include/VectorClass.h"

namespace SparseLinAlgPack {

VectorSpaceSerial::VectorSpaceSerial( size_type dim )
{
	initialize(dim);
}

void VectorSpaceSerial::initialize( size_type dim )
{
	dim_ = dim;
}

bool VectorSpaceSerial::is_compatible(const VectorSpace& a_vec_space ) const
{
	return this->dim() == a_vec_space.dim();
}

// Overridden from VectorSpace

index_type VectorSpaceSerial::dim() const
{
	return dim_;
}

VectorSpace::space_ptr_t
VectorSpaceSerial::clone() const
{
	namespace rcp = ReferenceCountingPack;
	return rcp::rcp( new VectorSpaceSerial( dim_	) );
}

VectorSpace::vec_mut_ptr_t
VectorSpaceSerial::create_member() const
{
	namespace rcp = ReferenceCountingPack;
	return rcp::rcp(new VectorWithOpMutableDense(dim_));
}

} // end namespace SparseLinAlgPack
