// /////////////////////////////////
// SparseVectorClass.cpp

#include "SparseLinAlgPack/include/SparseVectorClassDecl.h"

void SparseLinAlgPack::SparseVectorUtilityPack::assert_is_sorted(bool is_sorted)
{
	if(!is_sorted)
		throw NotSortedException("SparseVector***<> : The sparse vector is not sorted.");

}
