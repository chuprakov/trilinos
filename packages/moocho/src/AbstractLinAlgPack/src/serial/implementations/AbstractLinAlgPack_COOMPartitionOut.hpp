// //////////////////////////////////////////////////////////////////////////////////////
// COOMPartitionOut.h
//
// Output stream operators for Partition<> and TransposedPartition<>

#ifndef COOM_PARTITION_OUT_H
#define COOM_PARTITION_OUT_H

#include "COOMatrixTmplOutFunc.h"

namespace SparseLinAlgPack {

///
/** Partition<> output stream operator.
  *
  * This operator function calls the function output_COOM(os,part,0).
  */
template <class T_Indice, class T_Value>
inline std::ostream& operator<<(std::ostream& os
	, const COOMatrixPartitionedViewUtilityPack::Partition<T_Indice,T_Value>& part) {
	return output_COOM(os,part,0);
}

///
/** TransposedPartition<> output stream operator.
  *
  * This operator function calls the function output_COOM(os,trans_part,0).
  */
template <class T_Indice, class T_Value>
inline std::ostream& operator<<(std::ostream& os
	, const COOMatrixPartitionedViewUtilityPack::TransposedPartition<T_Indice,T_Value>& trans_part)
{
	return output_COOM(os,trans_part,0);	
}

}	// end namespace SparseLinAlgPack

#endif // VECTOROUT_H