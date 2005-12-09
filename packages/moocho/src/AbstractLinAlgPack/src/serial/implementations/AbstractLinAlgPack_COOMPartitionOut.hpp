// //////////////////////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_COOMPartitionOut.hpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.
//
// Output stream operators for Partition<> and TransposedPartition<>

#ifndef COOM_PARTITION_OUT_H
#define COOM_PARTITION_OUT_H

#include "AbstractLinAlgPack_COOMatrixTmplOutFunc.hpp"

namespace AbstractLinAlgPack {

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

}	// end namespace AbstractLinAlgPack

#endif // VECTOROUT_H
