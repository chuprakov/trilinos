// ////////////////////////////////////////////////////////////////////
// SpVectorView.h

#ifndef SP_VECTOR_VIEW_H
#define SP_VECTOR_VIEW_H

#include "SpVectorClass.h"

namespace AbstractLinAlgPack {

///
/** Create an RTOp_SubVector view object from a SpVectorSlice object.
 */
RTOp_SubVector sub_vec_view( const SpVectorSlice& sv );

} // end namespace AbstractLinAlgPack

#endif // SP_VECTOR_VIEW_H
