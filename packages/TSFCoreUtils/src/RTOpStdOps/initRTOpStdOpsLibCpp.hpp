// ////////////////////////////////////////////////
// initRTOpStdOpsLibCpp.hpp

#ifndef INIT_RTOP_STD_OPS_LIB_CPP_HPP
#define INIT_RTOP_STD_OPS_LIB_CPP_HPP

#include "RTOpServer.hpp"

namespace RTOpPack {

///
/** Initialize an <tt>RTOpPack::RTOpServer</tt> object for all of the
 * operators defined in this library.
 */
void initRTOpStdOpsLibCpp( RTOpPack::RTOpServer<RTOp_value_type> *op_server );

} // namespace RTOpPack

#endif // INIT_RTOP_STD_OPS_LIB_CPP_HPP
