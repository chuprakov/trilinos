// ////////////////////////////////////////////////////////////////////////
// test_nlp_first_order_direct.h

#ifndef TEST_NLP_FIRST_ORDER_DIRECT_H
#define TEST_NLP_FIRST_ORDER_DIRECT_H

#include <iosfwd>

#include "NLPInterfacePack/include/NLPInterfacePackTypes.h"

namespace OptionsFromStreamPack {
	class OptionsFromStream;
}

namespace NLPInterfacePack {

///
/** Test an NLPFirstOrderDirect object.
 *
 * ToDo: Finish documentation!
 */
bool test_nlp_first_order_direct(
	NLPFirstOrderDirect*                          nlp
	,OptionsFromStreamPack::OptionsFromStream*    options
	,std::ostream*                                out
	);

} // end namespace NLPInterfacePack

#endif // TEST_NLP_FIRST_ORDER_DIRECT_H
