// ///////////////////////////////////////////////////////////////////////
// NLPInterfacePackTypes.h

#ifndef NLP_INTERFACE_PACK_TYPES_H
#define NLP_INTERFACE_PACK_TYPES_H

#include "SparseLinAlgPack/include/SparseLinAlgPackTypes.h"

namespace NLPInterfacePack {

// types from SparseLinAlgPack
#include "SparseLinAlgPack/include/SparseLinAlgPackPublicTypes.ud"

// NLP interface classes

class NLP;
class NLPFirstOrderInfo;
class NLPSecondOrderInfo;
class NLPReduced;
class NLPFullToReduced;
class NLPDualCalc;

}	// end namespace NLPInterfacePack 

#endif // NLP_INTERFACE_PACK_TYPES_H