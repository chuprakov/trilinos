// ///////////////////////////////////////////////////////////////////////
// NLPInterfacePackTypes.h
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

#ifndef NLP_INTERFACE_PACK_TYPES_H
#define NLP_INTERFACE_PACK_TYPES_H

#include "AbstractLinAlgPack/include/AbstractLinAlgPackTypes.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h" // Needed for doxygen?

namespace NLPInterfacePack {

#include "AbstractLinAlgPack/include/AbstractLinAlgPackPublicTypes.ud"

// NLP interface classes

class NLP;
class NLPObjGradient;
class NLPFirstOrderDirect;
class NLPFirstOrderInfo;
class NLPSecondOrderInfo;
class NLPReduced;

// NLP testing classes

class NLPFirstDerivativesTester;
class NLPFirstDerivativesTesterSetOptions;
class NLPFirstOrderDirectTester;
class NLPFirstOrderDirectTesterSetOptions;
class NLPTester;
class NLPTesterSetOptions;

// Node implementation classes

class NLPFullToReduced;
class NLPDualCalc;

}	// end namespace NLPInterfacePack 

#endif // NLP_INTERFACE_PACK_TYPES_H
