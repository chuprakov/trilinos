// //////////////////////////////////////////////////////////////////////////////////
// LinAlgOpPack.hpp
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

#ifndef SLAP_LIN_ALG_OP_PACK_H
#define SLAP_LIN_ALG_OP_PACK_H

#include "DenseLinAlgPack/src/LinAlgOpPack.hpp"
#include "MatrixSymOpNonsingSerial.hpp"

namespace LinAlgOpPack {

using AbstractLinAlgPack::MatrixOpSerial;
using AbstractLinAlgPack::MatrixNonsingSerial;
using AbstractLinAlgPack::MatrixOpNonsingSerial;
using AbstractLinAlgPack::MatrixSymOpSerial;
using AbstractLinAlgPack::MatrixSymOpNonsingSerial;

using AbstractLinAlgPack::Vp_StMtV;
using AbstractLinAlgPack::V_InvMtV;
using AbstractLinAlgPack::Mp_StM;
using AbstractLinAlgPack::Mp_StMtM;
using AbstractLinAlgPack::M_StInvMtM;

} // end namespace LinAlgOpPack

#endif // SLAP_LIN_ALG_OP_PACK_H
