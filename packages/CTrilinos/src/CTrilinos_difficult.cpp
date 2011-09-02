
/*! @HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */


#include "CTrilinos_config.h"
#include "CTrilinos_difficult.hpp"


namespace CTrilinos {


#ifdef HAVE_CTRILINOS_IFPACK
Ifpack::EPrecType convert_to_difficult_enum( CT_EPrecType_E_t en )
{
    switch (en) {
    case CT_EPrecType_E_POINT_RELAXATION:
        return Ifpack::POINT_RELAXATION;
    case CT_EPrecType_E_POINT_RELAXATION_STAND_ALONE:
        return Ifpack::POINT_RELAXATION_STAND_ALONE;
    case CT_EPrecType_E_BLOCK_RELAXATION:
        return Ifpack::BLOCK_RELAXATION;
    case CT_EPrecType_E_BLOCK_RELAXATION_STAND_ALONE:
        return Ifpack::BLOCK_RELAXATION_STAND_ALONE;
    case CT_EPrecType_E_BLOCK_RELAXATION_STAND_ALONE_ILU:
        return Ifpack::BLOCK_RELAXATION_STAND_ALONE_ILU;
#ifdef HAVE_CTRILINOS_AMESOS
    case CT_EPrecType_E_BLOCK_RELAXATION_STAND_ALONE_AMESOS:
        return Ifpack::BLOCK_RELAXATION_STAND_ALONE_AMESOS;
    case CT_EPrecType_E_BLOCK_RELAXATION_AMESOS:
        return Ifpack::BLOCK_RELAXATION_AMESOS;
    case CT_EPrecType_E_AMESOS:
        return Ifpack::AMESOS;
    case CT_EPrecType_E_AMESOS_STAND_ALONE:
        return Ifpack::AMESOS_STAND_ALONE;
#endif /* HAVE_CTRILINOS_AMESOS */
    case CT_EPrecType_E_IC:
        return Ifpack::IC;
    case CT_EPrecType_E_IC_STAND_ALONE:
        return Ifpack::IC_STAND_ALONE;
    case CT_EPrecType_E_ICT:
        return Ifpack::ICT;
    case CT_EPrecType_E_ICT_STAND_ALONE:
        return Ifpack::ICT_STAND_ALONE;
    case CT_EPrecType_E_ILU:
        return Ifpack::ILU;
    case CT_EPrecType_E_ILU_STAND_ALONE:
        return Ifpack::ILU_STAND_ALONE;
    case CT_EPrecType_E_ILUT:
        return Ifpack::ILUT;
    case CT_EPrecType_E_ILUT_STAND_ALONE:
        return Ifpack::ILUT_STAND_ALONE;
#ifdef HAVE_IFPACK_SPARSKIT
    case CT_EPrecType_E_SPARSKIT:
        return Ifpack::SPARSKIT;
#endif /* HAVE_IFPACK_SPARSKIT */
#ifdef HAVE_IFPACK_HIPS
    case CT_EPrecType_E_HIPS:
        return Ifpack::HIPS;
#endif /* HAVE_IFPACK_HIPS */
#ifdef HAVE_HYPRE
    case CT_EPrecType_E_HYPRE:
        return Ifpack::HYPRE;
#endif /* HAVE_HYPRE */
    case CT_EPrecType_E_CHEBYSHEV:
        return Ifpack::CHEBYSHEV;
    default:
        throw CTrilinosMiscException("Likely error in preprocessor macros for Ifpack::EPrecType conversion.\n");
        break;
    }
}
#endif /* HAVE_CTRILINOS_IFPACK */


#ifdef HAVE_CTRILINOS_IFPACK
CT_EPrecType_E_t convert_from_difficult_enum( Ifpack::EPrecType en )
{
    switch (en) {
    case Ifpack::POINT_RELAXATION:
        return CT_EPrecType_E_POINT_RELAXATION;
    case Ifpack::POINT_RELAXATION_STAND_ALONE:
        return CT_EPrecType_E_POINT_RELAXATION_STAND_ALONE;
    case Ifpack::BLOCK_RELAXATION:
        return CT_EPrecType_E_BLOCK_RELAXATION;
    case Ifpack::BLOCK_RELAXATION_STAND_ALONE:
        return CT_EPrecType_E_BLOCK_RELAXATION_STAND_ALONE;
    case Ifpack::BLOCK_RELAXATION_STAND_ALONE_ILU:
        return CT_EPrecType_E_BLOCK_RELAXATION_STAND_ALONE_ILU;
#ifdef HAVE_CTRILINOS_AMESOS
    case Ifpack::BLOCK_RELAXATION_STAND_ALONE_AMESOS:
        return CT_EPrecType_E_BLOCK_RELAXATION_STAND_ALONE_AMESOS;
    case Ifpack::BLOCK_RELAXATION_AMESOS:
        return CT_EPrecType_E_BLOCK_RELAXATION_AMESOS;
    case Ifpack::AMESOS:
        return CT_EPrecType_E_AMESOS;
    case Ifpack::AMESOS_STAND_ALONE:
        return CT_EPrecType_E_AMESOS_STAND_ALONE;
#endif /* HAVE_CTRILINOS_AMESOS */
    case Ifpack::IC:
        return CT_EPrecType_E_IC;
    case Ifpack::IC_STAND_ALONE:
        return CT_EPrecType_E_IC_STAND_ALONE;
    case Ifpack::ICT:
        return CT_EPrecType_E_ICT;
    case Ifpack::ICT_STAND_ALONE:
        return CT_EPrecType_E_ICT_STAND_ALONE;
    case Ifpack::ILU:
        return CT_EPrecType_E_ILU;
    case Ifpack::ILU_STAND_ALONE:
        return CT_EPrecType_E_ILU_STAND_ALONE;
    case Ifpack::ILUT:
        return CT_EPrecType_E_ILUT;
    case Ifpack::ILUT_STAND_ALONE:
        return CT_EPrecType_E_ILUT_STAND_ALONE;
#ifdef HAVE_IFPACK_SPARSKIT
    case Ifpack::SPARSKIT:
        return CT_EPrecType_E_SPARSKIT;
#endif /* HAVE_IFPACK_SPARSKIT */
#ifdef HAVE_IFPACK_HIPS
    case Ifpack::HIPS:
        return CT_EPrecType_E_HIPS;
#endif /* HAVE_IFPACK_HIPS */
#ifdef HAVE_HYPRE
    case Ifpack::HYPRE:
        return CT_EPrecType_E_HYPRE;
#endif /* HAVE_HYPRE */
    case Ifpack::CHEBYSHEV:
        return CT_EPrecType_E_CHEBYSHEV;
    default:
        throw CTrilinosMiscException("Likely error in preprocessor macros for Ifpack::EPrecType conversion.\n");
        break;
    }
}
#endif /* HAVE_CTRILINOS_IFPACK */


} // namespace CTrilinos

