
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
#include "CTrilinos_utils.hpp"


namespace CTrilinos {


void pass_bool_out( const bool * pval, boolean * pvalout )
{
    pvalout = new boolean[1];
    *pvalout = (*pval ? TRUE : FALSE);
}


void pass_bool_in( const boolean * pval, bool * pvalout )
{
    pvalout = new bool[1];
    *pvalout = (*pval != FALSE ? true : false);
}


void pass_string_out( const std::string * const s, char *c[] )
{
    *c = new char[s->size()+1];
    strcpy(*c, s->c_str());
}


void pass_string_in( const char * const c[], std::string *s )
{
    s = new std::string(*c);
}


/* stringify the enum name */
std::string
enum2str( CTrilinos_Table_ID_t ty )
{
    std::string s;
    switch (ty) {
    case CT_Epetra_Distributor_ID:
        s.assign("CT_Epetra_Distributor_ID");
        break;
    case CT_Epetra_SerialComm_ID:
        s.assign("CT_Epetra_SerialComm_ID");
        break;
    case CT_Epetra_BLAS_ID:
        s.assign("CT_Epetra_BLAS_ID");
        break;
    case CT_Epetra_Comm_ID:
        s.assign("CT_Epetra_Comm_ID");
        break;
    case CT_Epetra_Operator_ID:
        s.assign("CT_Epetra_Operator_ID");
        break;
    case CT_Epetra_MultiVector_ID:
        s.assign("CT_Epetra_MultiVector_ID");
        break;
    case CT_Epetra_OffsetIndex_ID:
        s.assign("CT_Epetra_OffsetIndex_ID");
        break;
    case CT_Epetra_Object_ID:
        s.assign("CT_Epetra_Object_ID");
        break;
    case CT_Epetra_RowMatrix_ID:
        s.assign("CT_Epetra_RowMatrix_ID");
        break;
    case CT_Epetra_CompObject_ID:
        s.assign("CT_Epetra_CompObject_ID");
        break;
    case CT_Epetra_Directory_ID:
        s.assign("CT_Epetra_Directory_ID");
        break;
    case CT_Epetra_Flops_ID:
        s.assign("CT_Epetra_Flops_ID");
        break;
    case CT_Epetra_SrcDistObject_ID:
        s.assign("CT_Epetra_SrcDistObject_ID");
        break;
    case CT_Epetra_MpiComm_ID:
        s.assign("CT_Epetra_MpiComm_ID");
        break;
    case CT_Epetra_CrsMatrix_ID:
        s.assign("CT_Epetra_CrsMatrix_ID");
        break;
    case CT_Epetra_CrsGraph_ID:
        s.assign("CT_Epetra_CrsGraph_ID");
        break;
    case CT_Epetra_DistObject_ID:
        s.assign("CT_Epetra_DistObject_ID");
        break;
    case CT_Epetra_Vector_ID:
        s.assign("CT_Epetra_Vector_ID");
        break;
    case CT_Epetra_Export_ID:
        s.assign("CT_Epetra_Export_ID");
        break;
    case CT_Epetra_Map_ID:
        s.assign("CT_Epetra_Map_ID");
        break;
    case CT_Epetra_BlockMap_ID:
        s.assign("CT_Epetra_BlockMap_ID");
        break;
    case CT_Epetra_Import_ID:
        s.assign("CT_Epetra_Import_ID");
        break;
    case CT_Epetra_Time_ID:
        s.assign("CT_Epetra_Time_ID");
        break;
    case CT_Epetra_JadMatrix_ID:
        s.assign("CT_Epetra_JadMatrix_ID");
        break;
    case CT_Epetra_LinearProblem_ID:
        s.assign("CT_Epetra_LinearProblem_ID");
        break;
    case CT_Epetra_LAPACK_ID:
        s.assign("CT_Epetra_LAPACK_ID");
        break;
    case CT_Teuchos_CommandLineProcessor_ID:
        s.assign("CT_Teuchos_CommandLineProcessor_ID");
        break;
    case CT_Teuchos_ParameterList_ID:
        s.assign("CT_Teuchos_ParameterList_ID");
        break;
    case CT_Teuchos_ParameterEntry_ID:
        s.assign("CT_Teuchos_ParameterEntry_ID");
        break;
    case CT_Teuchos_any_ID:
        s.assign("CT_Teuchos_any_ID");
        break;
    case CT_Amesos_BaseSolver_ID:
        s.assign("CT_Amesos_BaseSolver_ID");
        break;
    case CT_Amesos_ID:
        s.assign("CT_Amesos_ID");
        break;
    case CT_Epetra_FECrsMatrix_ID:
        s.assign("CT_Epetra_FECrsMatrix_ID");
        break;
    case CT_Epetra_IntSerialDenseVector_ID:
        s.assign("CT_Epetra_IntSerialDenseVector_ID");
        break;
    case CT_Epetra_SerialDenseMatrix_ID:
        s.assign("CT_Epetra_SerialDenseMatrix_ID");
        break;
    case CT_AztecOO_ID:
        s.assign("CT_AztecOO_ID");
        break;
    case CT_AztecOO_StatusTest_ID:
        s.assign("CT_AztecOO_StatusTest_ID");
        break;
    case CT_AztecOO_StatusTestCombo_ID:
        s.assign("CT_AztecOO_StatusTestCombo_ID");
        break;
    case CT_AztecOO_StatusTestMaxIters_ID:
        s.assign("CT_AztecOO_StatusTestMaxIters_ID");
        break;
    case CT_AztecOO_StatusTestResNorm_ID:
        s.assign("CT_AztecOO_StatusTestResNorm_ID");
        break;
    case CT_Ifpack_ID:
        s.assign("CT_Ifpack_ID");
        break;
    case CT_Ifpack_Preconditioner_ID:
        s.assign("CT_Ifpack_Preconditioner_ID");
        break;
    case CT_Epetra_SerialDenseVector_ID:
        s.assign("CT_Epetra_SerialDenseVector_ID");
        break;
    case CT_Pliris_ID:
        s.assign("CT_Pliris_ID");
        break;
    default:
        s.assign("(unrecognized)");
        break;
    }

    return s;
}


} // namespace CTrilinos

