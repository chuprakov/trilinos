
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


#ifdef HAVE_CTRILINOS_GALERI


#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"
#include "CEpetra_LinearProblem_Cpp.hpp"
#include "CEpetra_RowMatrix_Cpp.hpp"
#include "CEpetra_CrsMatrix_Cpp.hpp"
#include "CGaleri_Utils_Cpp.hpp"
#include "CGaleri_Utils.h"
#include "Galeri_Utils.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;



} // namespace


//
// Definitions from CGaleri_Utils.h
//


extern "C" {


CT_Epetra_MultiVector_ID_t Galeri_Utils_CreateCartesianCoordinates ( 
  const char CoordType[], CT_Epetra_BlockMap_ID_t BlockMapID, 
  CT_Teuchos_ParameterList_ID_t ListID )
{
    return CEpetra::storeMultiVector(
        Galeri::CreateCartesianCoordinates(std::string(CoordType), 
        CEpetra::getConstBlockMap(BlockMapID).getRawPtr(), 
        *CTeuchos::getParameterList(ListID)));
}

void Galeri_Utils_Solve_LinearProblem ( 
  CT_Epetra_LinearProblem_ID_t ProblemID )
{
    Galeri::Solve(*CEpetra::getConstLinearProblem(ProblemID));
}

void Galeri_Utils_Solve_Matrix ( 
  CT_Epetra_RowMatrix_ID_t MatrixID, 
  CT_Epetra_MultiVector_ID_t LHSID, 
  CT_Epetra_MultiVector_ID_t RHSID )
{
    Galeri::Solve(CEpetra::getConstRowMatrix(MatrixID).getRawPtr(), 
        CEpetra::getConstMultiVector(LHSID).getRawPtr(), 
        CEpetra::getConstMultiVector(RHSID).getRawPtr());
}

double Galeri_Utils_ComputeNorm ( 
  CT_Epetra_MultiVector_ID_t LHSID, 
  CT_Epetra_MultiVector_ID_t RHSID )
{
    return Galeri::ComputeNorm(CEpetra::getConstMultiVector(
        LHSID).getRawPtr(), CEpetra::getConstMultiVector(
        RHSID).getRawPtr());
}

double Galeri_Utils_ComputeNorm_Matrix ( 
  CT_Epetra_RowMatrix_ID_t AID, CT_Epetra_MultiVector_ID_t LHSID, 
  CT_Epetra_MultiVector_ID_t RHSID )
{
    return Galeri::ComputeNorm(CEpetra::getConstRowMatrix(
        AID).getRawPtr(), CEpetra::getConstMultiVector(
        LHSID).getRawPtr(), CEpetra::getConstMultiVector(
        RHSID).getRawPtr());
}

const char * Galeri_Utils_toString_Int ( int x )
{
    return Galeri::toString(x).c_str();
}

const char * Galeri_Utils_toString_UInt ( unsigned int x )
{
    return Galeri::toString(x).c_str();
}

const char * Galeri_Utils_toString_Double ( double x )
{
    return Galeri::toString(x).c_str();
}

void Galeri_Utils_GetNeighboursCartesian2d ( 
  const int i, const int nx, const int ny, int * left, int * right, 
  int * lower, int * upper )
{
    Galeri::GetNeighboursCartesian2d(i, nx, ny, *left, *right, 
        *lower, *upper);
}

void Galeri_Utils_GetNeighboursCartesian2d_Both ( 
  const int i, const int nx, const int ny, int * left, int * right, 
  int * lower, int * upper, int * left2, int * right2, int * lower2, 
  int * upper2 )
{
    Galeri::GetNeighboursCartesian2d(i, nx, ny, *left, *right, 
        *lower, *upper, *left2, *right2, *lower2, *upper2);
}

void Galeri_Utils_GetNeighboursCartesian3d ( 
  const int i, const int nx, const int ny, const int nz, int * left, 
  int * right, int * lower, int * upper, int * below, int * above )
{
    Galeri::GetNeighboursCartesian3d(i, nx, ny, nz, *left, *right, 
        *lower, *upper, *below, *above);
}

void Galeri_Utils_PrintStencil2D ( 
  CT_Epetra_CrsMatrix_ID_t MatrixID, const int nx, const int ny, 
  int GID )
{
    Galeri::PrintStencil2D(CEpetra::getConstCrsMatrix(
        MatrixID).getRawPtr(), nx, ny, GID);
}


} // extern "C"


//
// Definitions from CGaleri_Utils_Cpp.hpp
//




#endif /* HAVE_CTRILINOS_GALERI */


