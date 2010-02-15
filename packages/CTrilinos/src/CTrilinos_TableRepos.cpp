
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


/*! @file CTrilinos_TableRepos.cpp
 * @brief Central table repository for CTrilinos. */


#include "CTrilinos_TableRepos.hpp"


namespace CTrilinos {


void TableRepos::getTable(Table<Epetra_Distributor> *& tab)
{
    tab = &tab_Epetra_Distributor;
}
void TableRepos::getTable(Table<Epetra_SerialComm> *& tab)
{
    tab = &tab_Epetra_SerialComm;
}
void TableRepos::getTable(Table<Epetra_BLAS> *& tab)
{
    tab = &tab_Epetra_BLAS;
}
void TableRepos::getTable(Table<Epetra_Comm> *& tab)
{
    tab = &tab_Epetra_Comm;
}
void TableRepos::getTable(Table<Epetra_Operator> *& tab)
{
    tab = &tab_Epetra_Operator;
}
void TableRepos::getTable(Table<Epetra_MultiVector> *& tab)
{
    tab = &tab_Epetra_MultiVector;
}
void TableRepos::getTable(Table<Epetra_OffsetIndex> *& tab)
{
    tab = &tab_Epetra_OffsetIndex;
}
void TableRepos::getTable(Table<Epetra_Object> *& tab)
{
    tab = &tab_Epetra_Object;
}
void TableRepos::getTable(Table<Epetra_RowMatrix> *& tab)
{
    tab = &tab_Epetra_RowMatrix;
}
void TableRepos::getTable(Table<Epetra_CompObject> *& tab)
{
    tab = &tab_Epetra_CompObject;
}
void TableRepos::getTable(Table<Epetra_Directory> *& tab)
{
    tab = &tab_Epetra_Directory;
}
void TableRepos::getTable(Table<Epetra_Flops> *& tab)
{
    tab = &tab_Epetra_Flops;
}
void TableRepos::getTable(Table<Epetra_SrcDistObject> *& tab)
{
    tab = &tab_Epetra_SrcDistObject;
}
#ifdef HAVE_MPI
void TableRepos::getTable(Table<Epetra_MpiComm> *& tab)
{
    tab = &tab_Epetra_MpiComm;
}
#endif /* HAVE_MPI */
void TableRepos::getTable(Table<Epetra_CrsMatrix> *& tab)
{
    tab = &tab_Epetra_CrsMatrix;
}
void TableRepos::getTable(Table<Epetra_CrsGraph> *& tab)
{
    tab = &tab_Epetra_CrsGraph;
}
void TableRepos::getTable(Table<Epetra_DistObject> *& tab)
{
    tab = &tab_Epetra_DistObject;
}
void TableRepos::getTable(Table<Epetra_Vector> *& tab)
{
    tab = &tab_Epetra_Vector;
}
void TableRepos::getTable(Table<Epetra_Export> *& tab)
{
    tab = &tab_Epetra_Export;
}
void TableRepos::getTable(Table<Epetra_Map> *& tab)
{
    tab = &tab_Epetra_Map;
}
void TableRepos::getTable(Table<Epetra_BlockMap> *& tab)
{
    tab = &tab_Epetra_BlockMap;
}
void TableRepos::getTable(Table<Epetra_Import> *& tab)
{
    tab = &tab_Epetra_Import;
}
void TableRepos::getTable(Table<Epetra_Time> *& tab)
{
    tab = &tab_Epetra_Time;
}
void TableRepos::getTable(Table<Epetra_JadMatrix> *& tab)
{
    tab = &tab_Epetra_JadMatrix;
}
void TableRepos::getTable(Table<Epetra_LinearProblem> *& tab)
{
    tab = &tab_Epetra_LinearProblem;
}
void TableRepos::getTable(Table<Epetra_LAPACK> *& tab)
{
    tab = &tab_Epetra_LAPACK;
}
void TableRepos::getTable(Table<Teuchos::ParameterList> *& tab)
{
    tab = &tab_Teuchos_ParameterList;
}
#ifdef HAVE_CTRILINOS_AMESOS
void TableRepos::getTable(Table<Amesos_BaseSolver> *& tab)
{
    tab = &tab_Amesos_BaseSolver;
}
#endif /* HAVE_CTRILINOS_AMESOS */
void TableRepos::getTable(Table<Epetra_FECrsMatrix> *& tab)
{
    tab = &tab_Epetra_FECrsMatrix;
}
void TableRepos::getTable(Table<Epetra_IntSerialDenseVector> *& tab)
{
    tab = &tab_Epetra_IntSerialDenseVector;
}
void TableRepos::getTable(Table<Epetra_SerialDenseMatrix> *& tab)
{
    tab = &tab_Epetra_SerialDenseMatrix;
}
#ifdef HAVE_CTRILINOS_AZTECOO
void TableRepos::getTable(Table<AztecOO_StatusTest> *& tab)
{
    tab = &tab_AztecOO_StatusTest;
}
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
void TableRepos::getTable(Table<AztecOO_StatusTestCombo> *& tab)
{
    tab = &tab_AztecOO_StatusTestCombo;
}
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
void TableRepos::getTable(Table<AztecOO_StatusTestMaxIters> *& tab)
{
    tab = &tab_AztecOO_StatusTestMaxIters;
}
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
void TableRepos::getTable(Table<AztecOO_StatusTestResNorm> *& tab)
{
    tab = &tab_AztecOO_StatusTestResNorm;
}
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
void TableRepos::getTable(Table<Ifpack_Preconditioner> *& tab)
{
    tab = &tab_Ifpack_Preconditioner;
}
#endif /* HAVE_CTRILINOS_IFPACK */

} // namespace CTrilinos


