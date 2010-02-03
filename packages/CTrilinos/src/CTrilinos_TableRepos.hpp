
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


/*! @file CTrilinos_TableRepos.hpp
 * @brief Central table repository for CTrilinos. */


#ifndef CTRILINOS_TABLEREPOS_HPP
#define CTRILINOS_TABLEREPOS_HPP


#include "CTrilinos_config.h"


#include <string>

#include "CEpetra_Distributor_Cpp.hpp"
#include "CEpetra_SerialComm_Cpp.hpp"
#include "CEpetra_BLAS_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_Operator_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_OffsetIndex_Cpp.hpp"
#include "CEpetra_Object_Cpp.hpp"
#include "CEpetra_RowMatrix_Cpp.hpp"
#include "CEpetra_CompObject_Cpp.hpp"
#include "CEpetra_Directory_Cpp.hpp"
#include "CEpetra_Flops_Cpp.hpp"
#include "CEpetra_SrcDistObject_Cpp.hpp"
#include "CEpetra_MpiComm_Cpp.hpp"
#include "CEpetra_CrsMatrix_Cpp.hpp"
#include "CEpetra_CrsGraph_Cpp.hpp"
#include "CEpetra_DistObject_Cpp.hpp"
#include "CEpetra_Vector_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_Time_Cpp.hpp"
#include "CEpetra_JadMatrix_Cpp.hpp"
#include "CEpetra_LinearProblem_Cpp.hpp"
#include "CEpetra_LAPACK_Cpp.hpp"
#include "CTeuchos_CommandLineProcessor_Cpp.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"
#include "CTeuchos_ParameterEntry_Cpp.hpp"
#include "CTeuchos_any_Cpp.hpp"
#include "CAmesos_BaseSolver_Cpp.hpp"
#include "CAmesos_Cpp.hpp"
#include "CEpetra_FECrsMatrix_Cpp.hpp"
#include "CEpetra_IntSerialDenseVector_Cpp.hpp"
#include "CEpetra_SerialDenseMatrix_Cpp.hpp"
#include "CAztecOO_Cpp.hpp"
#include "CAztecOO_StatusTest_Cpp.hpp"
#include "CAztecOO_StatusTestCombo_Cpp.hpp"
#include "CAztecOO_StatusTestMaxIters_Cpp.hpp"
#include "CAztecOO_StatusTestResNorm_Cpp.hpp"
#include "CIfpack_Cpp.hpp"
#include "CIfpack_Preconditioner_Cpp.hpp"

#include "CTrilinos_enums.h"
#include "CTrilinos_Table.hpp"
#include "CTrilinos_utils_templ.hpp"


#include "Teuchos_RCP.hpp"
namespace CTrilinos {


/*! Central table repository for storing Trilinos objects. */
class TableRepos
{
  public:

    /*! constructor */
    TableRepos();

    /*! destructor */
    ~TableRepos();

    /*! retrieve the object */
    template <class T, typename S>
    const Teuchos::RCP<T> get(S id);

    /*! retrieve the object */
    template <class T, typename S>
    const Teuchos::RCP<const T> getConst(S id);

    /*! store a non-const RCP to object of type T */
    template <class T, typename S>
    S store(T* pobj, bool owned = true);

    /*! store a const RCP to object of type T */
    template <class T, typename S>
    S store(const T* pobj, bool owned = false);

    /*! remove an object from the table and invalidate the id struct */
    template <typename S>
    void remove(S * id);

    /*! dump the table's contents but keep it's properties */
    void purge();

    /*! migrate RCP from one table to another table */
    template <typename S>
    void migrate(S *id, CTrilinos_Table_ID_t tab);

    /*! create an alias for the object in another table */
    template <typename S1, typename S2>
    S1 alias(S2 id, CTrilinos_Table_ID_t tab, bool keepold = true);

  private:
    /*! get the Ctrilinos::Table for the given type */
    template <class T>
    Table<T> & getTable();

    /*! create an alias for the object in another table */
    template <class T>
    CTrilinos_Universal_ID_t do_alias(Table<T> &tab, CTrilinos_Universal_ID_t &aid, bool keepold);

    /*! create an alias for the object in another table */
    template <class T>
    CTrilinos_Universal_ID_t do_alias_const(Table<T> &tab, CTrilinos_Universal_ID_t &aid, bool keepold);

    /*! retrieve the object */
    template <class T>
    const Teuchos::RCP<T> getPoly(CTrilinos_Universal_ID_t aid);

    /*! retrieve the object */
    template <class T>
    const Teuchos::RCP<const T> getConstPoly(CTrilinos_Universal_ID_t aid);

    Table<Epetra_Distributor> tab_Epetra_Distributor;
    Table<Epetra_SerialComm> tab_Epetra_SerialComm;
    Table<Epetra_BLAS> tab_Epetra_BLAS;
    Table<Epetra_Comm> tab_Epetra_Comm;
    Table<Epetra_Operator> tab_Epetra_Operator;
    Table<Epetra_MultiVector> tab_Epetra_MultiVector;
    Table<Epetra_OffsetIndex> tab_Epetra_OffsetIndex;
    Table<Epetra_Object> tab_Epetra_Object;
    Table<Epetra_RowMatrix> tab_Epetra_RowMatrix;
    Table<Epetra_CompObject> tab_Epetra_CompObject;
    Table<Epetra_Directory> tab_Epetra_Directory;
    Table<Epetra_Flops> tab_Epetra_Flops;
    Table<Epetra_SrcDistObject> tab_Epetra_SrcDistObject;
#ifdef HAVE_MPI
    Table<Epetra_MpiComm> tab_Epetra_MpiComm;
#endif /* HAVE_MPI */
    Table<Epetra_CrsMatrix> tab_Epetra_CrsMatrix;
    Table<Epetra_CrsGraph> tab_Epetra_CrsGraph;
    Table<Epetra_DistObject> tab_Epetra_DistObject;
    Table<Epetra_Vector> tab_Epetra_Vector;
    Table<Epetra_Export> tab_Epetra_Export;
    Table<Epetra_Map> tab_Epetra_Map;
    Table<Epetra_BlockMap> tab_Epetra_BlockMap;
    Table<Epetra_Import> tab_Epetra_Import;
    Table<Epetra_Time> tab_Epetra_Time;
    Table<Epetra_JadMatrix> tab_Epetra_JadMatrix;
    Table<Epetra_LinearProblem> tab_Epetra_LinearProblem;
    Table<Epetra_LAPACK> tab_Epetra_LAPACK;
    Table<Teuchos::CommandLineProcessor> tab_Teuchos_CommandLineProcessor;
    Table<Teuchos::ParameterList> tab_Teuchos_ParameterList;
    Table<Teuchos::ParameterEntry> tab_Teuchos_ParameterEntry;
    Table<Teuchos::any> tab_Teuchos_any;
#ifdef HAVE_CTRILINOS_AMESOS
    Table<Amesos_BaseSolver> tab_Amesos_BaseSolver;
#endif /* HAVE_CTRILINOS_AMESOS */
#ifdef HAVE_CTRILINOS_AMESOS
    Table<Amesos> tab_Amesos;
#endif /* HAVE_CTRILINOS_AMESOS */
    Table<Epetra_FECrsMatrix> tab_Epetra_FECrsMatrix;
    Table<Epetra_IntSerialDenseVector> tab_Epetra_IntSerialDenseVector;
    Table<Epetra_SerialDenseMatrix> tab_Epetra_SerialDenseMatrix;
#ifdef HAVE_CTRILINOS_AZTECOO
    Table<AztecOO> tab_AztecOO;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    Table<AztecOO_StatusTest> tab_AztecOO_StatusTest;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    Table<AztecOO_StatusTestCombo> tab_AztecOO_StatusTestCombo;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    Table<AztecOO_StatusTestMaxIters> tab_AztecOO_StatusTestMaxIters;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    Table<AztecOO_StatusTestResNorm> tab_AztecOO_StatusTestResNorm;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
    Table<Ifpack> tab_Ifpack;
#endif /* HAVE_CTRILINOS_IFPACK */
#ifdef HAVE_CTRILINOS_IFPACK
    Table<Ifpack_Preconditioner> tab_Ifpack_Preconditioner;
#endif /* HAVE_CTRILINOS_IFPACK */

    bool call_me_lazy;  /* I was too lazy to deal with the commas in the init list, so... */
};

TableRepos::TableRepos() :
    tab_Epetra_Distributor(CT_Epetra_Distributor_ID),
    tab_Epetra_SerialComm(CT_Epetra_SerialComm_ID),
    tab_Epetra_BLAS(CT_Epetra_BLAS_ID),
    tab_Epetra_Comm(CT_Epetra_Comm_ID),
    tab_Epetra_Operator(CT_Epetra_Operator_ID),
    tab_Epetra_MultiVector(CT_Epetra_MultiVector_ID),
    tab_Epetra_OffsetIndex(CT_Epetra_OffsetIndex_ID),
    tab_Epetra_Object(CT_Epetra_Object_ID),
    tab_Epetra_RowMatrix(CT_Epetra_RowMatrix_ID),
    tab_Epetra_CompObject(CT_Epetra_CompObject_ID),
    tab_Epetra_Directory(CT_Epetra_Directory_ID),
    tab_Epetra_Flops(CT_Epetra_Flops_ID),
    tab_Epetra_SrcDistObject(CT_Epetra_SrcDistObject_ID),
#ifdef HAVE_MPI
    tab_Epetra_MpiComm(CT_Epetra_MpiComm_ID),
#endif /* HAVE_MPI */
    tab_Epetra_CrsMatrix(CT_Epetra_CrsMatrix_ID),
    tab_Epetra_CrsGraph(CT_Epetra_CrsGraph_ID),
    tab_Epetra_DistObject(CT_Epetra_DistObject_ID),
    tab_Epetra_Vector(CT_Epetra_Vector_ID),
    tab_Epetra_Export(CT_Epetra_Export_ID),
    tab_Epetra_Map(CT_Epetra_Map_ID),
    tab_Epetra_BlockMap(CT_Epetra_BlockMap_ID),
    tab_Epetra_Import(CT_Epetra_Import_ID),
    tab_Epetra_Time(CT_Epetra_Time_ID),
    tab_Epetra_JadMatrix(CT_Epetra_JadMatrix_ID),
    tab_Epetra_LinearProblem(CT_Epetra_LinearProblem_ID),
    tab_Epetra_LAPACK(CT_Epetra_LAPACK_ID),
    tab_Teuchos_CommandLineProcessor(CT_Teuchos_CommandLineProcessor_ID),
    tab_Teuchos_ParameterList(CT_Teuchos_ParameterList_ID),
    tab_Teuchos_ParameterEntry(CT_Teuchos_ParameterEntry_ID),
    tab_Teuchos_any(CT_Teuchos_any_ID),
#ifdef HAVE_CTRILINOS_AMESOS
    tab_Amesos_BaseSolver(CT_Amesos_BaseSolver_ID),
#endif /* HAVE_CTRILINOS_AMESOS */
#ifdef HAVE_CTRILINOS_AMESOS
    tab_Amesos(CT_Amesos_ID),
#endif /* HAVE_CTRILINOS_AMESOS */
    tab_Epetra_FECrsMatrix(CT_Epetra_FECrsMatrix_ID),
    tab_Epetra_IntSerialDenseVector(CT_Epetra_IntSerialDenseVector_ID),
    tab_Epetra_SerialDenseMatrix(CT_Epetra_SerialDenseMatrix_ID),
#ifdef HAVE_CTRILINOS_AZTECOO
    tab_AztecOO(CT_AztecOO_ID),
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    tab_AztecOO_StatusTest(CT_AztecOO_StatusTest_ID),
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    tab_AztecOO_StatusTestCombo(CT_AztecOO_StatusTestCombo_ID),
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    tab_AztecOO_StatusTestMaxIters(CT_AztecOO_StatusTestMaxIters_ID),
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    tab_AztecOO_StatusTestResNorm(CT_AztecOO_StatusTestResNorm_ID),
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
    tab_Ifpack(CT_Ifpack_ID),
#endif /* HAVE_CTRILINOS_IFPACK */
#ifdef HAVE_CTRILINOS_IFPACK
    tab_Ifpack_Preconditioner(CT_Ifpack_Preconditioner_ID),
#endif /* HAVE_CTRILINOS_IFPACK */
    call_me_lazy(true)
{
}

TableRepos::~TableRepos()
{
}

template <class T, typename S>
S TableRepos::store(T* pobj, bool owned)
{
    return concreteType<S>(getTable<T>().store(pobj, owned));
}

template <class T, typename S>
S TableRepos::store(const T* pobj, bool owned)
{
    return concreteType<S>(getTable<T>().store(pobj, owned));
}

template <class T>
Table<T> & TableRepos::getTable()
{
    Teuchos::RCP<T> dummy = Teuchos::null;
    /* I know this is bad, but it's for catching COMPILE-TIME errors! */
    return dummy->No_Function_For_This_Table_Type;
}

template <class T, typename S>
const Teuchos::RCP<T> TableRepos::get(S id)
{
    if (id.is_const)
        throw std::string("const violation");

    CTrilinos_Universal_ID_t aid = abstractType<S>(id);

    /* Shortcut if stored in the most obvious table */
    Table<T> &tab_maybe = getTable<T>();
    if (tab_maybe.isType(aid.table)) {
        return tab_maybe.get<T>(aid);
    } else {
        return getPoly<T>(aid);
    }
}

template <class T>
const Teuchos::RCP<T> TableRepos::getPoly(CTrilinos_Universal_ID_t aid)
{
    switch (aid.table) {
        case CT_Epetra_Distributor_ID:
            return tab_Epetra_Distributor.get<T>(aid);
            break;
        case CT_Epetra_SerialComm_ID:
            return tab_Epetra_SerialComm.get<T>(aid);
            break;
        case CT_Epetra_BLAS_ID:
            return tab_Epetra_BLAS.get<T>(aid);
            break;
        case CT_Epetra_Comm_ID:
            return tab_Epetra_Comm.get<T>(aid);
            break;
        case CT_Epetra_Operator_ID:
            return tab_Epetra_Operator.get<T>(aid);
            break;
        case CT_Epetra_MultiVector_ID:
            return tab_Epetra_MultiVector.get<T>(aid);
            break;
        case CT_Epetra_OffsetIndex_ID:
            return tab_Epetra_OffsetIndex.get<T>(aid);
            break;
        case CT_Epetra_Object_ID:
            return tab_Epetra_Object.get<T>(aid);
            break;
        case CT_Epetra_RowMatrix_ID:
            return tab_Epetra_RowMatrix.get<T>(aid);
            break;
        case CT_Epetra_CompObject_ID:
            return tab_Epetra_CompObject.get<T>(aid);
            break;
        case CT_Epetra_Directory_ID:
            return tab_Epetra_Directory.get<T>(aid);
            break;
        case CT_Epetra_Flops_ID:
            return tab_Epetra_Flops.get<T>(aid);
            break;
        case CT_Epetra_SrcDistObject_ID:
            return tab_Epetra_SrcDistObject.get<T>(aid);
            break;
#ifdef HAVE_MPI
        case CT_Epetra_MpiComm_ID:
            return tab_Epetra_MpiComm.get<T>(aid);
            break;
#endif /* HAVE_MPI */
        case CT_Epetra_CrsMatrix_ID:
            return tab_Epetra_CrsMatrix.get<T>(aid);
            break;
        case CT_Epetra_CrsGraph_ID:
            return tab_Epetra_CrsGraph.get<T>(aid);
            break;
        case CT_Epetra_DistObject_ID:
            return tab_Epetra_DistObject.get<T>(aid);
            break;
        case CT_Epetra_Vector_ID:
            return tab_Epetra_Vector.get<T>(aid);
            break;
        case CT_Epetra_Export_ID:
            return tab_Epetra_Export.get<T>(aid);
            break;
        case CT_Epetra_Map_ID:
            return tab_Epetra_Map.get<T>(aid);
            break;
        case CT_Epetra_BlockMap_ID:
            return tab_Epetra_BlockMap.get<T>(aid);
            break;
        case CT_Epetra_Import_ID:
            return tab_Epetra_Import.get<T>(aid);
            break;
        case CT_Epetra_Time_ID:
            return tab_Epetra_Time.get<T>(aid);
            break;
        case CT_Epetra_JadMatrix_ID:
            return tab_Epetra_JadMatrix.get<T>(aid);
            break;
        case CT_Epetra_LinearProblem_ID:
            return tab_Epetra_LinearProblem.get<T>(aid);
            break;
        case CT_Epetra_LAPACK_ID:
            return tab_Epetra_LAPACK.get<T>(aid);
            break;
        case CT_Teuchos_ParameterList_ID:
            return tab_Teuchos_ParameterList.get<T>(aid);
            break;
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_BaseSolver_ID:
            return tab_Amesos_BaseSolver.get<T>(aid);
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
        case CT_Epetra_FECrsMatrix_ID:
            return tab_Epetra_FECrsMatrix.get<T>(aid);
            break;
        case CT_Epetra_IntSerialDenseVector_ID:
            return tab_Epetra_IntSerialDenseVector.get<T>(aid);
            break;
        case CT_Epetra_SerialDenseMatrix_ID:
            return tab_Epetra_SerialDenseMatrix.get<T>(aid);
            break;
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTest_ID:
            return tab_AztecOO_StatusTest.get<T>(aid);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestCombo_ID:
            return tab_AztecOO_StatusTestCombo.get<T>(aid);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestMaxIters_ID:
            return tab_AztecOO_StatusTestMaxIters.get<T>(aid);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestResNorm_ID:
            return tab_AztecOO_StatusTestResNorm.get<T>(aid);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
        case CT_Ifpack_Preconditioner_ID:
            return tab_Ifpack_Preconditioner.get<T>(aid);
            break;
#endif /* HAVE_CTRILINOS_IFPACK */
        default:
            throw std::string("invalid table id");
            break;
    }

    return Teuchos::null;
}

template <class T, typename S>
const Teuchos::RCP<const T> TableRepos::getConst(S id)
{
    if (!id.is_const)
        return get<T,S>(id);

    CTrilinos_Universal_ID_t aid = abstractType<S>(id);

    /* Shortcut if stored in the most obvious table */
    Table<T> &tab_maybe = getTable<T>();
    if (tab_maybe.isType(aid.table)) {
        return tab_maybe.getConst<T>(aid);
    } else {
        return getConstPoly<T>(aid);
    }
}

template <class T>
const Teuchos::RCP<const T> TableRepos::getConstPoly(CTrilinos_Universal_ID_t aid)
{
    switch (aid.table) {
        case CT_Epetra_Distributor_ID:
            return tab_Epetra_Distributor.getConst<T>(aid);
            break;
        case CT_Epetra_SerialComm_ID:
            return tab_Epetra_SerialComm.getConst<T>(aid);
            break;
        case CT_Epetra_BLAS_ID:
            return tab_Epetra_BLAS.getConst<T>(aid);
            break;
        case CT_Epetra_Comm_ID:
            return tab_Epetra_Comm.getConst<T>(aid);
            break;
        case CT_Epetra_Operator_ID:
            return tab_Epetra_Operator.getConst<T>(aid);
            break;
        case CT_Epetra_MultiVector_ID:
            return tab_Epetra_MultiVector.getConst<T>(aid);
            break;
        case CT_Epetra_OffsetIndex_ID:
            return tab_Epetra_OffsetIndex.getConst<T>(aid);
            break;
        case CT_Epetra_Object_ID:
            return tab_Epetra_Object.getConst<T>(aid);
            break;
        case CT_Epetra_RowMatrix_ID:
            return tab_Epetra_RowMatrix.getConst<T>(aid);
            break;
        case CT_Epetra_CompObject_ID:
            return tab_Epetra_CompObject.getConst<T>(aid);
            break;
        case CT_Epetra_Directory_ID:
            return tab_Epetra_Directory.getConst<T>(aid);
            break;
        case CT_Epetra_Flops_ID:
            return tab_Epetra_Flops.getConst<T>(aid);
            break;
        case CT_Epetra_SrcDistObject_ID:
            return tab_Epetra_SrcDistObject.getConst<T>(aid);
            break;
#ifdef HAVE_MPI
        case CT_Epetra_MpiComm_ID:
            return tab_Epetra_MpiComm.getConst<T>(aid);
            break;
#endif /* HAVE_MPI */
        case CT_Epetra_CrsMatrix_ID:
            return tab_Epetra_CrsMatrix.getConst<T>(aid);
            break;
        case CT_Epetra_CrsGraph_ID:
            return tab_Epetra_CrsGraph.getConst<T>(aid);
            break;
        case CT_Epetra_DistObject_ID:
            return tab_Epetra_DistObject.getConst<T>(aid);
            break;
        case CT_Epetra_Vector_ID:
            return tab_Epetra_Vector.getConst<T>(aid);
            break;
        case CT_Epetra_Export_ID:
            return tab_Epetra_Export.getConst<T>(aid);
            break;
        case CT_Epetra_Map_ID:
            return tab_Epetra_Map.getConst<T>(aid);
            break;
        case CT_Epetra_BlockMap_ID:
            return tab_Epetra_BlockMap.getConst<T>(aid);
            break;
        case CT_Epetra_Import_ID:
            return tab_Epetra_Import.getConst<T>(aid);
            break;
        case CT_Epetra_Time_ID:
            return tab_Epetra_Time.getConst<T>(aid);
            break;
        case CT_Epetra_JadMatrix_ID:
            return tab_Epetra_JadMatrix.getConst<T>(aid);
            break;
        case CT_Epetra_LinearProblem_ID:
            return tab_Epetra_LinearProblem.getConst<T>(aid);
            break;
        case CT_Epetra_LAPACK_ID:
            return tab_Epetra_LAPACK.getConst<T>(aid);
            break;
        case CT_Teuchos_ParameterList_ID:
            return tab_Teuchos_ParameterList.getConst<T>(aid);
            break;
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_BaseSolver_ID:
            return tab_Amesos_BaseSolver.getConst<T>(aid);
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
        case CT_Epetra_FECrsMatrix_ID:
            return tab_Epetra_FECrsMatrix.getConst<T>(aid);
            break;
        case CT_Epetra_IntSerialDenseVector_ID:
            return tab_Epetra_IntSerialDenseVector.getConst<T>(aid);
            break;
        case CT_Epetra_SerialDenseMatrix_ID:
            return tab_Epetra_SerialDenseMatrix.getConst<T>(aid);
            break;
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTest_ID:
            return tab_AztecOO_StatusTest.getConst<T>(aid);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestCombo_ID:
            return tab_AztecOO_StatusTestCombo.getConst<T>(aid);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestMaxIters_ID:
            return tab_AztecOO_StatusTestMaxIters.getConst<T>(aid);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestResNorm_ID:
            return tab_AztecOO_StatusTestResNorm.getConst<T>(aid);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
        case CT_Ifpack_Preconditioner_ID:
            return tab_Ifpack_Preconditioner.getConst<T>(aid);
            break;
#endif /* HAVE_CTRILINOS_IFPACK */
        default:
            throw std::string("invalid table id");
            break;
    }

    return Teuchos::null;
}

template <>
Table<Epetra_Distributor> & TableRepos::getTable<Epetra_Distributor>()
{
    return tab_Epetra_Distributor;
}

template <>
Table<Epetra_SerialComm> & TableRepos::getTable<Epetra_SerialComm>()
{
    return tab_Epetra_SerialComm;
}

template <>
Table<Epetra_BLAS> & TableRepos::getTable<Epetra_BLAS>()
{
    return tab_Epetra_BLAS;
}

template <>
Table<Epetra_Comm> & TableRepos::getTable<Epetra_Comm>()
{
    return tab_Epetra_Comm;
}

template <>
Table<Epetra_Operator> & TableRepos::getTable<Epetra_Operator>()
{
    return tab_Epetra_Operator;
}

template <>
Table<Epetra_MultiVector> & TableRepos::getTable<Epetra_MultiVector>()
{
    return tab_Epetra_MultiVector;
}

template <>
Table<Epetra_OffsetIndex> & TableRepos::getTable<Epetra_OffsetIndex>()
{
    return tab_Epetra_OffsetIndex;
}

template <>
Table<Epetra_Object> & TableRepos::getTable<Epetra_Object>()
{
    return tab_Epetra_Object;
}

template <>
Table<Epetra_RowMatrix> & TableRepos::getTable<Epetra_RowMatrix>()
{
    return tab_Epetra_RowMatrix;
}

template <>
Table<Epetra_CompObject> & TableRepos::getTable<Epetra_CompObject>()
{
    return tab_Epetra_CompObject;
}

template <>
Table<Epetra_Directory> & TableRepos::getTable<Epetra_Directory>()
{
    return tab_Epetra_Directory;
}

template <>
Table<Epetra_Flops> & TableRepos::getTable<Epetra_Flops>()
{
    return tab_Epetra_Flops;
}

template <>
Table<Epetra_SrcDistObject> & TableRepos::getTable<Epetra_SrcDistObject>()
{
    return tab_Epetra_SrcDistObject;
}

#ifdef HAVE_MPI
template <>
Table<Epetra_MpiComm> & TableRepos::getTable<Epetra_MpiComm>()
{
    return tab_Epetra_MpiComm;
}

#endif /* HAVE_MPI */
template <>
Table<Epetra_CrsMatrix> & TableRepos::getTable<Epetra_CrsMatrix>()
{
    return tab_Epetra_CrsMatrix;
}

template <>
Table<Epetra_CrsGraph> & TableRepos::getTable<Epetra_CrsGraph>()
{
    return tab_Epetra_CrsGraph;
}

template <>
Table<Epetra_DistObject> & TableRepos::getTable<Epetra_DistObject>()
{
    return tab_Epetra_DistObject;
}

template <>
Table<Epetra_Vector> & TableRepos::getTable<Epetra_Vector>()
{
    return tab_Epetra_Vector;
}

template <>
Table<Epetra_Export> & TableRepos::getTable<Epetra_Export>()
{
    return tab_Epetra_Export;
}

template <>
Table<Epetra_Map> & TableRepos::getTable<Epetra_Map>()
{
    return tab_Epetra_Map;
}

template <>
Table<Epetra_BlockMap> & TableRepos::getTable<Epetra_BlockMap>()
{
    return tab_Epetra_BlockMap;
}

template <>
Table<Epetra_Import> & TableRepos::getTable<Epetra_Import>()
{
    return tab_Epetra_Import;
}

template <>
Table<Epetra_Time> & TableRepos::getTable<Epetra_Time>()
{
    return tab_Epetra_Time;
}

template <>
Table<Epetra_JadMatrix> & TableRepos::getTable<Epetra_JadMatrix>()
{
    return tab_Epetra_JadMatrix;
}

template <>
Table<Epetra_LinearProblem> & TableRepos::getTable<Epetra_LinearProblem>()
{
    return tab_Epetra_LinearProblem;
}

template <>
Table<Epetra_LAPACK> & TableRepos::getTable<Epetra_LAPACK>()
{
    return tab_Epetra_LAPACK;
}

template <>
Table<Teuchos::CommandLineProcessor> & TableRepos::getTable<Teuchos::CommandLineProcessor>()
{
    return tab_Teuchos_CommandLineProcessor;
}

template <>
Table<Teuchos::ParameterList> & TableRepos::getTable<Teuchos::ParameterList>()
{
    return tab_Teuchos_ParameterList;
}

template <>
Table<Teuchos::ParameterEntry> & TableRepos::getTable<Teuchos::ParameterEntry>()
{
    return tab_Teuchos_ParameterEntry;
}

template <>
Table<Teuchos::any> & TableRepos::getTable<Teuchos::any>()
{
    return tab_Teuchos_any;
}

#ifdef HAVE_CTRILINOS_AMESOS
template <>
Table<Amesos_BaseSolver> & TableRepos::getTable<Amesos_BaseSolver>()
{
    return tab_Amesos_BaseSolver;
}

#endif /* HAVE_CTRILINOS_AMESOS */
#ifdef HAVE_CTRILINOS_AMESOS
template <>
Table<Amesos> & TableRepos::getTable<Amesos>()
{
    return tab_Amesos;
}

#endif /* HAVE_CTRILINOS_AMESOS */
template <>
Table<Epetra_FECrsMatrix> & TableRepos::getTable<Epetra_FECrsMatrix>()
{
    return tab_Epetra_FECrsMatrix;
}

template <>
Table<Epetra_IntSerialDenseVector> & TableRepos::getTable<Epetra_IntSerialDenseVector>()
{
    return tab_Epetra_IntSerialDenseVector;
}

template <>
Table<Epetra_SerialDenseMatrix> & TableRepos::getTable<Epetra_SerialDenseMatrix>()
{
    return tab_Epetra_SerialDenseMatrix;
}

#ifdef HAVE_CTRILINOS_AZTECOO
template <>
Table<AztecOO> & TableRepos::getTable<AztecOO>()
{
    return tab_AztecOO;
}

#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
template <>
Table<AztecOO_StatusTest> & TableRepos::getTable<AztecOO_StatusTest>()
{
    return tab_AztecOO_StatusTest;
}

#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
template <>
Table<AztecOO_StatusTestCombo> & TableRepos::getTable<AztecOO_StatusTestCombo>()
{
    return tab_AztecOO_StatusTestCombo;
}

#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
template <>
Table<AztecOO_StatusTestMaxIters> & TableRepos::getTable<AztecOO_StatusTestMaxIters>()
{
    return tab_AztecOO_StatusTestMaxIters;
}

#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
template <>
Table<AztecOO_StatusTestResNorm> & TableRepos::getTable<AztecOO_StatusTestResNorm>()
{
    return tab_AztecOO_StatusTestResNorm;
}

#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
template <>
Table<Ifpack> & TableRepos::getTable<Ifpack>()
{
    return tab_Ifpack;
}

#endif /* HAVE_CTRILINOS_IFPACK */
#ifdef HAVE_CTRILINOS_IFPACK
template <>
Table<Ifpack_Preconditioner> & TableRepos::getTable<Ifpack_Preconditioner>()
{
    return tab_Ifpack_Preconditioner;
}

#endif /* HAVE_CTRILINOS_IFPACK */
template <class T>
CTrilinos_Universal_ID_t TableRepos::do_alias(Table<T> &tab, CTrilinos_Universal_ID_t &aid, bool keepold)
{
    CTrilinos_Universal_ID_t newid;

    switch (aid.table) {
        case CT_Epetra_Distributor_ID:
            newid = tab.alias(tab_Epetra_Distributor.get<T>(aid));
            if (!keepold) tab_Epetra_Distributor.remove(&aid);
            break;
        case CT_Epetra_SerialComm_ID:
            newid = tab.alias(tab_Epetra_SerialComm.get<T>(aid));
            if (!keepold) tab_Epetra_SerialComm.remove(&aid);
            break;
        case CT_Epetra_BLAS_ID:
            newid = tab.alias(tab_Epetra_BLAS.get<T>(aid));
            if (!keepold) tab_Epetra_BLAS.remove(&aid);
            break;
        case CT_Epetra_Comm_ID:
            newid = tab.alias(tab_Epetra_Comm.get<T>(aid));
            if (!keepold) tab_Epetra_Comm.remove(&aid);
            break;
        case CT_Epetra_Operator_ID:
            newid = tab.alias(tab_Epetra_Operator.get<T>(aid));
            if (!keepold) tab_Epetra_Operator.remove(&aid);
            break;
        case CT_Epetra_MultiVector_ID:
            newid = tab.alias(tab_Epetra_MultiVector.get<T>(aid));
            if (!keepold) tab_Epetra_MultiVector.remove(&aid);
            break;
        case CT_Epetra_OffsetIndex_ID:
            newid = tab.alias(tab_Epetra_OffsetIndex.get<T>(aid));
            if (!keepold) tab_Epetra_OffsetIndex.remove(&aid);
            break;
        case CT_Epetra_Object_ID:
            newid = tab.alias(tab_Epetra_Object.get<T>(aid));
            if (!keepold) tab_Epetra_Object.remove(&aid);
            break;
        case CT_Epetra_RowMatrix_ID:
            newid = tab.alias(tab_Epetra_RowMatrix.get<T>(aid));
            if (!keepold) tab_Epetra_RowMatrix.remove(&aid);
            break;
        case CT_Epetra_CompObject_ID:
            newid = tab.alias(tab_Epetra_CompObject.get<T>(aid));
            if (!keepold) tab_Epetra_CompObject.remove(&aid);
            break;
        case CT_Epetra_Directory_ID:
            newid = tab.alias(tab_Epetra_Directory.get<T>(aid));
            if (!keepold) tab_Epetra_Directory.remove(&aid);
            break;
        case CT_Epetra_Flops_ID:
            newid = tab.alias(tab_Epetra_Flops.get<T>(aid));
            if (!keepold) tab_Epetra_Flops.remove(&aid);
            break;
        case CT_Epetra_SrcDistObject_ID:
            newid = tab.alias(tab_Epetra_SrcDistObject.get<T>(aid));
            if (!keepold) tab_Epetra_SrcDistObject.remove(&aid);
            break;
#ifdef HAVE_MPI
        case CT_Epetra_MpiComm_ID:
            newid = tab.alias(tab_Epetra_MpiComm.get<T>(aid));
            if (!keepold) tab_Epetra_MpiComm.remove(&aid);
            break;
#endif /* HAVE_MPI */
        case CT_Epetra_CrsMatrix_ID:
            newid = tab.alias(tab_Epetra_CrsMatrix.get<T>(aid));
            if (!keepold) tab_Epetra_CrsMatrix.remove(&aid);
            break;
        case CT_Epetra_CrsGraph_ID:
            newid = tab.alias(tab_Epetra_CrsGraph.get<T>(aid));
            if (!keepold) tab_Epetra_CrsGraph.remove(&aid);
            break;
        case CT_Epetra_DistObject_ID:
            newid = tab.alias(tab_Epetra_DistObject.get<T>(aid));
            if (!keepold) tab_Epetra_DistObject.remove(&aid);
            break;
        case CT_Epetra_Vector_ID:
            newid = tab.alias(tab_Epetra_Vector.get<T>(aid));
            if (!keepold) tab_Epetra_Vector.remove(&aid);
            break;
        case CT_Epetra_Export_ID:
            newid = tab.alias(tab_Epetra_Export.get<T>(aid));
            if (!keepold) tab_Epetra_Export.remove(&aid);
            break;
        case CT_Epetra_Map_ID:
            newid = tab.alias(tab_Epetra_Map.get<T>(aid));
            if (!keepold) tab_Epetra_Map.remove(&aid);
            break;
        case CT_Epetra_BlockMap_ID:
            newid = tab.alias(tab_Epetra_BlockMap.get<T>(aid));
            if (!keepold) tab_Epetra_BlockMap.remove(&aid);
            break;
        case CT_Epetra_Import_ID:
            newid = tab.alias(tab_Epetra_Import.get<T>(aid));
            if (!keepold) tab_Epetra_Import.remove(&aid);
            break;
        case CT_Epetra_Time_ID:
            newid = tab.alias(tab_Epetra_Time.get<T>(aid));
            if (!keepold) tab_Epetra_Time.remove(&aid);
            break;
        case CT_Epetra_JadMatrix_ID:
            newid = tab.alias(tab_Epetra_JadMatrix.get<T>(aid));
            if (!keepold) tab_Epetra_JadMatrix.remove(&aid);
            break;
        case CT_Epetra_LinearProblem_ID:
            newid = tab.alias(tab_Epetra_LinearProblem.get<T>(aid));
            if (!keepold) tab_Epetra_LinearProblem.remove(&aid);
            break;
        case CT_Epetra_LAPACK_ID:
            newid = tab.alias(tab_Epetra_LAPACK.get<T>(aid));
            if (!keepold) tab_Epetra_LAPACK.remove(&aid);
            break;
        case CT_Teuchos_ParameterList_ID:
            newid = tab.alias(tab_Teuchos_ParameterList.get<T>(aid));
            if (!keepold) tab_Teuchos_ParameterList.remove(&aid);
            break;
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_BaseSolver_ID:
            newid = tab.alias(tab_Amesos_BaseSolver.get<T>(aid));
            if (!keepold) tab_Amesos_BaseSolver.remove(&aid);
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
        case CT_Epetra_FECrsMatrix_ID:
            newid = tab.alias(tab_Epetra_FECrsMatrix.get<T>(aid));
            if (!keepold) tab_Epetra_FECrsMatrix.remove(&aid);
            break;
        case CT_Epetra_IntSerialDenseVector_ID:
            newid = tab.alias(tab_Epetra_IntSerialDenseVector.get<T>(aid));
            if (!keepold) tab_Epetra_IntSerialDenseVector.remove(&aid);
            break;
        case CT_Epetra_SerialDenseMatrix_ID:
            newid = tab.alias(tab_Epetra_SerialDenseMatrix.get<T>(aid));
            if (!keepold) tab_Epetra_SerialDenseMatrix.remove(&aid);
            break;
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTest_ID:
            newid = tab.alias(tab_AztecOO_StatusTest.get<T>(aid));
            if (!keepold) tab_AztecOO_StatusTest.remove(&aid);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestCombo_ID:
            newid = tab.alias(tab_AztecOO_StatusTestCombo.get<T>(aid));
            if (!keepold) tab_AztecOO_StatusTestCombo.remove(&aid);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestMaxIters_ID:
            newid = tab.alias(tab_AztecOO_StatusTestMaxIters.get<T>(aid));
            if (!keepold) tab_AztecOO_StatusTestMaxIters.remove(&aid);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestResNorm_ID:
            newid = tab.alias(tab_AztecOO_StatusTestResNorm.get<T>(aid));
            if (!keepold) tab_AztecOO_StatusTestResNorm.remove(&aid);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
        case CT_Ifpack_Preconditioner_ID:
            newid = tab.alias(tab_Ifpack_Preconditioner.get<T>(aid));
            if (!keepold) tab_Ifpack_Preconditioner.remove(&aid);
            break;
#endif /* HAVE_CTRILINOS_IFPACK */
        default:
            newid.table = CT_Invalid_ID;
            newid.index = -1;
            newid.is_const = aid.is_const;
            throw std::string("invalid table id or non-polymorphic class");
            break;
    }

    return newid;
}

template <class T>
CTrilinos_Universal_ID_t TableRepos::do_alias_const(Table<T> &tab, CTrilinos_Universal_ID_t &aid, bool keepold)
{
    CTrilinos_Universal_ID_t newid;

    switch (aid.table) {
        case CT_Epetra_Distributor_ID:
            newid = tab.alias(tab_Epetra_Distributor.getConst<T>(aid));
            if (!keepold) tab_Epetra_Distributor.remove(&aid);
            break;
        case CT_Epetra_SerialComm_ID:
            newid = tab.alias(tab_Epetra_SerialComm.getConst<T>(aid));
            if (!keepold) tab_Epetra_SerialComm.remove(&aid);
            break;
        case CT_Epetra_BLAS_ID:
            newid = tab.alias(tab_Epetra_BLAS.getConst<T>(aid));
            if (!keepold) tab_Epetra_BLAS.remove(&aid);
            break;
        case CT_Epetra_Comm_ID:
            newid = tab.alias(tab_Epetra_Comm.getConst<T>(aid));
            if (!keepold) tab_Epetra_Comm.remove(&aid);
            break;
        case CT_Epetra_Operator_ID:
            newid = tab.alias(tab_Epetra_Operator.getConst<T>(aid));
            if (!keepold) tab_Epetra_Operator.remove(&aid);
            break;
        case CT_Epetra_MultiVector_ID:
            newid = tab.alias(tab_Epetra_MultiVector.getConst<T>(aid));
            if (!keepold) tab_Epetra_MultiVector.remove(&aid);
            break;
        case CT_Epetra_OffsetIndex_ID:
            newid = tab.alias(tab_Epetra_OffsetIndex.getConst<T>(aid));
            if (!keepold) tab_Epetra_OffsetIndex.remove(&aid);
            break;
        case CT_Epetra_Object_ID:
            newid = tab.alias(tab_Epetra_Object.getConst<T>(aid));
            if (!keepold) tab_Epetra_Object.remove(&aid);
            break;
        case CT_Epetra_RowMatrix_ID:
            newid = tab.alias(tab_Epetra_RowMatrix.getConst<T>(aid));
            if (!keepold) tab_Epetra_RowMatrix.remove(&aid);
            break;
        case CT_Epetra_CompObject_ID:
            newid = tab.alias(tab_Epetra_CompObject.getConst<T>(aid));
            if (!keepold) tab_Epetra_CompObject.remove(&aid);
            break;
        case CT_Epetra_Directory_ID:
            newid = tab.alias(tab_Epetra_Directory.getConst<T>(aid));
            if (!keepold) tab_Epetra_Directory.remove(&aid);
            break;
        case CT_Epetra_Flops_ID:
            newid = tab.alias(tab_Epetra_Flops.getConst<T>(aid));
            if (!keepold) tab_Epetra_Flops.remove(&aid);
            break;
        case CT_Epetra_SrcDistObject_ID:
            newid = tab.alias(tab_Epetra_SrcDistObject.getConst<T>(aid));
            if (!keepold) tab_Epetra_SrcDistObject.remove(&aid);
            break;
#ifdef HAVE_MPI
        case CT_Epetra_MpiComm_ID:
            newid = tab.alias(tab_Epetra_MpiComm.getConst<T>(aid));
            if (!keepold) tab_Epetra_MpiComm.remove(&aid);
            break;
#endif /* HAVE_MPI */
        case CT_Epetra_CrsMatrix_ID:
            newid = tab.alias(tab_Epetra_CrsMatrix.getConst<T>(aid));
            if (!keepold) tab_Epetra_CrsMatrix.remove(&aid);
            break;
        case CT_Epetra_CrsGraph_ID:
            newid = tab.alias(tab_Epetra_CrsGraph.getConst<T>(aid));
            if (!keepold) tab_Epetra_CrsGraph.remove(&aid);
            break;
        case CT_Epetra_DistObject_ID:
            newid = tab.alias(tab_Epetra_DistObject.getConst<T>(aid));
            if (!keepold) tab_Epetra_DistObject.remove(&aid);
            break;
        case CT_Epetra_Vector_ID:
            newid = tab.alias(tab_Epetra_Vector.getConst<T>(aid));
            if (!keepold) tab_Epetra_Vector.remove(&aid);
            break;
        case CT_Epetra_Export_ID:
            newid = tab.alias(tab_Epetra_Export.getConst<T>(aid));
            if (!keepold) tab_Epetra_Export.remove(&aid);
            break;
        case CT_Epetra_Map_ID:
            newid = tab.alias(tab_Epetra_Map.getConst<T>(aid));
            if (!keepold) tab_Epetra_Map.remove(&aid);
            break;
        case CT_Epetra_BlockMap_ID:
            newid = tab.alias(tab_Epetra_BlockMap.getConst<T>(aid));
            if (!keepold) tab_Epetra_BlockMap.remove(&aid);
            break;
        case CT_Epetra_Import_ID:
            newid = tab.alias(tab_Epetra_Import.getConst<T>(aid));
            if (!keepold) tab_Epetra_Import.remove(&aid);
            break;
        case CT_Epetra_Time_ID:
            newid = tab.alias(tab_Epetra_Time.getConst<T>(aid));
            if (!keepold) tab_Epetra_Time.remove(&aid);
            break;
        case CT_Epetra_JadMatrix_ID:
            newid = tab.alias(tab_Epetra_JadMatrix.getConst<T>(aid));
            if (!keepold) tab_Epetra_JadMatrix.remove(&aid);
            break;
        case CT_Epetra_LinearProblem_ID:
            newid = tab.alias(tab_Epetra_LinearProblem.getConst<T>(aid));
            if (!keepold) tab_Epetra_LinearProblem.remove(&aid);
            break;
        case CT_Epetra_LAPACK_ID:
            newid = tab.alias(tab_Epetra_LAPACK.getConst<T>(aid));
            if (!keepold) tab_Epetra_LAPACK.remove(&aid);
            break;
        case CT_Teuchos_ParameterList_ID:
            newid = tab.alias(tab_Teuchos_ParameterList.getConst<T>(aid));
            if (!keepold) tab_Teuchos_ParameterList.remove(&aid);
            break;
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_BaseSolver_ID:
            newid = tab.alias(tab_Amesos_BaseSolver.getConst<T>(aid));
            if (!keepold) tab_Amesos_BaseSolver.remove(&aid);
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
        case CT_Epetra_FECrsMatrix_ID:
            newid = tab.alias(tab_Epetra_FECrsMatrix.getConst<T>(aid));
            if (!keepold) tab_Epetra_FECrsMatrix.remove(&aid);
            break;
        case CT_Epetra_IntSerialDenseVector_ID:
            newid = tab.alias(tab_Epetra_IntSerialDenseVector.getConst<T>(aid));
            if (!keepold) tab_Epetra_IntSerialDenseVector.remove(&aid);
            break;
        case CT_Epetra_SerialDenseMatrix_ID:
            newid = tab.alias(tab_Epetra_SerialDenseMatrix.getConst<T>(aid));
            if (!keepold) tab_Epetra_SerialDenseMatrix.remove(&aid);
            break;
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTest_ID:
            newid = tab.alias(tab_AztecOO_StatusTest.getConst<T>(aid));
            if (!keepold) tab_AztecOO_StatusTest.remove(&aid);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestCombo_ID:
            newid = tab.alias(tab_AztecOO_StatusTestCombo.getConst<T>(aid));
            if (!keepold) tab_AztecOO_StatusTestCombo.remove(&aid);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestMaxIters_ID:
            newid = tab.alias(tab_AztecOO_StatusTestMaxIters.getConst<T>(aid));
            if (!keepold) tab_AztecOO_StatusTestMaxIters.remove(&aid);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestResNorm_ID:
            newid = tab.alias(tab_AztecOO_StatusTestResNorm.getConst<T>(aid));
            if (!keepold) tab_AztecOO_StatusTestResNorm.remove(&aid);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
        case CT_Ifpack_Preconditioner_ID:
            newid = tab.alias(tab_Ifpack_Preconditioner.getConst<T>(aid));
            if (!keepold) tab_Ifpack_Preconditioner.remove(&aid);
            break;
#endif /* HAVE_CTRILINOS_IFPACK */
        default:
            newid.table = CT_Invalid_ID;
            newid.index = -1;
            newid.is_const = aid.is_const;
            throw std::string("invalid table id or non-polymorphic class");
            break;
    }

    return newid;
}

template <typename S1, typename S2>
S1 TableRepos::alias(S2 id, CTrilinos_Table_ID_t tab, bool keepold)
{
    CTrilinos_Universal_ID_t aid = abstractType<S2>(id);
    CTrilinos_Universal_ID_t newid;

    switch (tab) {
        case CT_Epetra_Distributor_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_Distributor, aid, keepold)
                                  : do_alias(tab_Epetra_Distributor, aid, keepold));
            break;
        case CT_Epetra_SerialComm_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_SerialComm, aid, keepold)
                                  : do_alias(tab_Epetra_SerialComm, aid, keepold));
            break;
        case CT_Epetra_BLAS_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_BLAS, aid, keepold)
                                  : do_alias(tab_Epetra_BLAS, aid, keepold));
            break;
        case CT_Epetra_Comm_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_Comm, aid, keepold)
                                  : do_alias(tab_Epetra_Comm, aid, keepold));
            break;
        case CT_Epetra_Operator_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_Operator, aid, keepold)
                                  : do_alias(tab_Epetra_Operator, aid, keepold));
            break;
        case CT_Epetra_MultiVector_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_MultiVector, aid, keepold)
                                  : do_alias(tab_Epetra_MultiVector, aid, keepold));
            break;
        case CT_Epetra_OffsetIndex_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_OffsetIndex, aid, keepold)
                                  : do_alias(tab_Epetra_OffsetIndex, aid, keepold));
            break;
        case CT_Epetra_Object_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_Object, aid, keepold)
                                  : do_alias(tab_Epetra_Object, aid, keepold));
            break;
        case CT_Epetra_RowMatrix_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_RowMatrix, aid, keepold)
                                  : do_alias(tab_Epetra_RowMatrix, aid, keepold));
            break;
        case CT_Epetra_CompObject_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_CompObject, aid, keepold)
                                  : do_alias(tab_Epetra_CompObject, aid, keepold));
            break;
        case CT_Epetra_Directory_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_Directory, aid, keepold)
                                  : do_alias(tab_Epetra_Directory, aid, keepold));
            break;
        case CT_Epetra_Flops_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_Flops, aid, keepold)
                                  : do_alias(tab_Epetra_Flops, aid, keepold));
            break;
        case CT_Epetra_SrcDistObject_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_SrcDistObject, aid, keepold)
                                  : do_alias(tab_Epetra_SrcDistObject, aid, keepold));
            break;
#ifdef HAVE_MPI
        case CT_Epetra_MpiComm_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_MpiComm, aid, keepold)
                                  : do_alias(tab_Epetra_MpiComm, aid, keepold));
            break;
#endif /* HAVE_MPI */
        case CT_Epetra_CrsMatrix_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_CrsMatrix, aid, keepold)
                                  : do_alias(tab_Epetra_CrsMatrix, aid, keepold));
            break;
        case CT_Epetra_CrsGraph_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_CrsGraph, aid, keepold)
                                  : do_alias(tab_Epetra_CrsGraph, aid, keepold));
            break;
        case CT_Epetra_DistObject_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_DistObject, aid, keepold)
                                  : do_alias(tab_Epetra_DistObject, aid, keepold));
            break;
        case CT_Epetra_Vector_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_Vector, aid, keepold)
                                  : do_alias(tab_Epetra_Vector, aid, keepold));
            break;
        case CT_Epetra_Export_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_Export, aid, keepold)
                                  : do_alias(tab_Epetra_Export, aid, keepold));
            break;
        case CT_Epetra_Map_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_Map, aid, keepold)
                                  : do_alias(tab_Epetra_Map, aid, keepold));
            break;
        case CT_Epetra_BlockMap_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_BlockMap, aid, keepold)
                                  : do_alias(tab_Epetra_BlockMap, aid, keepold));
            break;
        case CT_Epetra_Import_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_Import, aid, keepold)
                                  : do_alias(tab_Epetra_Import, aid, keepold));
            break;
        case CT_Epetra_Time_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_Time, aid, keepold)
                                  : do_alias(tab_Epetra_Time, aid, keepold));
            break;
        case CT_Epetra_JadMatrix_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_JadMatrix, aid, keepold)
                                  : do_alias(tab_Epetra_JadMatrix, aid, keepold));
            break;
        case CT_Epetra_LinearProblem_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_LinearProblem, aid, keepold)
                                  : do_alias(tab_Epetra_LinearProblem, aid, keepold));
            break;
        case CT_Epetra_LAPACK_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_LAPACK, aid, keepold)
                                  : do_alias(tab_Epetra_LAPACK, aid, keepold));
            break;
        case CT_Teuchos_ParameterList_ID:
            newid = (aid.is_const ? do_alias_const(tab_Teuchos_ParameterList, aid, keepold)
                                  : do_alias(tab_Teuchos_ParameterList, aid, keepold));
            break;
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_BaseSolver_ID:
            newid = (aid.is_const ? do_alias_const(tab_Amesos_BaseSolver, aid, keepold)
                                  : do_alias(tab_Amesos_BaseSolver, aid, keepold));
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
        case CT_Epetra_FECrsMatrix_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_FECrsMatrix, aid, keepold)
                                  : do_alias(tab_Epetra_FECrsMatrix, aid, keepold));
            break;
        case CT_Epetra_IntSerialDenseVector_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_IntSerialDenseVector, aid, keepold)
                                  : do_alias(tab_Epetra_IntSerialDenseVector, aid, keepold));
            break;
        case CT_Epetra_SerialDenseMatrix_ID:
            newid = (aid.is_const ? do_alias_const(tab_Epetra_SerialDenseMatrix, aid, keepold)
                                  : do_alias(tab_Epetra_SerialDenseMatrix, aid, keepold));
            break;
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTest_ID:
            newid = (aid.is_const ? do_alias_const(tab_AztecOO_StatusTest, aid, keepold)
                                  : do_alias(tab_AztecOO_StatusTest, aid, keepold));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestCombo_ID:
            newid = (aid.is_const ? do_alias_const(tab_AztecOO_StatusTestCombo, aid, keepold)
                                  : do_alias(tab_AztecOO_StatusTestCombo, aid, keepold));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestMaxIters_ID:
            newid = (aid.is_const ? do_alias_const(tab_AztecOO_StatusTestMaxIters, aid, keepold)
                                  : do_alias(tab_AztecOO_StatusTestMaxIters, aid, keepold));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestResNorm_ID:
            newid = (aid.is_const ? do_alias_const(tab_AztecOO_StatusTestResNorm, aid, keepold)
                                  : do_alias(tab_AztecOO_StatusTestResNorm, aid, keepold));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
        case CT_Ifpack_Preconditioner_ID:
            newid = (aid.is_const ? do_alias_const(tab_Ifpack_Preconditioner, aid, keepold)
                                  : do_alias(tab_Ifpack_Preconditioner, aid, keepold));
            break;
#endif /* HAVE_CTRILINOS_IFPACK */
        default:
            newid.table = CT_Invalid_ID;
            newid.index = -1;
            newid.is_const = aid.is_const;
            throw std::string("invalid table id or non-polymorphic class");
            break;
    }

    return concreteType<S1>(newid);
}

template <typename S>
void TableRepos::migrate(S *id, CTrilinos_Table_ID_t tab)
{
    S newid = alias<S,S>(*id, tab, false);
    *id = newid;
}

template <typename S>
void TableRepos::remove(S * id)
{
    CTrilinos_Universal_ID_t aid = abstractType<S>(*id);

    switch (aid.table) {
        case CT_Epetra_Distributor_ID:
            tab_Epetra_Distributor.remove(&aid);
            break;
        case CT_Epetra_SerialComm_ID:
            tab_Epetra_SerialComm.remove(&aid);
            break;
        case CT_Epetra_BLAS_ID:
            tab_Epetra_BLAS.remove(&aid);
            break;
        case CT_Epetra_Comm_ID:
            tab_Epetra_Comm.remove(&aid);
            break;
        case CT_Epetra_Operator_ID:
            tab_Epetra_Operator.remove(&aid);
            break;
        case CT_Epetra_MultiVector_ID:
            tab_Epetra_MultiVector.remove(&aid);
            break;
        case CT_Epetra_OffsetIndex_ID:
            tab_Epetra_OffsetIndex.remove(&aid);
            break;
        case CT_Epetra_Object_ID:
            tab_Epetra_Object.remove(&aid);
            break;
        case CT_Epetra_RowMatrix_ID:
            tab_Epetra_RowMatrix.remove(&aid);
            break;
        case CT_Epetra_CompObject_ID:
            tab_Epetra_CompObject.remove(&aid);
            break;
        case CT_Epetra_Directory_ID:
            tab_Epetra_Directory.remove(&aid);
            break;
        case CT_Epetra_Flops_ID:
            tab_Epetra_Flops.remove(&aid);
            break;
        case CT_Epetra_SrcDistObject_ID:
            tab_Epetra_SrcDistObject.remove(&aid);
            break;
#ifdef HAVE_MPI
        case CT_Epetra_MpiComm_ID:
            tab_Epetra_MpiComm.remove(&aid);
            break;
#endif /* HAVE_MPI */
        case CT_Epetra_CrsMatrix_ID:
            tab_Epetra_CrsMatrix.remove(&aid);
            break;
        case CT_Epetra_CrsGraph_ID:
            tab_Epetra_CrsGraph.remove(&aid);
            break;
        case CT_Epetra_DistObject_ID:
            tab_Epetra_DistObject.remove(&aid);
            break;
        case CT_Epetra_Vector_ID:
            tab_Epetra_Vector.remove(&aid);
            break;
        case CT_Epetra_Export_ID:
            tab_Epetra_Export.remove(&aid);
            break;
        case CT_Epetra_Map_ID:
            tab_Epetra_Map.remove(&aid);
            break;
        case CT_Epetra_BlockMap_ID:
            tab_Epetra_BlockMap.remove(&aid);
            break;
        case CT_Epetra_Import_ID:
            tab_Epetra_Import.remove(&aid);
            break;
        case CT_Epetra_Time_ID:
            tab_Epetra_Time.remove(&aid);
            break;
        case CT_Epetra_JadMatrix_ID:
            tab_Epetra_JadMatrix.remove(&aid);
            break;
        case CT_Epetra_LinearProblem_ID:
            tab_Epetra_LinearProblem.remove(&aid);
            break;
        case CT_Epetra_LAPACK_ID:
            tab_Epetra_LAPACK.remove(&aid);
            break;
        case CT_Teuchos_CommandLineProcessor_ID:
            tab_Teuchos_CommandLineProcessor.remove(&aid);
            break;
        case CT_Teuchos_ParameterList_ID:
            tab_Teuchos_ParameterList.remove(&aid);
            break;
        case CT_Teuchos_ParameterEntry_ID:
            tab_Teuchos_ParameterEntry.remove(&aid);
            break;
        case CT_Teuchos_any_ID:
            tab_Teuchos_any.remove(&aid);
            break;
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_BaseSolver_ID:
            tab_Amesos_BaseSolver.remove(&aid);
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_ID:
            tab_Amesos.remove(&aid);
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
        case CT_Epetra_FECrsMatrix_ID:
            tab_Epetra_FECrsMatrix.remove(&aid);
            break;
        case CT_Epetra_IntSerialDenseVector_ID:
            tab_Epetra_IntSerialDenseVector.remove(&aid);
            break;
        case CT_Epetra_SerialDenseMatrix_ID:
            tab_Epetra_SerialDenseMatrix.remove(&aid);
            break;
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_ID:
            tab_AztecOO.remove(&aid);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTest_ID:
            tab_AztecOO_StatusTest.remove(&aid);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestCombo_ID:
            tab_AztecOO_StatusTestCombo.remove(&aid);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestMaxIters_ID:
            tab_AztecOO_StatusTestMaxIters.remove(&aid);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestResNorm_ID:
            tab_AztecOO_StatusTestResNorm.remove(&aid);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
        case CT_Ifpack_ID:
            tab_Ifpack.remove(&aid);
            break;
#endif /* HAVE_CTRILINOS_IFPACK */
#ifdef HAVE_CTRILINOS_IFPACK
        case CT_Ifpack_Preconditioner_ID:
            tab_Ifpack_Preconditioner.remove(&aid);
            break;
#endif /* HAVE_CTRILINOS_IFPACK */
        default:
            throw std::string("invalid table id");
            break;
    }

    *id = concreteType<S>(aid);
}

void TableRepos::purge()
{
    tab_Epetra_Distributor.purge();
    tab_Epetra_SerialComm.purge();
    tab_Epetra_BLAS.purge();
    tab_Epetra_Comm.purge();
    tab_Epetra_Operator.purge();
    tab_Epetra_MultiVector.purge();
    tab_Epetra_OffsetIndex.purge();
    tab_Epetra_Object.purge();
    tab_Epetra_RowMatrix.purge();
    tab_Epetra_CompObject.purge();
    tab_Epetra_Directory.purge();
    tab_Epetra_Flops.purge();
    tab_Epetra_SrcDistObject.purge();
#ifdef HAVE_MPI
    tab_Epetra_MpiComm.purge();
#endif /* HAVE_MPI */
    tab_Epetra_CrsMatrix.purge();
    tab_Epetra_CrsGraph.purge();
    tab_Epetra_DistObject.purge();
    tab_Epetra_Vector.purge();
    tab_Epetra_Export.purge();
    tab_Epetra_Map.purge();
    tab_Epetra_BlockMap.purge();
    tab_Epetra_Import.purge();
    tab_Epetra_Time.purge();
    tab_Epetra_JadMatrix.purge();
    tab_Epetra_LinearProblem.purge();
    tab_Epetra_LAPACK.purge();
    tab_Teuchos_CommandLineProcessor.purge();
    tab_Teuchos_ParameterList.purge();
    tab_Teuchos_ParameterEntry.purge();
    tab_Teuchos_any.purge();
#ifdef HAVE_CTRILINOS_AMESOS
    tab_Amesos_BaseSolver.purge();
#endif /* HAVE_CTRILINOS_AMESOS */
#ifdef HAVE_CTRILINOS_AMESOS
    tab_Amesos.purge();
#endif /* HAVE_CTRILINOS_AMESOS */
    tab_Epetra_FECrsMatrix.purge();
    tab_Epetra_IntSerialDenseVector.purge();
    tab_Epetra_SerialDenseMatrix.purge();
#ifdef HAVE_CTRILINOS_AZTECOO
    tab_AztecOO.purge();
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    tab_AztecOO_StatusTest.purge();
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    tab_AztecOO_StatusTestCombo.purge();
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    tab_AztecOO_StatusTestMaxIters.purge();
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    tab_AztecOO_StatusTestResNorm.purge();
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
    tab_Ifpack.purge();
#endif /* HAVE_CTRILINOS_IFPACK */
#ifdef HAVE_CTRILINOS_IFPACK
    tab_Ifpack_Preconditioner.purge();
#endif /* HAVE_CTRILINOS_IFPACK */
}


} // namespace CTrilinos


#endif

