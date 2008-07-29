//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

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
Questions? Contact Alan Williams (william@sandia.gov)
                or Erik Boman    (egboman@sandia.gov)

************************************************************************
*/
//@HEADER

#include <Isorropia_EpetraOperator.hpp>
#ifdef HAVE_ISORROPIA_ZOLTAN
#include <Isorropia_Zoltan_Repartition.hpp>
#endif
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>

#endif

#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <ctype.h>

namespace Isorropia {

#ifdef HAVE_EPETRA

namespace Epetra {

Operator::
Operator(Teuchos::RefCountPtr<const Epetra_CrsGraph> input_graph,
	 const Teuchos::ParameterList& paramlist)
  : input_map_(),
    input_graph_(input_graph),
    input_matrix_(),
    costs_(),
    paramlist_(),
    operation_already_computed_(false)
{
  input_map_ = Teuchos::rcp(&(input_graph->RowMap()), false);
  paramlist_ = paramlist;
}

Operator::
Operator(Teuchos::RefCountPtr<const Epetra_CrsGraph> input_graph,
	 Teuchos::RefCountPtr<CostDescriber> costs,
	 const Teuchos::ParameterList& paramlist)
  : input_map_(),
    input_graph_(input_graph),
    input_matrix_(),
    costs_(costs),
    paramlist_(),
    operation_already_computed_(false)
{

  input_map_ = Teuchos::rcp(&(input_graph->RowMap()), false);
  paramlist_ = paramlist;
}

Operator::
Operator(Teuchos::RefCountPtr<const Epetra_RowMatrix> input_matrix,
	 const Teuchos::ParameterList& paramlist)
  : input_map_(),
    input_graph_(),
    input_matrix_(input_matrix),
    costs_(),
    paramlist_(),
    operation_already_computed_(false)
{
  input_map_ = Teuchos::rcp(&(input_matrix->RowMatrixRowMap()),false);
  paramlist_ = paramlist;
}

Operator::
Operator(Teuchos::RefCountPtr<const Epetra_RowMatrix> input_matrix,
	 Teuchos::RefCountPtr<CostDescriber> costs,
	 const Teuchos::ParameterList& paramlist)
  : input_map_(),
    input_graph_(),
    input_matrix_(input_matrix),
    costs_(costs),
    paramlist_(),
    operation_already_computed_(false)
{
  input_map_ = Teuchos::rcp(&(input_matrix->RowMatrixRowMap()),false);
  paramlist_ = paramlist;
}

Operator::~Operator()
{
}

void Operator::setParameters(const Teuchos::ParameterList& paramlist)
{
  paramlist_ = paramlist;
}

bool Operator::alreadyComputed() const
{
  return operation_already_computed_;
}

int Operator::getNewPropertyOfElem(int myElem) const
{
  std::map<int,int>::const_iterator iter = exports_.find(myElem);
  if (iter != exports_.end()) {
    return(iter->second);
  }

  return( input_graph_->RowMap().Comm().MyPID() );
}

int Operator::getNbrElemsWithProperty(int partition) const
{
  int myPart = input_map_->Comm().MyPID();
  if (partition != myPart) {
    throw Isorropia::Exception("Partitioner::numElemsInPartition not implemented for non-local partitions.");
  }

  return(myNewElements_.size());
}

void
Operator::getElemsWithProperty(int partition, int* elementList, int len) const
{
  int myPart = input_map_->Comm().MyPID();
  if (partition != myPart) {
    throw Isorropia::Exception("error in Epetra_Map::MyGlobalElements");
  }

  unsigned length = len;
  if (myNewElements_.size() < length) length = myNewElements_.size();

  for(unsigned i=0; i<length; ++i) {
    elementList[i] = myNewElements_[i];
  }
}

int
Operator::getNumberOfProperties() const {
  return (numberOfProperties_);
}


void Operator::stringToUpper(std::string &s, int &changed)
{
  std::string::iterator siter;
  changed = 0;

  for (siter = s.begin(); siter != s.end() ; siter++)
  {
    if (islower(*siter)){
      *siter = toupper(*siter);
      changed++;
    }
  }
}

void Operator::paramsToUpper(Teuchos::ParameterList &plist, int &changed)
{
  changed = 0;

  // get a list of all parameter names in the list

  std::vector<std::string> paramNames ;
  Teuchos::ParameterList::ConstIterator pIter;

  pIter = plist.begin();

  while (1){
    //////////////////////////////////////////////////////////////////////
    // Compiler considered this while statement an error
    // while ( pIter = plist.begin() ; pIter != plist.end() ; pIter++ ){
    // }
    //////////////////////////////////////////////////////////////////////
    if (pIter == plist.end()) break;
    const std::string & nm = plist.name(pIter);
    paramNames.push_back(nm);
    pIter++;
  }

  // Change parameter names and values to upper case
  
  for (int i=0; i < paramNames.size(); i++){

    std::string origName(paramNames[i]);
    int paramNameChanged = 0;
    stringToUpper(paramNames[i], paramNameChanged);

    if (plist.isSublist(origName)){
      Teuchos::ParameterList &sublist = plist.sublist(origName);

      int sublistChanged=0;
      paramsToUpper(sublist, sublistChanged);

      if (paramNameChanged){

        // this didn't work, so I need to remove the old sublist
        // and create a new one
        //
        //sublist.setName(paramNames[i]);

        Teuchos::ParameterList newlist(sublist);
        plist.remove(origName);
        plist.set(paramNames[i], newlist);
      }
    }
    else if (plist.isParameter(origName)){

      std::string paramVal(plist.get<std::string>(origName));

      int paramValChanged=0;
      stringToUpper(paramVal, paramValChanged);

      if (paramNameChanged || paramValChanged){
        if (paramNameChanged){
          plist.remove(origName);
        }
        plist.set(paramNames[i], paramVal);
        changed++;
      }
    }
  } // next parameter or sublist
}

} // namespace EPETRA

#endif //HAVE_EPETRA

}//namespace Isorropia

