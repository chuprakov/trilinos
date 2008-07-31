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


// #include <Isorropia_Zoltan_Repartition.hpp>
#ifdef HAVE_ISORROPIA_ZOLTAN

#endif
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Isorropia_EpetraZoltanLib.hpp>

#ifndef HAVE_MPI
#error "Isorropia_Zoltan requires MPI."
#endif


#ifdef HAVE_EPETRA
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#endif
#endif


#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <ctype.h>

/* TODO: clean up the code */

namespace Isorropia {

#ifdef HAVE_EPETRA

namespace Epetra {

ZoltanLibClass::ZoltanLibClass(Teuchos::RefCountPtr<const Epetra_CrsGraph> input_graph):
  Library(input_graph)
{
  setInputType("HYPERGRAPH");
}

ZoltanLibClass::ZoltanLibClass(Teuchos::RefCountPtr<const Epetra_CrsGraph> input_graph,
			  Teuchos::RefCountPtr<CostDescriber> costs):
  Library(input_graph, costs)
{
  setInputType("HYPERGRAPH");
}

ZoltanLibClass::ZoltanLibClass(Teuchos::RefCountPtr<const Epetra_RowMatrix> input_matrix):
  Library(input_matrix)
{
  setInputType("HYPERGRAPH");
}

ZoltanLibClass::ZoltanLibClass(Teuchos::RefCountPtr<const Epetra_RowMatrix> input_matrix,
			  Teuchos::RefCountPtr<CostDescriber> costs):
  Library(input_matrix, costs)
{
  setInputType("HYPERGRAPH");
}

int ZoltanLibClass::precompute()
{
  int err = 0, numChanged = 0;
  std::string str1("Isorropia::ZoltanLibClass::precompute ");
  std::string str2;
  MPI_Comm mpicomm;

  bool isHypergraph = true;
  std::string lb_method_str("LB_METHOD");
  if (zoltanParamList_.isParameter(lb_method_str)){
    std::string lb_meth = zoltanParamList_.get(lb_method_str, "HYPERGRAPH");
    if (lb_meth == "GRAPH"){
      isHypergraph = false;
      inputType_ = "GRAPH";
    }
  }


  Library::precompute();
//   if (!isHypergraph){
//     bool square = false;
//     bool symmetric = false;
//     if (input_graph_.get() != 0){
//       if (input_graph_->NumGlobalRows() == input_graph_->NumGlobalCols()){
// 	square = true;
// 	// TODO - is there a quick way to figure out if the graph is
// 	// symmetric?  I can't see a way to do it.  For now we let
// 	// Zoltan figure this out.
// 	symmetric = true;
//       }
//     }
//     else{
//       if (input_matrix_->NumGlobalRows() == input_matrix_->NumGlobalCols()){
// 	square = true;
// 	// TODO - is there a quick way to figure out if the matrix is
// 	// symmetric?  I can't see a way to do it.  For now we let
// 	// Zoltan figure this out.
// 	symmetric = true;
//       }
//     }
//     if (!square){
//       str2 = "LB_METHOD=GRAPH, matrix or graph must be square";
//       throw Isorropia::Exception(str1+str2);
//     }
//     if (!symmetric){
//       str2 = "LB_METHOD=GRAPH, matrix or graph must be symmetric";
//       throw Isorropia::Exception(str1+str2);
//     }
//   }

  computeCost();

//   preCheckPartition();


  if (input_matrix_.get() == 0) {
    queryObject_ =  Teuchos::rcp(new ZoltanLib::QueryObject(input_graph_, costs_, isHypergraph));
    const  Epetra_Comm &ecomm = input_graph_->RowMap().Comm();
#ifdef HAVE_MPI
    const Epetra_MpiComm &empicomm = dynamic_cast<const Epetra_MpiComm &>(ecomm);
    mpicomm = empicomm.Comm();

#else /* HAVE_MPI */
  return (-1);
#endif /* HAVE_MPI */
  }
  else {
    queryObject_ =  Teuchos::rcp(new ZoltanLib::QueryObject(input_matrix_, costs_, isHypergraph));
    const Epetra_Comm &ecomm = input_matrix_->RowMatrixRowMap().Comm();
#ifdef HAVE_MPI
    const Epetra_MpiComm &empicomm = dynamic_cast<const Epetra_MpiComm &>(ecomm);

    mpicomm = empicomm.Comm();
#else /* HAVE_MPI */
  return (-1);
#endif /* HAVE_MPI */
  }

  float version;
  int argcTmp=0;
  char *argvTmp[1];

  // create a Zoltan problem

  argvTmp[0] = NULL;
  Zoltan_Initialize(argcTmp, argvTmp, &version);

  zz_ = new Zoltan(mpicomm);

  if (zz_ == NULL){
    throw Isorropia::Exception("Error creating Zoltan object");
    return (-1);
  }

  // set problem parameters

  std::string dbg_level_str("DEBUG_LEVEL");
  if (!zoltanParamList_.isParameter(dbg_level_str)) {
    zoltanParamList_.set(dbg_level_str, "0");
  }

  std::string lb_meth = zoltanParamList_.get(lb_method_str, "HYPERGRAPH");

  if (!zoltanParamList_.isParameter(lb_method_str)) {
    zoltanParamList_.set(lb_method_str, lb_meth);  // set to HYPERGRAPH
  }

  // If no process set vertex weights, let Zoltan default to
  // unit weights for vertices

  if (queryObject_->haveVertexWeights()) {
    if (!zoltanParamList_.isParameter("OBJ_WEIGHT_DIM")) {
      zoltanParamList_.set("OBJ_WEIGHT_DIM", "1");
    }
  }

  // If no process set graph or hypergraph edge weights,
  // let Zoltan default to unit weights for edges

  if (queryObject_->haveGraphEdgeWeights() ||
      queryObject_->haveHypergraphEdgeWeights()) {
    if (!zoltanParamList_.isParameter("EDGE_WEIGHT_DIM")) {
      zoltanParamList_.set("EDGE_WEIGHT_DIM", "1");
    }
  }


  Teuchos::ParameterList::ConstIterator
    iter = zoltanParamList_.begin(),
    iter_end = zoltanParamList_.end();

  for(; iter != iter_end; ++iter) {
    const std::string& name = iter->first;
    const std::string& value = Teuchos::getValue<std::string>(iter->second);
    zz_->Set_Param(name, value);
  }

  // Set the query functions

  zz_->Set_Num_Obj_Fn(ZoltanLib::QueryObject::Number_Objects, (void *)queryObject_.get());
  zz_->Set_Obj_List_Fn(ZoltanLib::QueryObject::Object_List, (void *)queryObject_.get());

  int ierr;
  num_obj_ = ZoltanLib::QueryObject::Number_Objects((void *)queryObject_.get(), &ierr);

  if (lb_meth == "HYPERGRAPH"){
    zz_->Set_HG_Size_CS_Fn(ZoltanLib::QueryObject::HG_Size_CS, (void *)queryObject_.get());
    zz_->Set_HG_CS_Fn(ZoltanLib::QueryObject::HG_CS, (void *)queryObject_.get());
    zz_->Set_HG_Size_Edge_Wts_Fn(ZoltanLib::QueryObject::HG_Size_Edge_Weights,
				 (void *)queryObject_.get());
    zz_->Set_HG_Edge_Wts_Fn(ZoltanLib::QueryObject::HG_Edge_Weights, (void *)queryObject_.get());
  }
  else{
    zz_->Set_Num_Edges_Multi_Fn(ZoltanLib::QueryObject::Number_Edges_Multi, (void *)queryObject_.get());
    zz_->Set_Edge_List_Multi_Fn(ZoltanLib::QueryObject::Edge_List_Multi, (void *)queryObject_.get());
  }

}


void ZoltanLibClass::computeCost()
{
  int err = 0, numChanged = 0;
  std::string str1("Isorropia::ZoltanLibClass::computeCost ");
  std::string str2;


  const Epetra_Comm &comm = input_map_->Comm();
  int localProc = comm.MyPID();
  int nprocs = comm.NumProc();

    // If vertex/edge costs have been set, do a global operation to find
    // out how many weights were given.  Some processes may provide no
    // weights - they need to be informed that weights are being provided
    // by the application.  Do some sanity checks.

    err = 0;
    int gerr = 0;
    int base = input_map_->IndexBase();

    int numMyVWeights = 0;
    int numMyGWeights = 0;
    int numMyHGWeights = 0;
    int globalNumCols = 0;
    int myNZ = 0;
    int globalNZ = 0;
    int mySelfEdges = 0;
    int globalSelfEdges = 0;

    int myRows = input_map_->NumMyElements();
    int globalNumRows = input_map_->NumGlobalElements();

    if (input_graph_.get() == 0) {
      myNZ = input_matrix_->NumMyNonzeros();
      mySelfEdges = input_matrix_->NumMyDiagonals();
      globalNZ = input_matrix_->NumGlobalNonzeros();
      globalSelfEdges = input_matrix_->NumGlobalDiagonals();
      globalNumCols = input_matrix_->NumGlobalCols();
    }
    else{
      myNZ = input_graph_->NumMyNonzeros();
      mySelfEdges = input_graph_->NumMyDiagonals();
      globalNZ = input_graph_->NumGlobalNonzeros();
      globalSelfEdges = input_graph_->NumGlobalDiagonals();
      globalNumCols = input_graph_->NumGlobalCols();
    }

    if (costs_.get() != 0) {

      numMyVWeights = costs_->getNumVertices();

      if (costs_->haveGraphEdgeWeights()){
        for (int i=0; i<numMyVWeights; i++){
          int gid = input_map_->GID(i);
          if (gid >= base){
            numMyGWeights += costs_->getNumGraphEdges(gid);
          }
        }
      }
      numMyHGWeights = costs_->getNumHypergraphEdgeWeights();

      if ((numMyVWeights > 0) && (numMyVWeights != myRows)){
        str2 = "Number of my vertex weights != number of my rows";
        err = 1;
      }
      else if ((numMyGWeights > 0) && (numMyGWeights != (myNZ - mySelfEdges))){
        str2 = "Number of my graph edge weights != number of my nonzeros";
        err = 1;
      }
    }
    else{
      costs_ = Teuchos::rcp(new CostDescriber());
    }

    comm.SumAll(&err, &gerr ,1);

    if (gerr > 0){
      throw Isorropia::Exception(str1+str2);
    }

    int lval[4], gval[4];
    lval[0] = numMyVWeights;
    lval[1] = numMyGWeights;
    lval[2] = numMyHGWeights;

    comm.SumAll(lval, gval, 3);

    int numVWeights = gval[0];
    int numGWeights = gval[1];
    int numHGWeights = gval[2];

    if ((numVWeights > 0) && (numVWeights != globalNumRows)){
      str2 = "Number of vertex weights supplied by application != number of rows";
      throw Isorropia::Exception(str1+str2);
    }
    if ((numGWeights > 0) && (numGWeights != (globalNZ - globalSelfEdges))){
      str2 = "Number of graph edge weights supplied by application != number of edges";
      throw Isorropia::Exception(str1+str2);
    }
    if ((numHGWeights > 0) && (numHGWeights < globalNumCols)){
      str2 = "Number of hyperedge weights supplied by application < number of columns";
      throw Isorropia::Exception(str1+str2);
    }

    costs_->setNumGlobalVertexWeights(numVWeights);
    costs_->setNumGlobalGraphEdgeWeights(numGWeights);
    costs_->setNumGlobalHypergraphEdgeWeights(numHGWeights);

}


void ZoltanLibClass::preCheckPartition()
{

  int err = 0, numChanged = 0;
  std::string str1("Isorropia::ZoltanLibClass::precheckPartition ");
  std::string str2;

  const Epetra_Comm &comm = input_map_->Comm();
  int localProc = comm.MyPID();
  int nprocs = comm.NumProc();

  // Checks for Zoltan parameters NUM_GLOBAL_PARTITIONS and NUM_LOCAL_PARTITIONS.

  std::string gparts_str("NUM_GLOBAL_PARTITIONS");
  std::string lparts_str("NUM_LOCAL_PARTITIONS");
  std::string gparts("0");
  std::string lparts("0");

  if (zoltanParamList_.isParameter(gparts_str)){
    gparts = zoltanParamList_.get<std::string>(gparts_str);
  }
  if (zoltanParamList_.isParameter(lparts_str)){
    lparts = zoltanParamList_.get<std::string>(lparts_str);
  }

  int myGparts = atoi(gparts.c_str());
  int myLparts = atoi(lparts.c_str());
  int maxGparts, maxLparts, sumLparts;
  int fixGparts = -1;
  int fixLparts = -1;
  int numrows = input_map_->NumGlobalElements();

  int myParts[] = {myGparts, myLparts};
  int maxParts[2];

  comm.MaxAll(myParts, maxParts, 2);

  maxGparts = maxParts[0];
  maxLparts = maxParts[1];

  // Fix problem if the number of rows is less than the number
  // of processes.  We need to set NUM_GLOBAL_PARTITIONS to
  // the number of rows, unless it was already set to something
  // greater than 0 but less than the number of rows.

  if (numrows < nprocs){
    if ((maxGparts == 0) || (maxGparts > numrows)){
      maxGparts = numrows;
    }
  }

  if (maxLparts > 0)
    comm.SumAll(&myLparts, &sumLparts, 1);
  else
    sumLparts = 0;

  if ((maxGparts > 0) || (maxLparts > 0)){
    // One or more processes set NGP or NLP so we need to check
    // them for validity.

  if (maxGparts != myGparts){
    fixGparts = maxGparts;    // all procs should have same NGP
  }

  if (maxGparts > 0){
    if (maxGparts > numrows){
      // This is an error because we can't split rows among partitions
      str2 = "NUM_GLOBAL_PARTITIONS exceeds number of rows (objects to be partitioned)";
      throw Isorropia::Exception(str1+str2);
    }

    if ((sumLparts > 0) && (sumLparts != maxGparts)){
      // This is an error because local partitions must sum to number of global
      str2 = "NUM_GLOBAL_PARTITIONS not equal to sum of NUM_LOCAL_PARTITIONS";
      throw Isorropia::Exception(str1+str2);
    }

    if ((sumLparts == 0) && (maxGparts < nprocs)){
      // Set NUM_LOCAL_PARTITIONS to 1 or 0, because Zoltan will divide
      // a partition across 2 or more processes when the number of
      // partitions is less than the number of processes.  This doesn't
      // work for Epetra matrices, where rows are not owned by more than
      // one process.

      fixLparts = (localProc < maxGparts) ? 1 : 0;
    }
  }
  else if (maxLparts > 0){

    // Set NUM_GLOBAL_PARTITIONS to sum of local partitions.  It's possible
    // that Zoltan does this already, but just to be safe...

    fixGparts = sumLparts;
  }
  if (fixGparts > 0){
    std::ostringstream os;
    os << fixGparts;
    std::string s = os.str();
    zoltanParamList_.set(gparts_str, s);
  }
  if (fixLparts > 0){
    std::ostringstream os;
    os << fixLparts;
    std::string s = os.str();
    zoltanParamList_.set(lparts_str, s);
  }
  }
}

int ZoltanLibClass::
repartition(Teuchos::ParameterList& zoltanParamList,
	    std::vector<int>& myNewElements,
	    std::map<int,int>& exports,
	    std::map<int,int>& imports)
{
  zoltanParamList_ = zoltanParamList;

  precompute();
  preCheckPartition();

  std::string lb_approach_str("LB_APPROACH");
  if (!zoltanParamList_.isParameter(lb_approach_str)) {
    zoltanParamList_.set(lb_approach_str, "PARTITION");
  }

  //Generate Load Balance
  int changes, num_gid_entries, num_lid_entries, num_import, num_export;
  ZOLTAN_ID_PTR import_global_ids, import_local_ids;
  ZOLTAN_ID_PTR export_global_ids, export_local_ids;
  int * import_procs, * export_procs;
  int *import_to_part, *export_to_part;

  int err = zz_->LB_Partition(changes, num_gid_entries, num_lid_entries,
   num_import, import_global_ids, import_local_ids, import_procs, import_to_part,
   num_export, export_global_ids, export_local_ids, export_procs, export_to_part );

  if (err != ZOLTAN_OK){
    throw Isorropia::Exception("Error computing partitioning with Zoltan");
    return -1;
  }

  //Generate New Element List
  int numMyElements = queryObject_->RowMap().NumMyElements();
  std::vector<int> elementList( numMyElements );
  queryObject_->RowMap().MyGlobalElements( &elementList[0] );

  int newNumMyElements = numMyElements - num_export + num_import;
  myNewElements.resize( newNumMyElements );

  for( int i = 0; i < num_export; ++i ) {
    exports[export_global_ids[i]] = export_procs[i];
  }

  for( int i = 0; i < num_import; ++i ) {
    imports[import_global_ids[i]] = import_procs[i];
  }

  //Add unmoved indices to new list
  int loc = 0;
  for( int i = 0; i < numMyElements; ++i ) {
    if( !exports.count( elementList[i] ) ) {
      myNewElements[loc++] = elementList[i];
    }
  }

  //Add imports to end of list
  for( int i = 0; i < num_import; ++i ) {
    myNewElements[loc+i] = import_global_ids[i];
  }

  //Free Zoltan Data
  zz_->LB_Free_Part(&import_global_ids, &import_local_ids,
                     &import_procs, &import_to_part);
  zz_->LB_Free_Part(&export_global_ids, &export_local_ids,
                     &export_procs, &export_to_part);

  postcompute();
}

int ZoltanLibClass::
color(Teuchos::ParameterList& zoltanParamList,
      std::vector<int>& myNewElements)
{
  zoltanParamList_ = zoltanParamList;

  precompute();

  //Generate Load Balance
  int changes, num_gid_entries, num_lid_entries;
  ZOLTAN_ID_PTR import_global_ids, import_local_ids;
  int *colors = new int(num_obj_);

  int err = zz_->Color(num_gid_entries, num_lid_entries, num_obj_, import_global_ids, import_local_ids, colors);

  if (err != ZOLTAN_OK){
    throw Isorropia::Exception("Error computing partitioning with Zoltan");
    return -1;
  }

  for (int i = 0 ; i < num_obj_ ; ++i) {
    myNewElements[i] = colors[i];
  }

  delete[] colors;
}

int ZoltanLibClass::
order(Teuchos::ParameterList& zoltanParamList,
      std::vector<int>& myNewElements)
{
  zoltanParamList_ = zoltanParamList;

  precompute();

  //Generate Load Balance
  int changes, num_gid_entries, num_lid_entries;
  ZOLTAN_ID_PTR import_global_ids, import_local_ids;
  int *rank = new int(num_obj_);

  int err = zz_->Order(num_gid_entries, num_lid_entries, num_obj_, import_global_ids, import_local_ids, rank, NULL);

  if (err != ZOLTAN_OK){
    throw Isorropia::Exception("Error computing partitioning with Zoltan");
    return -1;
  }

  for (int i = 0 ; i < num_obj_ ; ++i) {
    myNewElements[i] = rank[i];
  }

  delete[] rank;

}



int ZoltanLibClass::postcompute()
{
  if (zz_)
    delete zz_;
  zz_ = NULL;

  return (0);
}




} // namespace EPETRA

#endif //HAVE_EPETRA

}//namespace Isorropia

