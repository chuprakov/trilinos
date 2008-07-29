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

#include <Isorropia_EpetraPartitioner.hpp>
#ifdef HAVE_ISORROPIA_ZOLTAN
#include <Isorropia_Zoltan_Repartition.hpp>
#include <Isorropia_EpetraZoltanLib.hpp>
#endif
#include <Isorropia_EpetraInternalPartitioner.hpp>
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


  /** Constructor that accepts an Epetra_CrsGraph object, called by
        API function create_partitioner().

     \param input_graph Matrix-graph object for which a new partitioning
        is to be computed. A Teuchos::RefCountPtr is used here because a
        reference to the input object may be held by this object after
        this constructor completes and returns.

     \param paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. No reference to this input
        object is held after this constructor completes.<br>
  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the balancing. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  Refer to the Zoltan users guide for specific parameters that Zoltan
  recognizes. A couple of important ones are "LB_METHOD" (valid values
  include "GRAPH", "HYPERGRAPH"), "DEBUG_LEVEL" (valid values are
  0 to 10, default is 1), etc.

     \param compute_partitioning_now Optional argument defaults to true.
        If true, the method compute_partitioning() will be called before
        this constructor returns.
  */
Partitioner::Partitioner(Teuchos::RefCountPtr<const Epetra_CrsGraph> input_graph,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_graph, paramlist)
{
  if (compute_partitioning_now)
    compute_partitioning(true);
}

  /** Constructor that accepts an Epetra_CrsGraph object and a CostDescriber, called by
        API function create_partitioner().

     \param input_graph Matrix-graph object for which a new partitioning
        is to be computed. A Teuchos::RefCountPtr is used here because a
        reference to the input object may be held by this object after
        this constructor completes and returns.

     \param costs CostDescriber object which allows for user-specified
       weights of varying types to be provided to the partitioner.

     \param paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. No reference to this input
        object is held after this constructor completes.<br>
  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the balancing. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  Refer to the Zoltan users guide for specific parameters that Zoltan
  recognizes. A couple of important ones are "LB_METHOD" (valid values
  include "GRAPH", "HYPERGRAPH"), "DEBUG_LEVEL" (valid values are
  0 to 10, default is 1), etc.

     \param compute_partitioning_now Optional argument defaults to true.
        If true, the method compute_partitioning() will be called before
        this constructor returns.
  */
Partitioner::Partitioner(Teuchos::RefCountPtr<const Epetra_CrsGraph> input_graph,
			 Teuchos::RefCountPtr<CostDescriber> costs,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_graph, costs, paramlist)
{
  if (compute_partitioning_now)
    compute_partitioning(true);
}


  /**
     Constructor that accepts an Epetra_RowMatrix object, called by
       API function create_partitioner().

     \param input_matrix Matrix object for which a new partitioning is
        to be computed. A Teuchos::RefCountPtr is used here because a
        reference to the input object may be held by this object after
        this constructor completes and returns.

     \param paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. No reference to this input
        object is held after this constructor completes.<br>
  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the balancing. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  Refer to the Zoltan users guide for specific parameters that Zoltan
  recognizes. A couple of important ones are "LB_METHOD" (valid values
  include "GRAPH", "HYPERGRAPH"), "DEBUG_LEVEL" (valid values are
  0 to 10, default is 1), etc.

     \param compute_partitioning_now Optional argument defaults to true.
        If true, the method compute_partitioning() will be called before
        this constructor returns.
  */
Partitioner::Partitioner(Teuchos::RefCountPtr<const Epetra_RowMatrix> input_matrix,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_matrix, paramlist)
{
  if (compute_partitioning_now)
    compute_partitioning(true);
}


  /**
     Constructor that accepts an Epetra_RowMatrix object and a
     CostDescriber, called by API function create_partitioner(). 

     \param input_matrix Matrix object for which a new partitioning is
        to be computed. A Teuchos::RefCountPtr is used here because a
        reference to the input object may be held by this object after
        this constructor completes and returns.

     \param costs CostDescriber object which allows for user-specified
       weights of varying types to be provided to the partitioner.

     \param paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. No reference to this input
        object is held after this constructor completes.<br>
  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the balancing. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  Refer to the Zoltan users guide for specific parameters that Zoltan
  recognizes. A couple of important ones are "LB_METHOD" (valid values
  include "GRAPH", "HYPERGRAPH"), "DEBUG_LEVEL" (valid values are
  0 to 10, default is 1), etc.

     \param compute_partitioning_now Optional argument defaults to true.
        If true, the method compute_partitioning() will be called before
        this constructor returns.
  */
Partitioner::Partitioner(Teuchos::RefCountPtr<const Epetra_RowMatrix> input_matrix,
			 Teuchos::RefCountPtr<CostDescriber> costs,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_matrix, costs, paramlist)
{
  if (compute_partitioning_now)
    compute_partitioning(true);
}


  /** Destructor */
Partitioner::~Partitioner(){}

  /** setParameters() is an internal Partitioner method which handles
      the parameters from a Teuchos::ParameterList object. 

      The input
      ParameterList object is copied into an internal ParameterList
      attribute, and no reference to the input object is held after
      this function returns. (Thus, the input paramlist object may be
      altered or destroyed as soon as this method returns.)<br>
  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the balancing. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  Refer to the Zoltan users guide for specific parameters that Zoltan
  recognizes. A couple of important ones are "LB_METHOD" (valid values
  include "GRAPH", "HYPERGRAPH"), "DEBUG_LEVEL" (valid values are
  0 to 10, default is 1), etc.
   */

  /**  compute_partitioning is an internal method that computes 
       a rebalanced partitioning for the data in the object
      that this class was constructed with.

      \param force_repartitioning Optional argument defaults to false. By
         default, compute_partitioning() only does anything the first time
         it is called, and subsequent repeated calls are no-ops. If the user's
         intent is to re-compute the partitioning (e.g., if parameters
         or other inputs have been changed), then setting this flag to
         true will force a new partitioning to be computed.
   */
void Partitioner::
compute_partitioning(bool force_repartitioning)
{

  bool use_zoltan = false;
  Teuchos::ParameterList& sublist = paramlist_;

  std::string partitioning_method_str("PARTITIONING_METHOD");
  std::string partitioning_method =
    paramlist_.get(partitioning_method_str, "UNSPECIFIED");

  std::string zoltan("ZOLTAN");

  if (alreadyComputed() && !force_repartitioning)
    return;

#ifdef HAVE_ISORROPIA_ZOLTAN
  if (partitioning_method != "SIMPLE_LINEAR") {
    use_zoltan = true;
  }

  if (use_zoltan) {
    if (input_graph_.get() != 0)
      lib_ = Teuchos::rcp(new ZoltanLibClass(input_graph_, costs_));
    else
      lib_ = Teuchos::rcp(new ZoltanLibClass(input_matrix_, costs_));
  }

  sublist = (paramlist_.sublist(zoltan));
#else /* HAVE_ISORROPIA_ZOLTAN */
  if (paramlist_.isSublist(zoltan)) {
    throw Isorropia::Exception("Zoltan requested, but Zoltan not enabled.");
  }
#endif /* HAVE_ISORROPIA_ZOLTAN */

  if (use_zoltan == false) {
    if (input_graph_.get() == 0)
      lib_ = Teuchos::rcp(new InternalPartitioner(input_matrix_, costs_));
    else
      lib_ = Teuchos::rcp(new InternalPartitioner(input_graph_, costs_));
  }

  lib_->repartition(paramlist_, myNewElements_, exports_, imports_);
  operation_already_computed_ = true;
}

void Partitioner::
compute(bool force_repartitioning)
{
  compute_partitioning(force_repartitioning);
}

  /** An internal method which determines whether the 
      method compute_partitioning() has already been
      called on this class instance.
  */
bool Partitioner::partitioning_already_computed() const {
  return (alreadyComputed());
}

  /** An internal method which returns the new partition ID for a given element that
     resided locally in the old partitioning.
  */
int Partitioner::newPartitionNumber(int myElem) const
{
  return (getNewPropertyOfElem(myElem));
}

  /** An internal method which returns the number of elements in a given partition.

      (Currently only implemented for the case where 'partition' is local.)
  */
int Partitioner::numElemsInPartition(int partition) const
{
  return (getNbrElemsWithProperty(partition));
}

  /** An internal method which fills caller-allocated list (of length len) with the
      global element ids to be located in the given partition.

      (Currently only implemented for the case where 'partition' is local.)
  */
void Partitioner::elemsInPartition(int partition, int* elementList, int len) const {
  return (getElemsWithProperty(partition, elementList, len));
}

} // namespace EPETRA

#endif //HAVE_EPETRA

}//namespace Isorropia

