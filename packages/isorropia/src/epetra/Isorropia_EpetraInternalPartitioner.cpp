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

#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Isorropia_EpetraInternalPartitioner.hpp>


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

InternalPartitioner::InternalPartitioner(Teuchos::RefCountPtr<const Epetra_CrsGraph> input_graph):
  Library(input_graph)
 {
   setInputType("HYPERGRAPH");
 }

InternalPartitioner::InternalPartitioner(Teuchos::RefCountPtr<const Epetra_CrsGraph> input_graph,
			  Teuchos::RefCountPtr<CostDescriber> costs):
  Library(input_graph, costs)
{
   setInputType("HYPERGRAPH");
}

InternalPartitioner::InternalPartitioner(Teuchos::RefCountPtr<const Epetra_RowMatrix> input_matrix):
  Library(input_matrix)
{
   setInputType("HYPERGRAPH");
}

InternalPartitioner::InternalPartitioner(Teuchos::RefCountPtr<const Epetra_RowMatrix> input_matrix,
			  Teuchos::RefCountPtr<CostDescriber> costs):
  Library(input_matrix, costs)
{
   setInputType("HYPERGRAPH");
}

InternalPartitioner::~InternalPartitioner() {}

int InternalPartitioner::precompute()
{
  std::string str1("Isorropia::InternalPartitioner::precompute ");
  std::string str2;


  Library::precompute();

  return (0);
}




int InternalPartitioner::
repartition(Teuchos::ParameterList& paramList,
	    std::vector<int>& myNewElements,
	    std::map<int,int>& exports,
	    std::map<int,int>& imports)
{
  precompute();

  int nrows = input_map_->NumMyElements();

  if (nrows && costs_.get()){
    if (costs_->haveVertexWeights()){
      std::map<int, float> vwgts;
      costs_->getVertexWeights(vwgts);

      double *vals = new double [nrows];

      for (int i=0; i<nrows; i++){
	int gid = input_map_->GID(i);
	std::map<int, float>::iterator iter = vwgts.find(gid);
	if (iter == vwgts.end()){
	  throw Isorropia::Exception("error 1 in simple linear repartitioning");
	}
	vals[i] = (double)iter->second;
      }
      weights_ = Teuchos::rcp(new Epetra_Vector(Copy, *input_map_, vals));

      delete [] vals;
    }
  }

  if (nrows && !weights_.get()){
    if (input_graph_.get() != 0) {
      weights_ = Teuchos::rcp(create_row_weights_nnz(*input_graph_));
    }
    else {
      weights_ = Teuchos::rcp(create_row_weights_nnz(*input_matrix_));
    }
  }

  int err = Isorropia::Epetra::repartition(*input_map_,
					   *weights_,
					   myNewElements,
					   exports, imports);
  if (err != 0) {
    throw Isorropia::Exception("error 2 in simple linear repartitioning");
  }

  return (err);
}

int InternalPartitioner::
color(Teuchos::ParameterList& paramList,
      std::vector<int>& myNewElements)
{
  throw Isorropia::Exception("Coloring only available in Zoltan");
  return (-1);
}

int InternalPartitioner::
order(Teuchos::ParameterList& paramList,
      std::vector<int>& myNewElements)
{
  throw Isorropia::Exception("Ordering only available in Zoltan");
  return (-1);
}


int InternalPartitioner::postcompute()
{
  return (0);
};


} // namespace EPETRA

#endif //HAVE_EPETRA

}//namespace Isorropia

