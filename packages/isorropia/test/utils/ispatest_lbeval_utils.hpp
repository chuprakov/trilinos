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

#ifndef _ispatest_lbeval_utils_hpp_
#define _ispatest_lbeval_utils_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <vector>

#ifdef HAVE_EPETRA

#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>
#include <Isorropia_EpetraCostDescriber.hpp>

/*  Load balance evaluation (As calculated in Zoltan_LB_Eval)

  Given an Epetra graph and possible vertex and hyperedge weights,
  calculate the hypergraph balance, cutn and cutl.  We think of
  the rows of the graph as the objects (vertices) to be partitioned
  and the columns of the graph as the hyperedges.

  ====================================

  Definition of hypergraph balance:

  Suppose Size_p is target size of partition p.  (Sum of all Size_p is 1.0)
  (For now Size_p is always 1/p in isorropia - we don't have an interface
  to set different desired partition sizes.)

  Let W_p be the total row weight on process p and WT be the sum of the
  row weights on all processes.  Then

  imbalance on process p =   |W_p - Size_p*WT| / Size_p*WT

  balance = (1 + maximum imbalance over all processes p)

  ====================================

  Definition of hypergraph cutn and cutl:

  Suppose Cut_k is the number of cuts in column k.  (So column k is in 
  Cut_k + 1 different partitions.)  Suppose W_k is the weight of column k.

  cutn = Sum over all k of W_K * ((Cut_k > 0) ? 1 : 0)

  cutl = Sum over all k of Cut_k * W_k

  ====================================

  TODO explain metrics computed in compute_graph_metrics
*/
namespace ispatest {

/** Compute Zoltan-style hypergraph metrics given a partitioned
    CrsGraph and a CostDescriber (weight) object.
 */
int compute_hypergraph_metrics(const Epetra_CrsGraph &graph,
            Isorropia::Epetra::CostDescriber &costs,
            double &balance, double &cutn, double &cutl);

/** Compute Zoltan-style hypergraph metrics given a partitioned
    RowMatrix and a CostDescriber (weight) object.
 */
int compute_hypergraph_metrics(const Epetra_RowMatrix &matrix,
            Isorropia::Epetra::CostDescriber &costs,
            double &balance, double &cutn, double &cutl);

/** Compute Zoltan-style hypergraph metrics given a partitioned
    RowMatrix and a CostDescriber (weight) object.
 */
int compute_hypergraph_metrics(const Epetra_BlockMap &rowmap, 
            const Epetra_BlockMap &colmap,
            int numGlobalColumns,
            Isorropia::Epetra::CostDescriber &costs,
            double &balance, double &cutn, double &cutl);

/** Compute graph metrics given an Epetra_RowMatrix
  */
int compute_graph_metrics(const Epetra_RowMatrix &matrix,
            Isorropia::Epetra::CostDescriber &costs,
            double &balance, int &numCuts, double &cutWgt, double &cutn, double &cutl);

/** Compute graph metrics given an Epetra_CrsGraph
  */
int compute_graph_metrics(const Epetra_CrsGraph &graph,
            Isorropia::Epetra::CostDescriber &costs,
            double &balance, int &numCuts, double &cutWgt, double &cutn, double &cutl);

/** Compute graph metrics given a row map, a column map, and a vector with one
    element for each row.  The element is a vector containing the column local ID 
    for each non zero in that row.
  */
int compute_graph_metrics(const Epetra_BlockMap &rowmap,
                          const Epetra_BlockMap &colmap,
                          std::vector<std::vector<int> > &rows,
                          Isorropia::Epetra::CostDescriber &costs,
            double &balance, int &numCuts, double &cutWgt, double &cutn, double &cutl);


/** Print out a distributed RowMatrix.  This only works for small test
    matrices of 1s and 0s, and 10 or fewer processes.
  */

void show_matrix(const char *txt, const Epetra_RowMatrix &matrix, const Epetra_Comm &comm);

/** Print out a distributed CrsGraph.  This only works for small test
    matrices of 1s and 0s and 10 or fewer processes.
  */

void show_matrix(const char *txt, const Epetra_CrsGraph &graph, const Epetra_Comm &comm);

/** Print out a distributed LinearProblem.  This only works for small test
    matrices of 1s and 0s and 10 or fewer processes.
  */

void show_matrix(const char *txt, const Epetra_LinearProblem &problem, const Epetra_Comm &comm);

}//namespace ispatest

#endif //HAVE_EPTERA

#endif

