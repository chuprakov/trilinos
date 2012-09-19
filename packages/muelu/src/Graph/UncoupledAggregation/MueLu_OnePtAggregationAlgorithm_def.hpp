// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * MueLu_OnePtAggregationAlgorithm_def.hpp
 *
 *  Created on: Sep 18, 2012
 *      Author: Tobias Wiesner
 */

#ifndef MUELU_ONEPTAGGREGATIONALGORITHM_DEF_HPP_
#define MUELU_ONEPTAGGREGATIONALGORITHM_DEF_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_OnePtAggregationAlgorithm.hpp"

#include "MueLu_Graph.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
OnePtAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::OnePtAggregationAlgorithm(RCP<const FactoryBase> const &graphFact)
{
}

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
LocalOrdinal OnePtAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildAggregates(Graph const & graph, Aggregates & aggregates, Teuchos::ArrayRCP<unsigned int> & aggStat, Teuchos::ArrayRCP<unsigned int> & coarse_aggStat) const {
  Monitor m(*this, "Coarsen Uncoupled (OnePtAggregationAlgorithm)");

  const LocalOrdinal nRows = graph.GetNodeNumVertices();
  const int myRank = graph.GetComm()->getRank();

  // vertex ids for output
  Teuchos::ArrayRCP<LocalOrdinal> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<LocalOrdinal> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);

  // some internal variables
  LocalOrdinal nLocalAggregates = aggregates.GetNumAggregates();    // number of local aggregates on current proc
  LocalOrdinal iNode1  = 0;        // current node

  // main loop over all local rows of grpah(A)
  while (iNode1 < nRows) {

  //std::cout << aggStat[iNode1] << std::endl;
    if (aggStat[iNode1] == NodeStats::ONEPT) {

      aggregates.SetIsRoot(iNode1);    // mark iNode1 as root node for new aggregate 'ag'
      Aggregate ag;
      ag.list.push_back(iNode1);
      ag.index = nLocalAggregates++;

      // prepare aggStat array for next coarser level
      coarse_aggStat[ag.index] = NodeStats::ONEPT;       // mark coarse node as "1pt aggregate"
      // finalize aggregate
      for(size_t k=0; k<ag.list.size(); k++) {
        aggStat[ag.list[k]] = NodeStats::AGGREGATED;
        vertex2AggId[ag.list[k]] = ag.index;
        procWinner[ag.list[k]] = myRank;
      }
    }

    iNode1++;
  } // end while

  // update aggregate object
  aggregates.SetNumAggregates(nLocalAggregates);

  // print aggregation information
  this->PrintAggregationInformation("Phase 0 (1pt aggregates):", graph, aggregates, aggStat);

  // collect some local information
  LO nLocalAggregated    = 0;
  LO nLocalNotAggregated = 0;
  for (LO i = 0; i < nRows; i++) {
    if (aggStat[i] == NodeStats::AGGREGATED) nLocalAggregated++;
    else nLocalNotAggregated++;
  }

  return nLocalNotAggregated;
}

} // end namespace


#endif /* MUELU_ONEPTAGGREGATIONALGORITHM_DEF_HPP_ */
