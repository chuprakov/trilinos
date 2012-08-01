/*
 * MueLu_CheapAggregationAlgorithm_def.hpp
 *
 *  Created on: Jul 25, 2012
 *      Author: wiesner
 */

#ifndef MUELU_CHEAPAGGREGATIONALGORITHM_DEF_HPP_
#define MUELU_CHEAPAGGREGATIONALGORITHM_DEF_HPP_

#include <queue>

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_CheapAggregationAlgorithm_decl.hpp"

#include "MueLu_Graph.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_LinkedList.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
CheapAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CheapAggregationAlgorithm(RCP<const FactoryBase> const &graphFact)
: ordering_(NATURAL), minNodesPerAggregate_(3), maxNeighAlreadySelected_(0)
  { }

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
LocalOrdinal CheapAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Phase1(Graph const & graph, Aggregates & aggregates, Teuchos::ArrayRCP<NodeState> & aggStat) const {
  Monitor m(*this, "Coarsen Uncoupled (CheapAggregationAlgorithm): Phase 1");

  RCP<const Map> gBoundaryNodeMap = graph.GetBoundaryNodeMap();

  std::string orderingType;
  switch (ordering_) {
  case NATURAL:
    orderingType = "Natural";
    break;
  case RANDOM:
      orderingType = "Random";
      break;
  case GRAPH:
    orderingType = "Graph";
    break;
  default:
    break;
  }

  const LocalOrdinal nRows = graph.GetNodeNumVertices();

  // vertex ids for output
  Teuchos::ArrayRCP<LocalOrdinal> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<LocalOrdinal> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);

  // some internal variables
  LocalOrdinal nLocalAggregates = 0;   // number of local aggregates on current proc
  std::queue<LocalOrdinal> graph_ordering_inodes; // inodes for graph ordering
  LocalOrdinal iNode2 = 0;        // local iteration variable
  LocalOrdinal iNode1  = 0;        // current node
  Teuchos::ArrayRCP<LO> randomVector;
  
  if ( ordering_ == RANDOM )       /* random ordering */
    {
      //TODO: could be stored in a class that respect interface of LinkedList

      randomVector = Teuchos::arcp<LO>(nRows); //size_t or int ?-> to be propagated
      for (my_size_t i = 0; i < nRows; ++i) randomVector[i] = i;
      RandomReorder(randomVector);
    } 

  // main loop over all local rows of grpah(A)
  while (iNode2 < nRows) {

    // pick the next node to aggregate
    if      (ordering_ == NATURAL) iNode1 = iNode2++;
    else if (ordering_ == RANDOM ) iNode1 = randomVector[iNode2++];
    else if (ordering_ == GRAPH) {
      // if there are no nodes for graph ordering scheme
      if(graph_ordering_inodes.size() == 0) {
        // add exactly one ready node for graph ordering aggregates
        for(LocalOrdinal jnode=0; jnode<nRows; jnode++) {
          if(aggStat[jnode] == READY) {
            graph_ordering_inodes.push(jnode);
            break;
          }
        }
      }
      if(graph_ordering_inodes.size()==0) break; // there's no ready node any more -> end phase 1
      iNode1 = graph_ordering_inodes.front(); // take next node from graph ordering queue
      graph_ordering_inodes.pop();           // delete this node in list
    }

    // consider iNode1 only if it is in READY state
    if(aggStat[iNode1] == READY) {
      // build new aggregate
      Aggregate ag;
      ag.list.push_back(iNode1);

      // extract column information from graph for current row on current proc
      Teuchos::ArrayView<const LocalOrdinal> neighOfINode = graph.getNeighborVertices(iNode1);

      LocalOrdinal cnt_neighbours = 0;  // number of possible neighbour nodes for current new aggregate

      // build tentative new aggregate
      for( typename Teuchos::ArrayView<const LocalOrdinal>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
        LocalOrdinal index = *it;

        // note: this is uncoupled coarsening
        // only column indices of current proc are allowed
        if(graph.isLocalNeighborVertex(index)) {
          // check status of current neighbor node
          if(aggStat[index]==READY ||
              aggStat[index]==NOTSEL ) {
            ag.list.push_back(index); // add neighbor node to current aggregate
          }
          else if(aggStat[index]!=BDRY){
            cnt_neighbours++;
          }
        } // end if: current column index belongs to this proc
      } // end for: loop over all columns in current row

      // if there are too many neighbours aggregated or the number of nodes
      // in the new aggregate is too few, don't do this one.

      // check if aggregate ag is acceptable
      if((cnt_neighbours > GetMaxNeighAlreadySelected() ) || // aggregate is too close to bdry nodes
          (ag.list.size() < (unsigned int) GetMinNodesPerAggregate())) {      // not enough nodes in new aggregate
        // failed to build a new aggregate
        ag.list.clear();
        aggStat[iNode1] = NOTSEL; // node not valid to be supernode, mark it as NOTSEL
        if(ordering_ == GRAPH) {
          // even though the aggregate around iNode1 is not perfect, we try the ndoes where iNode1 is connected to
          // loop over all column indices
          for (typename Teuchos::ArrayView<const LocalOrdinal>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
            LocalOrdinal index = *it;
            if(graph.isLocalNeighborVertex(index) && aggStat[index] == READY) // if node connected to iNode1 is not aggregated
              graph_ordering_inodes.push(index);
          }
        }
      } else {
        // accept new aggregate
        aggregates.SetIsRoot(iNode1);    // mark iNode1 as root node for new aggregate 'ag'
        ag.index = nLocalAggregates++;       // aggregate accepted, increment aggregate counter
        vertex2AggId[iNode1] = ag.index;
        procWinner[iNode1] = graph.GetComm()->getRank();
        std::cout << "build new aggregate of size " << ag.list.size() << " nodes" << std::endl;
        std::cout << "nodes: ";
        for (unsigned int k=0; k<ag.list.size(); k++)
          std::cout << ag.list[k] << " ";
        std::cout << std::endl;

        if(ag.list.size() == 2) {
          std::cout << "!!!!!!!!!!!!!!!!!!!!!! built a 1pt aggregate !!!!!!!!!!!!!!!!!!"<< std::endl;
        }
        for (unsigned int k=0; k<ag.list.size(); k++) {
          aggStat[ag.list[k]] = SELECTED;
          vertex2AggId[ag.list[k]] = ag.index;
          procWinner[ag.list[k]] = graph.GetComm()->getRank();
          if(ordering_ == GRAPH) {
            Teuchos::ArrayView<const LocalOrdinal> neighOfJNode = graph.getNeighborVertices(ag.list[k]);
            for(typename Teuchos::ArrayView<const LocalOrdinal>::const_iterator it = neighOfJNode.begin(); it!=neighOfJNode.end(); ++it) {
              LocalOrdinal index = *it;
              if(graph.isLocalNeighborVertex(index) && aggStat[index] == READY)
                graph_ordering_inodes.push(index);
            }
          } // end GRAPH
        } // end loop over all nodes in aggregate
      } // end if accept aggs or decline aggs


    } // end if aggStat[iNode1] == READY

  } // end while

  // update aggregate object
  aggregates.SetNumAggregates(nLocalAggregates);

  // clean up
  if(graph_ordering_inodes.size()>0){
    for(unsigned int k=0; k<graph_ordering_inodes.size(); k++)
      graph_ordering_inodes.pop();
  }

  // verbose
  {
    const RCP<const Teuchos::Comm<int> > & comm = graph.GetComm();

    if(IsPrint(Warnings0)) {
      GO localReady  = 0;
      GO globalReady = 0;

      // compute number of local nodes with status "ready"
      for(LO i=0; i<nRows; ++i)
        if(aggStat[i] == READY) localReady++;
      // compute global number of nodes with status "ready"
      sumAll(comm, (GO)localReady, globalReady);
      if(globalReady > 0)
        GetOStream(Warnings0,0) << "Aggregation (UC): Phase 1 (WARNING) " << globalReady << " unaggregated nodes left (status READY)" << std::endl;
    }

    if(IsPrint(Statistics1)) {
      LO localSelected  = 0;
      GO globalSelected = 0;
      GO globalNRows    = 0;
      for(LO i=0; i<nRows; ++i)
        if(aggStat[i] == SELECTED) localSelected++;
      sumAll(comm, (GO)localSelected, globalSelected);
      sumAll(comm, (GO)nRows, globalNRows);
      GetOStream(Statistics1, 0) << "Aggregation (UC): Phase 1: Nodes aggregated = " << globalSelected << " out of " << globalNRows << " nodes" << std::endl;
      GO nAggregatesGlobal = 0;
      sumAll(comm, (GO)nLocalAggregates, nAggregatesGlobal);
      GetOStream(Statistics1, 0) << "Aggregation (UC): Phase 1: Total aggregates = " << nAggregatesGlobal << std::endl;
    }
  }

  // collect some local information
  LO nLocalSelected    = 0;
  LO nLocalBdry        = 0;
  LO nLocalNotSelected = 0;
  LO nLocalReady       = 0;
  for (LO i = 0; i < nRows; i++) {
    if      (aggStat[i] == SELECTED) nLocalSelected++;
    else if (aggStat[i] == BDRY)     nLocalBdry++;
    else if (aggStat[i] == NOTSEL)   nLocalNotSelected++;
    else if (aggStat[i] == READY)    nLocalReady++;
  }

  return nLocalReady + nLocalNotSelected;
}


template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
LocalOrdinal CheapAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Phase1a(Graph const & graph, Aggregates & aggregates, Teuchos::ArrayRCP<NodeState> & aggStat, Teuchos::ArrayRCP<NodeState> & coarse_aggStat) const {
  Monitor m(*this, "Coarsen Uncoupled (CheapAggregationAlgorithm): Phase 1a");

  RCP<const Map> gBoundaryNodeMap = graph.GetBoundaryNodeMap();

  std::string orderingType;
  switch (ordering_) {
  case NATURAL:
    orderingType = "Natural";
    break;
  case RANDOM:
    orderingType = "Random";
    break;
  case GRAPH:
    orderingType = "Graph";
    break;
  default:
    break;
  }

  const LocalOrdinal nRows = graph.GetNodeNumVertices();

  // vertex ids for output
  Teuchos::ArrayRCP<LocalOrdinal> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<LocalOrdinal> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);

  // some internal variables
  LocalOrdinal nLocalAggregates = 0;   // number of local aggregates on current proc
  std::queue<LocalOrdinal> graph_ordering_inodes; // inodes for graph ordering
  LocalOrdinal iNode2 = 0;        // local iteration variable
  LocalOrdinal iNode1  = 0;        // current node
  Teuchos::ArrayRCP<LO> randomVector;
  
  if ( ordering_ == RANDOM )       /* random ordering */
    {
      //TODO: could be stored in a class that respect interface of LinkedList

      randomVector = Teuchos::arcp<LO>(nRows); //size_t or int ?-> to be propagated
      for (my_size_t i = 0; i < nRows; ++i) randomVector[i] = i;
      RandomReorder(randomVector);
    } 
    
  // main loop over all local rows of grpah(A)
  while (iNode2 < nRows) {

    // pick the next node to aggregate
    if      (ordering_ == NATURAL) iNode1 = iNode2++;
    else if (ordering_ == RANDOM ) iNode1 = randomVector[iNode2++];
    else if (ordering_ == GRAPH) {
      // if there are no nodes for graph ordering scheme
      if(graph_ordering_inodes.size() == 0) {
        // add exactly one ready node for graph ordering aggregates
        for(LocalOrdinal jnode=0; jnode<nRows; jnode++) {
          if(aggStat[jnode] == READY ||
              aggStat[jnode] == READY_1PT) {
            graph_ordering_inodes.push(jnode);
            break;
          }
        }
      }
      if(graph_ordering_inodes.size()==0) break; // there's no ready node any more -> end phase 1
      iNode1 = graph_ordering_inodes.front(); // take next node from graph ordering queue
      graph_ordering_inodes.pop();           // delete this node in list
    }

    // consider iNode1 only if it is in READY state
    if(aggStat[iNode1] == READY) {
      // build new aggregate
      Aggregate ag;
      ag.list.push_back(iNode1);

      // extract column information from graph for current row on current proc
      Teuchos::ArrayView<const LocalOrdinal> neighOfINode = graph.getNeighborVertices(iNode1);

      LocalOrdinal cnt_neighbours = 0;  // number of possible neighbour nodes for current new aggregate

      // build tentative new aggregate
      for( typename Teuchos::ArrayView<const LocalOrdinal>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
        // note: this is uncoupled coarsening
        // only column indices of current proc are allowed
        if(graph.isLocalNeighborVertex(*it)) {
          // check status of current neighbor node
          if(aggStat[*it]==READY ||
              aggStat[*it]==NOTSEL ) {
            ag.list.push_back(*it); // add neighbor node to current aggregate
          } else if(aggStat[*it]!=READY_1PT) {
            cnt_neighbours++;
          }
        } // end if: current column index belongs to this proc
      } // end for: loop over all columns in current row

      // if there are too many neighbours aggregated or the number of nodes
      // in the new aggregate is too few, don't do this one.

      // check if aggregate ag is acceptable
      if((cnt_neighbours > GetMaxNeighAlreadySelected() ) || // aggregate is too close to bdry nodes
         (ag.list.size() < (unsigned int) GetMinNodesPerAggregate())) {      // not enough nodes in new aggregate
        // failed to build a new aggregate
        ag.list.clear();
        aggStat[iNode1] = NOTSEL; // node not valid to be supernode, mark it as NOTSEL
        if(ordering_ == GRAPH) {
          // even though the aggregate around iNode1 is not perfect, we try the ndoes where iNode1 is connected to
          // loop over all column indices
          for (typename Teuchos::ArrayView<const LocalOrdinal>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
            if(graph.isLocalNeighborVertex(*it) &&
                (aggStat[*it] == READY ||
                 aggStat[*it] == READY_1PT)) // if node connected to iNode1 is not aggregated
              graph_ordering_inodes.push(*it);
          }
        }
      } else {
        // accept new aggregate
        aggregates.SetIsRoot(iNode1);    // mark iNode1 as root node for new aggregate 'ag'
        ag.index = nLocalAggregates++;       // aggregate accepted, increment aggregate counter
        vertex2AggId[iNode1] = ag.index;
        procWinner[iNode1] = graph.GetComm()->getRank();
        /*std::cout << "build new aggregate of size " << ag.list.size() << " nodes" << std::endl;
        std::cout << "nodes: ";
        for (unsigned int k=0; k<ag.list.size(); k++)
          std::cout << ag.list[k] << " ";
        std::cout << std::endl;*/

        if(ag.list.size() == 2) {
          std::cout << "!!!!!!!!!!!!!!!!!!!!!! built a 1pt aggregate !!!!!!!!!!!!!!!!!!"<< std::endl;
        }
        for (unsigned int k=0; k<ag.list.size(); k++) {
          aggStat[ag.list[k]] = SELECTED;
          coarse_aggStat[ag.index] = READY; // mark aggregate id to be a valid READY node on the next coarser grid
          vertex2AggId[ag.list[k]] = ag.index;
          procWinner[ag.list[k]] = graph.GetComm()->getRank();
          if(ordering_ == GRAPH) {
            Teuchos::ArrayView<const LocalOrdinal> neighOfJNode = graph.getNeighborVertices(ag.list[k]);
            for(typename Teuchos::ArrayView<const LocalOrdinal>::const_iterator it = neighOfJNode.begin(); it!=neighOfJNode.end(); ++it) {
              if(graph.isLocalNeighborVertex(*it) &&
                  (aggStat[*it] == READY||
                   aggStat[*it] == READY_1PT))
                graph_ordering_inodes.push(*it);
            }
          } // end GRAPH
        } // end loop over all nodes in aggregate
      } // end if accept aggs or decline aggs


    } // end if aggStat[iNode1] == READY
    else if (aggStat[iNode1] == READY_1PT) {
      // TODO
      aggregates.SetIsRoot(iNode1);    // mark iNode1 as root node for new aggregate 'ag'
      Aggregate ag;
      ag.list.push_back(iNode1);
      ag.index = nLocalAggregates++;

      coarse_aggStat[ag.index] = READY_1PT;

      // finalize aggregate
      for(size_t k=0; k<ag.list.size(); k++) {
        aggStat[ag.list[k]] = SELECTED_1PT;
        vertex2AggId[ag.list[k]] = ag.index;
        procWinner[ag.list[k]] = graph.GetComm()->getRank();
      }
    }

  } // end while

  // update aggregate object
  aggregates.SetNumAggregates(nLocalAggregates);

  // clean up
  if(graph_ordering_inodes.size()>0){
    for(unsigned int k=0; k<graph_ordering_inodes.size(); k++)
      graph_ordering_inodes.pop();
  }

  // verbose
  {
    const RCP<const Teuchos::Comm<int> > & comm = graph.GetComm();

    if(IsPrint(Warnings0)) {
      GO localReady  = 0;
      GO globalReady = 0;

      // compute number of local nodes with status "ready"
      for(LO i=0; i<nRows; ++i)
        if(aggStat[i] == READY || aggStat[i] == READY_1PT) localReady++;
      // compute global number of nodes with status "ready"
      sumAll(comm, (GO)localReady, globalReady);
      if(globalReady > 0)
        GetOStream(Warnings0,0) << "Aggregation (UC): Phase 1 (WARNING) " << globalReady << " unaggregated nodes left (status READY)" << std::endl;
    }

    if(IsPrint(Statistics1)) {
      LO localSelected  = 0;
      GO globalSelected = 0;
      GO globalNRows    = 0;
      for(LO i=0; i<nRows; ++i)
        if(aggStat[i] == SELECTED || aggStat[i] == SELECTED_1PT) localSelected++;
      sumAll(comm, (GO)localSelected, globalSelected);
      sumAll(comm, (GO)nRows, globalNRows);
      GetOStream(Statistics1, 0) << "Aggregation (UC): Phase 1: Nodes aggregated = " << globalSelected << " out of " << globalNRows << " nodes" << std::endl;
      GO nAggregatesGlobal = 0;
      sumAll(comm, (GO)nLocalAggregates, nAggregatesGlobal);
      GetOStream(Statistics1, 0) << "Aggregation (UC): Phase 1: Total aggregates = " << nAggregatesGlobal << std::endl;
    }
  }

  // collect some local information
  LO nLocalSelected    = 0;
  LO nLocalBdry        = 0;
  LO nLocalNotSelected = 0;
  LO nLocalReady       = 0;
  for (LO i = 0; i < nRows; i++) {
    if      (aggStat[i] == SELECTED || aggStat[i] == SELECTED_1PT ) nLocalSelected++;
    else if (aggStat[i] == BDRY)     nLocalBdry++;
    else if (aggStat[i] == NOTSEL)   nLocalNotSelected++;
    else if (aggStat[i] == READY || aggStat[i] == READY_1PT )    nLocalReady++;
  }

  return nLocalReady + nLocalNotSelected;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
LocalOrdinal CheapAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Phase2_maxlink(Graph const & graph, Aggregates & aggregates, Teuchos::ArrayRCP<NodeState> & aggStat) const {
  Monitor m(*this, "Coarsen Uncoupled (CheapAggregationAlgorithm): Phase 2 [max_link]");

  // vertex ids for output
  Teuchos::ArrayRCP<LocalOrdinal> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<LocalOrdinal> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);

  const LocalOrdinal nRows = graph.GetNodeNumVertices();

  // loop over all local rows
  for (LocalOrdinal iNode=0; iNode<nRows; iNode++) {
    if(aggStat[iNode] == NOTSEL ||  // this is a non-aggregated node
        aggStat[iNode] == READY) {  // this is not good. There should be no ready nodes any more
      LocalOrdinal selected_aggregate = -1;

      std::map<LocalOrdinal,LocalOrdinal> aggid2cntconnections;

      Teuchos::ArrayView<const LocalOrdinal> neighOfINode = graph.getNeighborVertices(iNode);
      for(typename Teuchos::ArrayView<const LocalOrdinal>::const_iterator it = neighOfINode.begin(); it!=neighOfINode.end(); ++it) {
        LocalOrdinal index = *it;
        if(graph.isLocalNeighborVertex(index) && aggStat[index] == SELECTED) {
          LocalOrdinal aggid = vertex2AggId[index]; // get (local) aggregate id
          if(aggid2cntconnections.count(aggid)>0) aggid2cntconnections[aggid] += 1;
          else aggid2cntconnections[aggid] = 1;
        } // if SELECTED

      } // loop over all columns

      // find aggregate id with most connections
      LocalOrdinal maxcnt = 0;
      for(typename std::map<LocalOrdinal,LocalOrdinal>::const_iterator it=aggid2cntconnections.begin(); it!=aggid2cntconnections.end();++it) {
        if(maxcnt < it->second) {
          maxcnt = it->second;
          selected_aggregate = it->first;
        }
      }

      // add node iNode to aggregate
      if(selected_aggregate != -1) {
        aggStat[iNode] = SELECTED;
        vertex2AggId[iNode] = selected_aggregate;
        procWinner[iNode] = graph.GetComm()->getRank();
      }
    } // end if aggState == NOTSEL...
  } // end loop over all local rows

  // verbose
  {
    const RCP<const Teuchos::Comm<int> > & comm = graph.GetComm();

    if(IsPrint(Warnings0)) {
      GO localReady  = 0;
      GO globalReady = 0;

      // compute number of local nodes with status "ready"
      for(LO i=0; i<nRows; ++i)
        if(aggStat[i] == READY) localReady++;
      // compute global number of nodes with status "ready"
      sumAll(comm, (GO)localReady, globalReady);
      if(globalReady > 0)
        GetOStream(Warnings0,0) << "Aggregation (UC): Phase 2 [max_link] (WARNING) " << globalReady << " unaggregated nodes left (status READY)" << std::endl;
    }

    if(IsPrint(Statistics1)) {
      LO localSelected  = 0;
      GO globalSelected = 0;
      GO globalNRows    = 0;
      for(LO i=0; i<nRows; ++i)
        if(aggStat[i] == SELECTED) localSelected++;
      sumAll(comm, (GO)localSelected, globalSelected);
      sumAll(comm, (GO)nRows, globalNRows);
      GetOStream(Statistics1, 0) << "Aggregation (UC): Phase 2 [max_link]: Nodes aggregated = " << globalSelected << " out of " << globalNRows << " nodes" << std::endl;
    }
  }

  // collect some local information
  LO nLocalSelected    = 0;
  LO nLocalBdry        = 0;
  LO nLocalNotSelected = 0;
  LO nLocalReady       = 0;
  for (LO i = 0; i < nRows; i++) {
    if      (aggStat[i] == SELECTED || aggStat[i] == SELECTED_1PT ) nLocalSelected++;
    else if (aggStat[i] == BDRY)     nLocalBdry++;
    else if (aggStat[i] == NOTSEL)   nLocalNotSelected++;
    else if (aggStat[i] == READY || aggStat[i] == READY_1PT )    nLocalReady++;
  }

  return nLocalReady + nLocalNotSelected;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
LocalOrdinal CheapAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Phase3(Graph const & graph, Aggregates & aggregates, Teuchos::ArrayRCP<NodeState> & aggStat) const {
  Monitor m(*this, "Coarsen Uncoupled (CheapAggregationAlgorithm): Phase 3");

  // form new aggregates from non-aggregated nodes

  // vertex ids for output
  Teuchos::ArrayRCP<LocalOrdinal> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<LocalOrdinal> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);

  const LocalOrdinal nRows = graph.GetNodeNumVertices();
  LocalOrdinal nLocalAggregates = aggregates.GetNumAggregates(); // return number of local aggregates on current proc

  // loop over all local rows
  for (LocalOrdinal iNode=0; iNode<nRows; iNode++) {
    if(aggStat[iNode] == NOTSEL ||  // this is a non-aggregated node
        aggStat[iNode] == READY) {  // this is not good. There should be no ready nodes any more
      Aggregate ag;
      ag.list.push_back(iNode);
      ag.index = nLocalAggregates++;

      Teuchos::ArrayView<const LocalOrdinal> neighOfINode = graph.getNeighborVertices(iNode);
      for(typename Teuchos::ArrayView<const LocalOrdinal>::const_iterator it = neighOfINode.begin(); it!=neighOfINode.end(); ++it) {
        LocalOrdinal index = *it;
        if(graph.isLocalNeighborVertex(index) && aggStat[index] != SELECTED && aggStat[index] != BDRY) {
          ag.list.push_back(index);
        } // if !SELECTED && !BDRY
      } // loop over all columns

      // finalize aggregate
      for(size_t k=0; k<ag.list.size(); k++) {
        aggStat[ag.list[k]] = SELECTED;
        vertex2AggId[ag.list[k]] = ag.index;
        procWinner[ag.list[k]] = graph.GetComm()->getRank();
      }
    } // end if NOTSEL
  }   // end for

  // update aggregate object
  aggregates.SetNumAggregates(nLocalAggregates);

  // verbose
  {
    const RCP<const Teuchos::Comm<int> > & comm = graph.GetComm();

    if(IsPrint(Warnings0)) {
      GO localReady  = 0;
      GO globalReady = 0;

      // compute number of local nodes with status "ready"
      for(LO i=0; i<nRows; ++i)
        if(aggStat[i] == READY) localReady++;
      // compute global number of nodes with status "ready"
      sumAll(comm, (GO)localReady, globalReady);
      if(globalReady > 0)
        GetOStream(Warnings0,0) << "Aggregation (UC): Phase 3 (WARNING) " << globalReady << " unaggregated nodes left (status READY)" << std::endl;
    }

    if(IsPrint(Statistics1)) {
      LO localSelected  = 0;
      GO globalSelected = 0;
      GO globalNRows    = 0;
      for(LO i=0; i<nRows; ++i)
        if(aggStat[i] == SELECTED) localSelected++;
      sumAll(comm, (GO)localSelected, globalSelected);
      sumAll(comm, (GO)nRows, globalNRows);
      GetOStream(Statistics1, 0) << "Aggregation (UC): Phase 3: Nodes aggregated = " << globalSelected << " out of " << globalNRows << " nodes" << std::endl;
      GO nAggregatesGlobal = 0;
      sumAll(comm, (GO)nLocalAggregates, nAggregatesGlobal);
      GetOStream(Statistics1, 0) << "Aggregation (UC): Phase 3: Total aggregates = " << nAggregatesGlobal << std::endl;
    }
  }

  // collect some local information
  LO nLocalSelected    = 0;
  LO nLocalBdry        = 0;
  LO nLocalNotSelected = 0;
  LO nLocalReady       = 0;
  for (LO i = 0; i < nRows; i++) {
    if      (aggStat[i] == SELECTED || aggStat[i] == SELECTED_1PT ) nLocalSelected++;
    else if (aggStat[i] == BDRY)     nLocalBdry++;
    else if (aggStat[i] == NOTSEL)   nLocalNotSelected++;
    else if (aggStat[i] == READY || aggStat[i] == READY_1PT )    nLocalReady++;
  }

  return nLocalReady + nLocalNotSelected;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
LocalOrdinal CheapAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Phase4(Graph const & graph, Aggregates & aggregates, Teuchos::ArrayRCP<NodeState> & aggStat) const {
  Monitor m(*this, "Coarsen Uncoupled (CheapAggregationAlgorithm): Phase 4");

  // vertex ids for output
  Teuchos::ArrayRCP<LocalOrdinal> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<LocalOrdinal> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);

  const LocalOrdinal nRows = graph.GetNodeNumVertices();
  LocalOrdinal nLocalAggregates = aggregates.GetNumAggregates(); // return number of local aggregates on current proc

  // loop over all local rows
  for (LocalOrdinal iNode=0; iNode<nRows; iNode++) {
    if(aggStat[iNode] == BDRY ) {
      Aggregate ag;
      ag.list.push_back(iNode);
      ag.index = nLocalAggregates++;

      Teuchos::ArrayView<const LocalOrdinal> neighOfINode = graph.getNeighborVertices(iNode);
      for(typename Teuchos::ArrayView<const LocalOrdinal>::const_iterator it = neighOfINode.begin(); it!=neighOfINode.end(); ++it) {
        LocalOrdinal index = *it;
        if(graph.isLocalNeighborVertex(index) && aggStat[index] != SELECTED && aggStat[index] != BDRY) {
          ag.list.push_back(index);
        } // if !SELECTED && !BDRY
      } // loop over all columns

      // finalize aggregate
      for(size_t k=0; k<ag.list.size(); k++) {
        aggStat[ag.list[k]] = SELECTED;
        vertex2AggId[ag.list[k]] = ag.index;
        procWinner[ag.list[k]] = graph.GetComm()->getRank();
      }
    } // end if BDRY
  }   // end for

  // update aggregate object
  aggregates.SetNumAggregates(nLocalAggregates);

  // verbose
  {
    const RCP<const Teuchos::Comm<int> > & comm = graph.GetComm();

    if(IsPrint(Warnings0)) {
      GO localReady  = 0;
      GO globalReady = 0;

      // compute number of local nodes with status "ready"
      for(LO i=0; i<nRows; ++i)
        if(aggStat[i] == READY) localReady++;
      // compute global number of nodes with status "ready"
      sumAll(comm, (GO)localReady, globalReady);
      if(globalReady > 0)
        GetOStream(Warnings0,0) << "Aggregation (UC): Phase 4 (WARNING) " << globalReady << " unaggregated nodes left (status READY)" << std::endl;
    }

    if(IsPrint(Statistics1)) {
      LO localSelected  = 0;
      GO globalSelected = 0;
      GO globalNRows    = 0;
      for(LO i=0; i<nRows; ++i)
        if(aggStat[i] == SELECTED) localSelected++;
      sumAll(comm, (GO)localSelected, globalSelected);
      sumAll(comm, (GO)nRows, globalNRows);
      GetOStream(Statistics1, 0) << "Aggregation (UC): Phase 4: Nodes aggregated = " << globalSelected << " out of " << globalNRows << " nodes" << std::endl;
      GO nAggregatesGlobal = 0;
      sumAll(comm, (GO)nLocalAggregates, nAggregatesGlobal);
      GetOStream(Statistics1, 0) << "Aggregation (UC): Phase 4: Total aggregates = " << nAggregatesGlobal << std::endl;
    }
  }

  // collect some local information
  LO nLocalSelected    = 0;
  LO nLocalBdry        = 0;
  LO nLocalNotSelected = 0;
  LO nLocalReady       = 0;
  for (LO i = 0; i < nRows; i++) {
    if      (aggStat[i] == SELECTED || aggStat[i] == SELECTED_1PT ) nLocalSelected++;
    else if (aggStat[i] == BDRY)     nLocalBdry++;
    else if (aggStat[i] == NOTSEL)   nLocalNotSelected++;
    else if (aggStat[i] == READY || aggStat[i] == READY_1PT )    nLocalReady++;
  }

  return nLocalReady + nLocalNotSelected;
}

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  void CheapAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RandomReorder(Teuchos::ArrayRCP<LO> list) const {
    //TODO: replace int
    int n = list.size();
    for(int i=0; i<n-1; i++) {
      std::swap(list[i], list[RandomOrdinal(i,n-1)]);
    }
  } 

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  int CheapAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RandomOrdinal(int min, int max) const {
    return min + static_cast<int>((max-min+1) * (static_cast<double>(std::rand()) / (RAND_MAX + 1.0)));
  }

} // end namespace


#endif /* MUELU_CHEAPAGGREGATIONALGORITHM_DEF_HPP_ */
