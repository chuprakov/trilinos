/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person    */
/* and disclaimer.                                                           */
/* ************************************************************************* */
/* ************************************************************************* */
/* ************************************************************************* */
/* Functions to aggregate vertices in a graph                                */
/* ************************************************************************* */
/* Inital author : Ray Tuminaro                                              */
/* Date          : Jan, 2010                                                 */
/* ************************************************************************* */

#define CLEAN_DEBUG

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <vector>

#include <Cthulhu_VectorFactory.hpp>
//#include <Cthulhu_ConfigDefs.hpp> // CombineMode

#include "MueLu_AggregationOptions.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Graph.hpp"

#include "MueLu_AggAlgorithm2Comm.hpp"

// MPI helper
#define sumAll(rcpComm, in, out) \
  Teuchos::reduceAll<int>(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));
#define minAll(rcpComm, in, out) \
  Teuchos::reduceAll<int>(*rcpComm, Teuchos::REDUCE_MIN, in, Teuchos::outArg(out));
#define maxAll(rcpComm, in, out) \
  Teuchos::reduceAll<int>(*rcpComm, Teuchos::REDUCE_MAX, in, Teuchos::outArg(out));

using namespace std;
using namespace MueLu;
using Teuchos::ArrayView;
typedef ArrayView<const int>::const_iterator iter;

int MueLu_PrintLevel();

#define MUELU_NOSCORE       -100  /* indicates that a quality score has not  */
                                  /* yet been assigned when determining to   */
                                  /* which existing aggregate a vertex       */
                                  /* should be assigned.                     */

#define MUELU_DISTONE_VERTEX_WEIGHT 100  /* Weights associated with all      */
                                  /* vertices that have a direct connection  */
                                  /* to the aggregate root.                  */

#define INCR_SCALING 3            /* Determines how much of a penalty should */
                                  /* be deduced from a score during Phase 5  */
                                  /* for each Phase 5 vertex already added   */
                                  /* to this aggregate. Specifically the     */
                                  /* penalty associated with aggregate y is  */
                                  /*   max (INCR_SCALING*NNewVtx,            */
                                  /*        UnpenalizedScore*(1-             */
                                  /*              MUELU_PENALTYFACTOR))      */
                                  /* where NNewVtx is the number of phase 5  */
                                  /* vertices already assigned to y.         */

#define MUELU_PENALTYFACTOR .30   /* Determines maximum allowable            */
                                  /* percentage of a score that can be       */
                                  /* deducted from this score for having     */
                                  /* already enlargened an aggregate to      */
                                  /* which we are contemplated adding another*/
                                  /* vertex.  Should be between 0 and 1.     */

int MueLu_RootCandidates(int nVertices, Teuchos::ArrayRCP<int> vertex2AggId, const Graph<int,int> graph,
                         Teuchos::ArrayRCP<int> &candidates, int &nCandidates, int &nCandidatesGlobal);

int MueLu_RemoveSmallAggs(Aggregates<int,int> & aggregates, int min_size,
                          RCP<Cthulhu::Vector<double> > & weights_, const AggAlgorithm2Comm & myWidget);

// Take a partially aggregated graph and complete the aggregation. This is
// typically needed to take care of vertices that are left over after
// creating a bunch of ideal aggregates (root plus immediate neighbors).
//
// On input, the structure Aggregates describes already aggregated vertices.
// The field procWinners[] indicates the processor owning the aggregate to
// which a vertex is "definitively" assigned. If on entry 
// procWinners[i] == MUELU_UNASSIGNED, MueLu_ArbitrateAndCommunicate() 
// will arbitrate and decide which processor's aggregate really has
// the vertex. If only one processor claims ownership (as in
// the Uncoupled case), no real arbitration is needed. Otherwise,
// random arbitration is done.
//
// This cleanup has many phases:
//   
//   Phase 1b: Invoke MueLu_ArbitrateAndCommunicate() to ensure that
//             all processors have the same view of aggregated vertices
//             (e.g., to which aggregate they have been assigend and
//             which processor owns that aggregate).
//   Phase 2:  Check for vertices (local or nonlocal) which are Adjacent
//             to root nodes. Tentatively assign these to the aggregate
//             associated with the root. Arbitrate any cases where 
//             several processors claim the same vertex for one of 
//             their aggregates via MueLu_ArbitrateAndCommunicate().
//   Phase 3:  Try to create new aggregates if it looks like there are
//             root node candidates which have many unaggregated neighbors.
//             This decision to make a new aggregate is based on only local
//             information. However, the new aggregate will be tentatively
//             assigned any unaggregated ghost vertices. Arbitration is again
//             done by MueLu_ArbitrateAndCommunicate() where local vertices
//             use a weight[] = 2 and ghost vertices have weight[] = 1.
//             The basic idea is that after arbitration, each aggregate
//             is guaranteed to keep all local vertices assigned in
//             this phase. Thus, by basing the aggregation creation logic 
//             on local information, we are guarantee to have a sufficiently
//             large aggregation. The only local vertices that might be
//             assigned to another processor's aggregates are unclaimed
//             during this phase of the aggregation.
//   Phase 5:  Sweep new points into existing aggregates. Each processor tries
//             to assign any (whether it is a ghost or local) unaggregated
//             vertex that it has to an aggregate that it owns. In other words,
//             processor p attempts to assign vertex v to aggregate y where
//             y is owned by p but v may be a ghost vertex (and so really 
//             assigned to another processor). Deciding which aggregate
//             a vertex is assigned to is done by scoring. Roughly, we want 
//
//                  a) larger scores for y if v is is close (graph distance)
//                     to y's root.
//                  b) larger scores for y if v has direct connections to 
//                     several different vertices already assigned to y.
//                  c) lower scores for y if several vertices have already
//                     been swept into y during this phase.
//
//             Some care must be taken for vertices that are shared (either
//             local vertices that are sent to other processors or ghost
//             vertices) in that each processor sharing the vertex
//             will attempt to assign it to different aggregates. 
//             MueLu_ArbitrateAndCommunicate() is again used for arbitration
//             with the score being given as the weight. 
//
//             The main tricky thing occurs when v is tentatively added to y.
//             When assigning vprime to y, the assumed connection with v should
//             not be the sole basis of this decisioin if there is some chance
//             that v might be lost in arbitration. This could actually lead to
//             vprime being disconnected from the rest of the aggregate.  This
//             is by building a list of shared ids and checking that there is
//             at least one vertex in the scoring that 
//             is either not shared or has already been definitively
//             assigned to this processor's aggregate (i.e. have been assigned
//             to a local aggregate and have been through arbitration).
//             
//             Scoring is done by first giving a mark to vertices that have been
//             already been assigned to aggregates. This mark essentially
//             reflects the distance of this point to the root. Specifically,
//
//               mark(v) <-- MUELU_DISTONE_VERTEX_WEIGHT if v assigned to 
//                                                        aggregate prior
//                                                        to this phase.
//
//               mark(v) <-- max(mark(vk))/2              otherwise
//
//             where max(mark(vk)) considers all vertices definitively
//             assigned to y that have direct connections to v.
//
//             Finally,
//               score(vtilde,y)<--sum(mark(vkhat)) - AggregateIncrementPenalty
//
//             where vtilde is an unaggregated vertex being considered for
//             assignment in aggregate y and vkhat are all vertices in y
//             with a direct connection to vtilde. AggregateIncrementPenalty
//             is equal to 
//                 max (INCR_SCALING*NNewVtx,
//                      sum(mark(vkhat))*(1-MUELU_PENALTYFACTOR))
//             where NNewVtx is the number of phase 5 vertices already
//             assigned to y.
//
//             One last wrinkle, is that we have wrapped the whole 
//             scoring/assigning of vertices around a big loop that
//             looks something like
//                for ( Threshold = big; Threshold >= 0; Reduce(Threshold)){
//                         .
//                         .
//                         .
//                     MueLu_ArbitrateAndCommunicate() i
//                }
//
//             New vertices are swept into aggregates only if their best
//             score is >= a Threshold.  This encourages only very good
//             vertices to be assigned first followed by somewhat less
//             well connected ones in later iterations of the loop.
//             It also helps minimize the number of exclusions that would
//             occur to address the issue mentioned above where we don't want
//             to make assignment decisions based on connections to vertices
//             that might be later lost in arbitration.
//   Phase 6:  Aggregate remaining vertices and avoid small aggregates (e.g.,
//             singletons) at all costs. Typically, most everything should
//             be aggregated by Phase's 1-5.  One way that we could still have 
//             unaggegated vertices is if processor p was never assigned a
//             root node (perhaps because the number of local vertices on p
//             is less than minNodesPerAggregate) and additionally p has
//             local ids which are not shared with any other processors (so
//             that no other processor's aggregate can claim these vertices).
//             
//             Phase 6 looks at the first unassigned vertex and all of its
//             local unassigned neighbors and makes a new aggregate. If this
//             aggregate has at least minNodesPerAggregate vertices, 
//             we continue this process of creating new aggregates by 
//             examining other unassigned vertices. If the new aggregate
//             is too small, we try add the next unassigned vertex
//             and its neighbors to the same newly created aggregate. 
//             Once again, we check the size of this new aggregate to
//             decide whether other unassigned vertices should be added
//             to this aggregate or used to create a new aggregate. 
//             If the last newly created aggregate (of course there may be just
//             one newly created aggregate) is too small, we then see if
//             there is at least one other aggregate that this processor owns.
//             If so, we merge the small newly created aggregate with aggregate
//             0. If not, we accept the fact that a small aggregate has been
//             created.

//  
// One final note about the use of MueLu_ArbitrateAndCommunicate(). No
// arbitration occurs (which means the procWinner[] is not set as well) for a
// global shared id if and only if all weights on all processors corresponding
// to this id is zero. Thus, the general idea is that any global id where we
// want arbitration to occur should have at least one associated weight on 
// one processor which is nonzero. Any global id where we don't want any
// arbitration should have all weights set to 0.
//
// Note: procWinners is also set to MyPid() by MueLu_ArbitrateAndCommunicate()
// for any nonshared gid's with a nonzero weight.
//

int MueLu_AggregateLeftOvers(const AggregationOptions &aggOptions, 
                             Aggregates<int,int> &aggregates,
                             const std::string & label,
                             const Graph<int,int> &graph)
{
  int nVertices    = graph.GetNodeNumVertices();
  int exp_nRows    = nVertices + graph.GetNodeNumGhost();
  int myPid        = graph.GetComm()->getRank();
  double printFlag = aggOptions.GetPrintFlag();
  int nAggregates  = aggregates.GetNumAggregates();

  int minNodesPerAggregate = aggOptions.GetMinNodesPerAggregate();

  const RCP<const Cthulhu::Map<int,int> > nonUniqueMap = aggregates.GetMap();
  const RCP<const Cthulhu::Map<int,int> > uniqueMap = graph.GetDomainMap();

  AggAlgorithm2Comm myWidget(uniqueMap, nonUniqueMap);

  // Pull stuff out of vectors

  Teuchos::ArrayRCP<int>  vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<int>   procWinner  = aggregates.GetProcWinner()->getDataNonConst(0);

  RCP<Cthulhu::Vector<double> > weights_ = Cthulhu::VectorFactory<double>::Build(nonUniqueMap); //TODO
  Teuchos::ArrayRCP<double> weights = weights_->getDataNonConst(0); //TODO

  // Aggregated vertices not "definitively" assigned to processors are
  // arbitrated by MueLu_ArbitrateAndCommunicate(). There is some
  // additional logic to prevent losing root nodes in arbitration.
  
  weights_->putScalar(0.);
  for (size_t i=0;i<nonUniqueMap->getNodeNumElements();i++) {
    if (procWinner[i] == MUELU_UNASSIGNED) {
      if (vertex2AggId[i] != MUELU_UNAGGREGATED) {
        weights[i] = 1.;
        if (aggregates.IsRoot(i)) weights[i] = 2.;
      }
    }
  }

  myWidget.MueLu_ArbitrateAndCommunicate(*weights_, aggregates, true);
  weights_->putScalar(0.); //All tentatively assigned vertices are now definitive

  // Tentatively assign any vertex (ghost or local) which neighbors a root
  // to the aggregate associated with the root.

  for (int i = 0; i < nVertices; i++) { 
    if ( aggregates.IsRoot(i) && (procWinner[i] == myPid) ) {

      // neighOfINode is the neighbor node list of node 'i'.
      ArrayView<const int> neighOfINode = graph.getNeighborVertices(i);
      
      for (iter it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
        int colj = *it;
        if (vertex2AggId[colj] == MUELU_UNAGGREGATED) {
          weights[colj]= 1.;
          vertex2AggId[colj] = vertex2AggId[i];
        }
      }
    }
  }

  myWidget.MueLu_ArbitrateAndCommunicate(*weights_, aggregates, true);
  weights_->putScalar(0.); // All tentatively assigned vertices are now definitive

  // Record the number of aggregated vertices
  int total_phase_one_aggregated = 0;
  {
    int phase_one_aggregated = 0;
    for (int i = 0; i < nVertices; i++) {
      if (vertex2AggId[i] != MUELU_UNAGGREGATED)
        phase_one_aggregated++;
    }

    sumAll(graph.GetComm(), phase_one_aggregated, total_phase_one_aggregated);
  
  }
  
  int total_nVertices = 0;
  sumAll(graph.GetComm(), nVertices, total_nVertices);

  /* Among unaggregated points, see if we can make a reasonable size    */
  /* aggregate out of it. We do this by looking at neighbors and seeing */
  /* how many are unaggregated and on my processor. Loosely,            */
  /* base the number of new aggregates created on the percentage of     */
  /* unaggregated nodes.                                                */

  double factor = 1.;
  factor = ((double) total_phase_one_aggregated)/((double)(total_nVertices + 1));
  factor = pow(factor, aggOptions.GetPhase3AggCreation());

  for (int i = 0; i < nVertices; i++) {
    if (vertex2AggId[i] == MUELU_UNAGGREGATED) 
      {

        // neighOfINode is the neighbor node list of node 'iNode'.
        ArrayView<const int> neighOfINode = graph.getNeighborVertices(i);
        int rowi_N = neighOfINode.size();

        int nonaggd_neighbors = 0;
        for (iter it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
          int colj = *it;
          if (vertex2AggId[colj] == MUELU_UNAGGREGATED && colj < nVertices)
            nonaggd_neighbors++;
        }
        if (  (nonaggd_neighbors > minNodesPerAggregate) &&
              (((double) nonaggd_neighbors)/((double) rowi_N) > factor))
          {
            vertex2AggId[i] = (nAggregates)++;
            for (iter it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
              int colj = *it;
              if (vertex2AggId[colj]==MUELU_UNAGGREGATED) {
                vertex2AggId[colj] = vertex2AggId[i];
                if (colj < nVertices) weights[colj] = 2.;
                else                  weights[colj] = 1.;
              }
            }
            aggregates.SetIsRoot(i);
            weights[i] = 2.;
          }
      }
  } /*for (i = 0; i < nVertices; i++)*/
  
  myWidget.MueLu_ArbitrateAndCommunicate(*weights_, aggregates, true);
  weights_->putScalar(0.);//All tentatively assigned vertices are now definitive


  if ( printFlag < MueLu_PrintLevel()) {
    int Nphase1_agg = nAggregates;
    int total_aggs;

    sumAll(graph.GetComm(), Nphase1_agg, total_aggs);

    if (myPid == 0) {
      printf("Aggregation(%s) : Phase 1 - nodes aggregated = %d \n",label.c_str(),
             total_phase_one_aggregated);
      printf("Aggregation(%s) : Phase 1 - total aggregates = %d\n",label.c_str(), total_aggs);
    }
    int i = nAggregates - Nphase1_agg;

    { int ii; sumAll(graph.GetComm(),i,ii); i = ii; }
    if ( myPid == 0 ) {
      printf("Aggregation(%s) : Phase 3 - additional aggregates = %d\n",label.c_str(), i);
    }
  }

  // Determine vertices that are not shared by setting Temp to all ones
  // and doing MueLu_NonUnique2NonUnique(..., ADD). This sums values of all
  // local copies associated with each Gid. Thus, sums > 1 are shared.

  RCP<Cthulhu::Vector<double> > temp_ = Cthulhu::VectorFactory<double>::Build(nonUniqueMap);

  RCP<Cthulhu::Vector<double> > tempOutput_ = Cthulhu::VectorFactory<double>::Build(nonUniqueMap);
  Teuchos::ArrayRCP<double> tempOutput = tempOutput_->getDataNonConst(0);

  temp_->putScalar(1.);  
  tempOutput_->putScalar(0.); 

  myWidget.MueLu_NonUnique2NonUnique(*temp_, *tempOutput_, Cthulhu::ADD);
   
  vector<bool> gidNotShared(exp_nRows);
  for (int i = 0; i < exp_nRows; i++) {
    if (tempOutput[i] > 1.) gidNotShared[i] = false; 
    else  gidNotShared[i] = true; 
  }

  // Phase 4. 

  double nAggregatesTarget;
  int    nAggregatesGlobal, minNAggs, maxNAggs;

  nAggregatesTarget = ((double)  uniqueMap->getGlobalNumElements())*
    (((double) uniqueMap->getGlobalNumElements())/
     ((double) graph.GetGlobalNumEdges()));

  sumAll(graph.GetComm(), nAggregates, nAggregatesGlobal);

  minAll(graph.GetComm(), nAggregates, minNAggs);
  maxAll(graph.GetComm(), nAggregates, maxNAggs);

  //
  // Only do this phase if things look really bad. THIS
  // CODE IS PRETTY EXPERIMENTAL 
  //
#define MUELU_PHASE4BUCKETS 6
  if ((nAggregatesGlobal < graph.GetComm()->getSize()) &&
      (2.5*nAggregatesGlobal < nAggregatesTarget) &&
      (minNAggs ==0) && (maxNAggs <= 1)) {

    //TODO TODO
//     Epetra_Util util;
//     util.SetSeed( (unsigned int) myPid*2 + (int) (11*rand()));
//     k = (int)ceil( (10.*myPid)/graph.GetComm()->getSize());
//     for (i = 0; i < k+7; i++) util.SetSeed( (unsigned int) util.RandomInt() );
//    temp_->SetSeed( (unsigned int) util.RandomInt() );
    temp_->randomize(); 

    Teuchos::ArrayRCP<double> temp = temp_->getDataNonConst(0);

    // build a list of candidate root nodes (vertices not adjacent
    // to aggregated vertices)

    int nCandidates = 0, nCandidatesGlobal;

    Teuchos::ArrayRCP<int> candidates = Teuchos::arcp<int>(nVertices+1);

    double priorThreshold = 0.;
    for (int kkk = 0 ; kkk < MUELU_PHASE4BUCKETS; kkk++) {
      MueLu_RootCandidates(nVertices, vertex2AggId, graph,
                           candidates, nCandidates, nCandidatesGlobal);


      double nTargetNewGuys =  nAggregatesTarget - nAggregatesGlobal;
      double threshold      =  priorThreshold + (1. - priorThreshold)*
        nTargetNewGuys/(nCandidatesGlobal + .001);
   

      threshold = (threshold*(kkk+1.))/((double) MUELU_PHASE4BUCKETS);
      priorThreshold = threshold;

      for (int k = 0; k < nCandidates ; k++ ) { 
        int i = candidates[k];                  
        if ((vertex2AggId[i] == MUELU_UNAGGREGATED) && 
            (fabs(temp[i])  < threshold)) {
          // Note: priorThreshold <= fabs(temp[i]) <= 1

          // neighOfINode is the neighbor node list of node 'iNode'.
          ArrayView<const int> neighOfINode = graph.getNeighborVertices(i);

          if (neighOfINode.size() > minNodesPerAggregate) { //TODO: check if this test is exactly was we want to do
            int count = 0;
            for (iter it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
              int Adjacent    = *it;
              // This might not be true if someone close to i
              // is chosen as a root via fabs(temp[]) < Threshold
              if (vertex2AggId[Adjacent] == MUELU_UNAGGREGATED){
                count++;
                vertex2AggId[Adjacent] = nAggregates;
                weights[Adjacent] = 1.;
              }
            }
            if (count >= minNodesPerAggregate) {
              vertex2AggId[i] = nAggregates++;
              weights[i] = 2.;
              aggregates.SetIsRoot(i);
            }
            else { // undo things
              for (iter it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
                int Adjacent    = *it;
                if (vertex2AggId[Adjacent] == nAggregates){
                  vertex2AggId[Adjacent] = MUELU_UNAGGREGATED;
                  weights[Adjacent] = 0.;
                }
              }
            }

          }
        }
      }

      myWidget.MueLu_ArbitrateAndCommunicate(*weights_, aggregates, true);
      weights_->putScalar(0.); // All tentatively assigned vertices are now definitive
      sumAll(graph.GetComm(), nAggregates, nAggregatesGlobal);

      // check that there are no aggregates sizes below minNodesPerAggregate
     
      aggregates.SetNumAggregates(nAggregates);

      MueLu_RemoveSmallAggs(aggregates, minNodesPerAggregate, 
			    weights_, myWidget);
      
      nAggregates = aggregates.GetNumAggregates();
    }   // one possibility
  }

  // Initialize things for Phase 5. This includes building the transpose
  // of the matrix ONLY for transposed rows that correspond to unaggregted
  // ghost vertices. Further, the transpose is only a local transpose. 
  // Nonzero edges which exist on other processors are not represented.

  Teuchos::ArrayRCP<int> Mark = Teuchos::arcp<int>(exp_nRows+1);
  Teuchos::ArrayRCP<int> agg_incremented = Teuchos::arcp<int>(nAggregates+1);
  Teuchos::ArrayRCP<int> SumOfMarks = Teuchos::arcp<int>(nAggregates+1);

  for (int i = 0; i < exp_nRows; i++)   Mark[i] = MUELU_DISTONE_VERTEX_WEIGHT;
  for (int i = 0; i < nAggregates; i++) agg_incremented[i] = 0;
  for (int i = 0; i < nAggregates; i++) SumOfMarks[i] = 0;

  // Grab the transpose matrix graph for unaggregated ghost vertices.
  //     a) count the number of nonzeros per row in the transpose

  vector<int> RowPtr(exp_nRows+1-nVertices);
  for (int i = nVertices; i < exp_nRows ;  i++) RowPtr[i-nVertices] = 0;
  for (int i = 0; i < nVertices;  i++) {

    // neighOfINode is the neighbor node list of node 'iNode'.
    ArrayView<const int> neighOfINode = graph.getNeighborVertices(i);

    for (iter it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
      int j = *it;
      if ( (j >= nVertices) && (vertex2AggId[j] == MUELU_UNAGGREGATED)){
        RowPtr[j-nVertices]++;
      }
    }
  }

  //     b) Convert RowPtr[i] to point to 1st first nnz spot in row i.

  int iSum = 0, iTemp;
  for (int i = nVertices; i < exp_nRows ;  i++) {
    iTemp = RowPtr[i-nVertices];
    RowPtr[i-nVertices] = iSum;
    iSum += iTemp;
  }
  RowPtr[exp_nRows-nVertices] = iSum;
  vector<int> cols(iSum+1);
   
  //     c) Traverse matrix and insert entries in proper location.
  for (int i = 0; i < nVertices;  i++) {

    // neighOfINode is the neighbor node list of node 'iNode'.
    ArrayView<const int> neighOfINode = graph.getNeighborVertices(i);
    
    for (iter it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
      int j = *it;
      if ( (j >= nVertices) && (vertex2AggId[j] == MUELU_UNAGGREGATED)){
        cols[RowPtr[j-nVertices]++] = i;
      }
    }
  }

  //     d) RowPtr[i] points to beginning of row i+1 so shift by one location.
  for (int i = exp_nRows; i > nVertices;  i--)
    RowPtr[i-nVertices] = RowPtr[i-1-nVertices];
  RowPtr[0] = 0;
  
  int bestScoreCutoff;
  int thresholds[10] = {300,200,100,50,25,13,7,4,2,0};

  // Stick unaggregated vertices into existing aggregates as described above. 
 
  bool cannotLoseAllFriends; // Used to address possible loss of vertices in 
  // arbitration of shared nodes discussed above.
  for (int kk = 0; kk < 10; kk += 2) {
    bestScoreCutoff = thresholds[kk];
    for (int i = 0; i < exp_nRows; i++) {
      if (vertex2AggId[i] == MUELU_UNAGGREGATED) {

        // neighOfINode is the neighbor node list of node 'iNode'.
        ArrayView<const int> neighOfINode;

        // Grab neighboring vertices which is either in graph for local ids
        // or sits in transposed fragment just constructed above for ghosts.
        if (i < nVertices) {
          neighOfINode = graph.getNeighborVertices(i);
        }
        else {
          int *rowi_col = NULL, rowi_N;
          rowi_col = &(cols[RowPtr[i-nVertices]]);
          rowi_N   = RowPtr[i+1-nVertices] - RowPtr[i-nVertices];
          
          neighOfINode = ArrayView<const int>(rowi_col, rowi_N);
        }
        for (iter it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
          int Adjacent    = *it;
          int AdjacentAgg = vertex2AggId[Adjacent];

          //Adjacent is aggregated and either I own the aggregate
          // or I could own the aggregate after arbitration.
          if ((AdjacentAgg != MUELU_UNAGGREGATED) && 
              ((procWinner[Adjacent] == myPid) ||     
               (procWinner[Adjacent] == MUELU_UNASSIGNED))){
            SumOfMarks[AdjacentAgg] += Mark[Adjacent];
          }
        }
        int best_score = MUELU_NOSCORE;
        int best_agg; // TODO: init ?
        int BestMark;
        for (iter it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
          int Adjacent    = *it;
          int AdjacentAgg = vertex2AggId[Adjacent];
          //Adjacent is unaggregated, has some value and no
          //other processor has definitively claimed him
          if ((AdjacentAgg != MUELU_UNAGGREGATED) && 
              (SumOfMarks[AdjacentAgg] != 0) &&
              ((procWinner[Adjacent] == myPid) ||
               (procWinner[Adjacent] == MUELU_UNASSIGNED ))) {

            // first figure out the penalty associated with
            // AdjacentAgg having already been incremented 
            // during this phase, then compute score.

            double penalty = (double) (INCR_SCALING*agg_incremented[AdjacentAgg]);
            if (penalty > MUELU_PENALTYFACTOR*((double)SumOfMarks[AdjacentAgg]))
              penalty = MUELU_PENALTYFACTOR*((double)SumOfMarks[AdjacentAgg]);
            int score = SumOfMarks[AdjacentAgg]- ((int) floor(penalty));

            if (score > best_score) { 
              best_agg             = AdjacentAgg; 
              best_score           = score;
              BestMark             = Mark[Adjacent];
              cannotLoseAllFriends = false;
 
              // This address issue mentioned above by checking whether
              // Adjacent could be lost in arbitration. weight==0 means that
              // Adjacent was not set during this loop of Phase 5 (and so it
              // has already undergone arbitration). GidNotShared == true 
              // obviously implies that Adjacent cannot be lost to arbitration
              if ((weights[Adjacent]== 0.) || (gidNotShared[Adjacent] == true))
                cannotLoseAllFriends = true;
            }
            // Another vertex within current best aggregate found.
            // We should have (best_score == score). We need to see
            // if we can improve BestMark and cannotLoseAllFriends.
            else if (best_agg == AdjacentAgg) {
              if ((weights[Adjacent]== 0.) || (gidNotShared[Adjacent] == true))
                cannotLoseAllFriends = true;
              if (Mark[Adjacent] > BestMark) BestMark = Mark[Adjacent];
            }
          }
        }
        // Clean up
         for (iter it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
          int Adjacent    = *it;
          int AdjacentAgg = vertex2AggId[Adjacent];
          if (AdjacentAgg >= 0) SumOfMarks[AdjacentAgg] = 0;
        }
        // Tentatively assign vertex to best_agg. 
        if ( (best_score >= bestScoreCutoff) && (cannotLoseAllFriends)) { 
          vertex2AggId[i] = best_agg;
          weights[i] = best_score;
          agg_incremented[best_agg]++;
          Mark[i] = (int) ceil(   ((double) BestMark)/2.);
        }
      }
    }

    myWidget.MueLu_ArbitrateAndCommunicate(*weights_, aggregates, true);
    weights_->putScalar(0.); // All tentatively assigned vertices are now definitive
  }

  // Phase 6: Aggregate remain unaggregated vertices and try at all costs
  //          to avoid small aggregates.
  //          One case where we can find ourselves in this situation
  //          is if all vertices vk adjacent to v have already been
  //          put in other processor's aggregates and v does not have
  //          a direct connection to a local vertex in any of these
  //          aggregates.

  int count = 0;
  int Nleftover = 0, Nsingle = 0;
  for (int i = 0; i < nVertices; i++) { 
    if ((vertex2AggId[i] == MUELU_UNAGGREGATED) ) {
      Nleftover++;

      // neighOfINode is the neighbor node list of node 'iNode'.
      ArrayView<const int> neighOfINode = graph.getNeighborVertices(i);

      // We don't want too small of an aggregate. So lets see if there is an
      // unaggregated neighbor that we can also put with this vertex

      vertex2AggId[i] = nAggregates;
      weights[i] = 1.;
      if (count == 0) aggregates.SetIsRoot(i);
      count++;
      for (iter it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
        int j = *it;
        if ((j != i)&&(vertex2AggId[j] == MUELU_UNAGGREGATED)&&
            (j < nVertices)) {
          vertex2AggId[j] = nAggregates;
          weights[j] = 1.;
          count++;
        }
      }
      if ( count >= minNodesPerAggregate) {
        nAggregates++; 
        count = 0;
      }
    }
  }
  // We have something which is under minNodesPerAggregate when 
  if (count != 0) {
    // Can stick small aggregate with 0th aggregate?
    if (nAggregates > 0) {
      for (int i = 0; i < nVertices; i++) {
        if ((vertex2AggId[i] == nAggregates) && (procWinner[i] == myPid)){
          vertex2AggId[i] = 0;
          aggregates.SetIsRoot(i,false);
        }
      }
    }
    else {
      Nsingle++;
      nAggregates++;
    }
  }

  myWidget.MueLu_ArbitrateAndCommunicate(*weights_, aggregates, false);


  if (printFlag < MueLu_PrintLevel()) {
    { int total_Nsingle=0;   sumAll(graph.GetComm(), Nsingle, total_Nsingle);     Nsingle = total_Nsingle; }
    { int total_Nleftover=0; sumAll(graph.GetComm(), Nleftover, total_Nleftover); Nleftover = total_Nleftover; }
    int total_aggs;
    sumAll(graph.GetComm(), nAggregates, total_aggs);
    if ( myPid == 0 ) {
      printf("Aggregation(%s) : Phase 3 - total aggregates = %d\n",label.c_str(), total_aggs);
      printf("Aggregation(%s) : Phase 6 - leftovers = %d and singletons = %d\n",label.c_str() ,Nleftover, Nsingle);
    }
  }

  aggregates.SetNumAggregates(nAggregates);

  return 0;
}

// build a list of candidate root nodes (vertices not adjacent to already
// aggregated vertices)
int MueLu_RootCandidates(int nVertices, Teuchos::ArrayRCP<int> vertex2AggId, const Graph<int,int> graph,
                         Teuchos::ArrayRCP<int> &candidates, int &nCandidates, int &nCandidatesGlobal) {

  //TODO: interface : ArrayRCP useless ?

  nCandidates = 0;
 
  for (int i = 0; i < nVertices; i++ ) {
    if (vertex2AggId[i] == MUELU_UNAGGREGATED) {
      bool noAggdNeighbors = true;

      // neighOfINode is the neighbor node list of node 'iNode'.
      ArrayView<const int> neighOfINode = graph.getNeighborVertices(i);

      for (iter it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
        int adjacent    = *it;
        if (vertex2AggId[adjacent] != MUELU_UNAGGREGATED) 
          noAggdNeighbors = false;
      }
      if (noAggdNeighbors == true) candidates[nCandidates++] = i;
    }
  }

  sumAll(graph.GetComm(), nCandidates, nCandidatesGlobal);

  return 0;
}

// Compute sizes of all the aggregates.
int MueLu_ComputeAggSizes(Aggregates<int,int> & aggregates, Teuchos::ArrayRCP<int> & aggSizes)
{
  int myPid = aggregates.GetMap()->getComm()->getRank();

  int nAggregates = aggregates.GetNumAggregates();

  Teuchos::ArrayRCP<int> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);
  Teuchos::ArrayRCP<int> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  int size = procWinner.size();

  for (int i = 0; i < nAggregates; i++) aggSizes[i] = 0;
  for (int k = 0; k < size; k++ ) {
    if (procWinner[k] == myPid) aggSizes[vertex2AggId[k]]++;
  }

  return 0; //TODO
}

int MueLu_RemoveSmallAggs(Aggregates<int,int> & aggregates, int min_size,
                          RCP<Cthulhu::Vector<double> > & weights_, const AggAlgorithm2Comm & myWidget)
{
  int myPid = aggregates.GetMap()->getComm()->getRank();
  
  int nAggregates = aggregates.GetNumAggregates();

  Teuchos::ArrayRCP<int> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);
  Teuchos::ArrayRCP<int> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  int size = procWinner.size();

  Teuchos::ArrayRCP<int> AggInfo = Teuchos::arcp<int>(nAggregates+1);

  Teuchos::ArrayRCP<double> weights = weights_->getDataNonConst(0);

  MueLu_ComputeAggSizes(aggregates, AggInfo);

  // Make a list of all aggregates indicating New AggId
  // Use AggInfo array for this.

  int NewNAggs = 0; 
  for (int i = 0; i < nAggregates; i++) {
    if ( AggInfo[i] < min_size) { 
      AggInfo[i] =  MUELU_UNAGGREGATED;
    }
    else AggInfo[i] = NewNAggs++;
  }

  for (int k = 0; k < size; k++ ) {
    if (procWinner[k] == myPid) {
      if (vertex2AggId[k] !=  MUELU_UNAGGREGATED) {
        vertex2AggId[k] = AggInfo[vertex2AggId[k]];
        weights[k] = 1.;
      }
      if (vertex2AggId[k] ==  MUELU_UNAGGREGATED) 
        aggregates.SetIsRoot(k,false);
    }
  }
  nAggregates = NewNAggs;

  myWidget.MueLu_ArbitrateAndCommunicate(*weights_, aggregates, true);
  weights_->putScalar(0.); // All tentatively assigned vertices are now definitive

  // procWinner is not set correctly for aggregates which have 
  // been eliminated
  for (int i = 0; i < size; i++) {
    if (vertex2AggId[i] == MUELU_UNAGGREGATED) 
      procWinner[i] = MUELU_UNASSIGNED;
  }
  aggregates.SetNumAggregates(nAggregates);

  return 0; //TODO
}
