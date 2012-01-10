// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_OrderingProblem.hpp

  This file contains the OrderingProblem class, which derives from 
  the Problem class.
*/

#ifndef _ZOLTAN2_ORDERINGPROBLEM_HPP_
#define _ZOLTAN2_ORDERINGPROBLEM_HPP_

#include <Zoltan2_Problem.hpp>
//#include <Zoltan2_OrderingAlgorithms.hpp> // TODO: Fix include path?
#include "algorithms/order/Zoltan2_OrderingAlgorithms.hpp"
#include <Zoltan2_OrderingSolution.hpp>

#include <Zoltan2_GraphModel.hpp>
#include <string>
#ifdef HAVE_ZOLTAN2_OVIS
#include <ovis.h>
#endif

using Teuchos::rcp_dynamic_cast;
using namespace std;

namespace Zoltan2{

////////////////////////////////////////////////////////////////////////
template<typename Adapter>
class OrderingProblem : public Problem<Adapter>
{
public:

  typedef typename Adapter::gid_t gid_t;
  typedef typename Adapter::lno_t lno_t;

  // Destructor
  virtual ~OrderingProblem() {};

  //! Constructor with InputAdapter Interface
#ifdef HAVE_ZOLTAN2_MPI
  OrderingProblem(Adapter *A, ParameterList *p, MPI_Comm comm=MPI_COMM_WORLD) 
                      : Problem<Adapter>(A, p, comm) 
#else
  OrderingProblem(Adapter *A, ParameterList *p) : Problem<Adapter>(A, p) 
#endif
  {
    HELLO;
    createOrderingProblem();
  };

  // Other methods
  void solve();
  // virtual void redistribute();

  OrderingSolution<gid_t, lno_t> *getSolution() {
    return solution_.getRawPtr();
  };

private:
  void createOrderingProblem();

  RCP<OrderingSolution<gid_t, lno_t> > solution_;

};

////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void OrderingProblem<Adapter>::solve()
{
  HELLO;

  size_t nVtx = this->graphModel_->getLocalNumVertices();

  // TODO: Assuming one MPI process now. nVtx = ngids = nlids
  this->solution_ = rcp(new OrderingSolution<gid_t, lno_t>(nVtx, nVtx));

  // Determine which algorithm to use based on defaults and parameters.
  // For now, assuming RCM.
  // Need some exception handling here, too.

  string method = this->params_->template get<string>("ORDER_METHOD", "RCM");
  typedef typename Adapter::base_adapter_t base_adapter_t;

  // TODO: Ignore case
  if (method.compare("RCM") == 0)
  {
      AlgRCM<base_adapter_t>(this->graphModel_, this->solution_, this->params_,
                      this->comm_);
  }
  else if (method.compare("random") == 0)
  {
      AlgRandom<base_adapter_t>(this->identifierModel_, this->solution_, this->params_,
                      this->comm_);
  }
  else if (method.compare("Minimum_Degree") == 0)
  {
      string pkg = this->params_->template get<string>("ORDER_PACKAGE", "AMD");
      if (pkg.compare("AMD") == 0)
      {
          AlgAMD<base_adapter_t>(this->graphModel_, this->solution_, this->params_,
                          this->comm_);
      }
  }
}

////////////////////////////////////////////////////////////////////////
//template <typename Adapter>
//void OrderingProblem<Adapter>::redistribute()
//{
//  HELLO;
//}

////////////////////////////////////////////////////////////////////////
//! createOrderingProblem 
//  Method with common functionality for creating a OrderingProblem.
//  Individual constructors do appropriate conversions of input, etc.
//  This method does everything that all constructors must do.

template <typename Adapter>
void OrderingProblem<Adapter>::createOrderingProblem()
{
  HELLO;
//  cout << __func__ << " input adapter type " 
//       << this->inputAdapter_->inputAdapterType() << " " 
//       << this->inputAdapter_->inputAdapterName() << endl;

#ifdef HAVE_ZOLTAN2_OVIS
  ovis_enabled(this->comm_->getRank());
#endif

  // Finalize parameters.  If the Problem wants to set or
  // change any parameters, do it before this call.

  this->env_->commitParameters();

  // Determine which parameters are relevant here.
  // For now, assume parameters similar to Zoltan:
  //   MODEL = graph, hypergraph, geometric, ids
  //   ALGORITHM = rcm, random, amd

  ModelType modelType = GraphModelType;

  typedef typename Adapter::base_adapter_t base_adapter_t;

  // Select Model based on parameters and InputAdapter type
  switch (modelType) {

  case GraphModelType:
    this->graphModel_ = rcp(new GraphModel<base_adapter_t>(
      this->baseInputAdapter_, this->envConst_, this->comm_, false, true));

    this->generalModel_ = rcp_implicit_cast<const Model<base_adapter_t> >(
      this->graphModel_);

    break;



  case IdentifierModelType:
    this->identifierModel_ = rcp(new IdentifierModel<base_adapter_t>(
      this->baseInputAdapter_, this->envConst_, this->comm_, false));

    this->generalModel_ = rcp_implicit_cast<const Model<base_adapter_t> >(
      this->identifierModel_);

    break;

  case HypergraphModelType:
  case GeometryModelType:
    cout << __func__ << " Model type " << modelType << " not yet supported." 
         << endl;
    break;

  default:
    cout << __func__ << " Invalid model" << modelType << endl;
    break;
  }
}

}
#endif
