/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
// 
// ***********************************************************************
// 
// @HEADER

*/

// Teko includes
#include "Teko_SmootherPreconditionerFactory.hpp"

#include "Teko_PreconditionerInverseFactory.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

namespace Teko {

// A small utility class to help distinguish smoothers
class SmootherRequestMesg : public RequestMesg {
public:
   SmootherRequestMesg(unsigned int block)
      : RequestMesg("__smoother_request_message__")
      , block_(block) {}

   unsigned int getBlock() const 
   { return block_; }

private:
   unsigned int block_;
};

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


SmootherParentLinearOp::SmootherParentLinearOp(const LinearOp & smootherOp,Teuchos::RCP<RequestHandler> & rh)
   : smootherOp_(smootherOp)
{
   setRequestHandler(rh); 
   getRequestHandler()->addRequestCallback(Teuchos::rcp(this,false));
}

/** @brief Perform a matrix vector multiply with this implicitly
  * defined blocked operator. 
  *
  * The <code>apply</code> function takes one vector as input 
  * and applies a linear operator. The result
  * is returned in \f$y\f$. If this operator is reprsented as \f$M\f$ then
  * \f$ y = \alpha M x + \beta y \f$
  *
  * @param[in]     x 
  * @param[in,out] y 
  * @param[in]     alpha (default=1)
  * @param[in]     beta  (default=0)
  */
void SmootherParentLinearOp::implicitApply(const MultiVector & x, MultiVector & y,
                                           const double alpha, const double beta) const
{
   initialGuess_ = toBlockedMultiVector(Teko::deepcopy(y));

   Teko::applyOp(smootherOp_,x,y,alpha,beta);
}

//! Set the request handler with pointers to the appropriate callbacks
void SmootherParentLinearOp::setRequestHandler(const Teuchos::RCP<RequestHandler> & rh)
{ 
   Teko_DEBUG_SCOPE("SmootherParentLinearOp::setRequestHandler",10);
   requestHandler_ = rh;
}

//! Get the request handler with pointers to the appropriate callbacks
Teuchos::RCP<RequestHandler> SmootherParentLinearOp::getRequestHandler() const
{ 
   Teko_DEBUG_SCOPE("SmootherParentLinearOp::getRequestHandler",10);
   return requestHandler_;
}

MultiVector SmootherParentLinearOp::request(const RequestMesg & rm)
{
   TEUCHOS_ASSERT(SmootherParentLinearOp::handlesRequest(rm));

   // by above assertion this should always succeed!
   const SmootherRequestMesg & srm 
         = Teuchos::dyn_cast<const SmootherRequestMesg>(rm);

   return Teko::getBlock(srm.getBlock(),initialGuess_);
}

bool SmootherParentLinearOp::handlesRequest(const RequestMesg & rm)
{
   bool handled = false;

   try {
      // try to dynamic cast to a smoother request type
      Teuchos::dyn_cast<const SmootherRequestMesg>(rm);

      // dynamic cast succeeded, mesg will be handled
      handled = true;
   }
   catch( const std::bad_cast &e ) {
      // failed dynamic cast: will not handle mesg
      handled = true; 
   }

   return handled;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

SmootherLinearOp::SmootherLinearOp(const LinearOp & A, const LinearOp & invM,unsigned int applications,bool useDestAsInitialGuess)
   : A_(A), invM_(invM), applications_(applications), initialGuessType_(Unspecified), requestMesg_(Teuchos::null)
{
   // set initial guess type
   initialGuessType_ = useDestAsInitialGuess ? DestAsInitialGuess : NoInitialGuess;
}

SmootherLinearOp::SmootherLinearOp(const LinearOp & A, const LinearOp & invM,unsigned int applications,unsigned int block)
   : A_(A), invM_(invM), applications_(applications), initialGuessType_(RequestInitialGuess), requestMesg_(Teuchos::null)
{
   requestMesg_ = Teuchos::rcp(new SmootherRequestMesg(block));
}

void SmootherLinearOp::implicitApply(const MultiVector & b, MultiVector & x,
                                     const double alpha, const double beta) const
{
   using Teuchos::RCP;

   MultiVector residual = deepcopy(b); // residual = b
   MultiVector scrap = deepcopy(b);    // scrap = b
   MultiVector error;                  // location for initial guess

   // construct initial guess: required to assign starting point for destination
   // vector appropriately
   switch(initialGuessType_) {
   case RequestInitialGuess:
      // get initial guess from request handler
      error = deepcopy(getRequestHandler()->request<MultiVector>(*requestMesg_));
      Thyra::assign<double>(x.ptr(),*error); // x = error (initial guess)
      break;
   case DestAsInitialGuess:
      error = deepcopy(x);    // error = x
      break;
   case NoInitialGuess:
      Thyra::assign<double>(x.ptr(),0.0); // x = 0
      error = deepcopy(x);                // error = x
      break;
   case Unspecified:
   default:
      TEUCHOS_ASSERT(false);
   }

   for(unsigned int current=0;current<applications_;++current) {
      // compute current residual
      Teko::applyOp(A_,error,scrap);
      Teko::update(-1.0,scrap,1.0,residual); // residual = residual-A*error

      // compute appoximate correction using residual
      Thyra::assign(error.ptr(),0.0); // set error to zero
      Teko::applyOp(invM_,residual,error);

      // update solution with error
      Teko::update(1.0,error,1.0,x); // x = x+error
   }
}

//! Set the request handler with pointers to the appropriate callbacks
void SmootherLinearOp::setRequestHandler(const Teuchos::RCP<RequestHandler> & rh)
{ 
   Teko_DEBUG_SCOPE("SmootherLinearOp::setRequestHandler",10);
   requestHandler_ = rh;
}

//! Get the request handler with pointers to the appropriate callbacks
Teuchos::RCP<RequestHandler> SmootherLinearOp::getRequestHandler() const
{ 
   Teko_DEBUG_SCOPE("SmootherLinearOp::getRequestHandler",10);
   return requestHandler_;
}

LinearOp buildSmootherLinearOp(const LinearOp & A,const LinearOp & invM,unsigned int applications,bool useDestAsInitialGuess)
{
   return Teuchos::rcp(new SmootherLinearOp(A,invM,applications,useDestAsInitialGuess));
}

LinearOp buildSmootherLinearOp(const LinearOp & A,const LinearOp & invM,unsigned int applications,unsigned int initialGuessBlock)
{
   return Teuchos::rcp(new SmootherLinearOp(A,invM,applications,initialGuessBlock));
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

//! Default constructor, for use with the AutoClone class.
SmootherPreconditionerFactory::SmootherPreconditionerFactory()
   : sweepCount_(0), initialGuessType_(Unspecified), initialGuessBlock_(0), precFactory_(Teuchos::null) 
{ }

/** \brief Function that is called to build the preconditioner
  *        for the linear operator that is passed in.
  */
LinearOp SmootherPreconditionerFactory::buildPreconditionerOperator(LinearOp & lo,PreconditionerState & state) const
{
   TEST_FOR_EXCEPTION(precFactory_==Teuchos::null,std::runtime_error,
                      "ERROR: Teko::SmootherPreconditionerFactory::buildPreconditionerOperator requires that a "
                   << "preconditioner factory has been set. Currently it is null!");

   // preconditions
   TEUCHOS_ASSERT(sweepCount_>0);
   TEUCHOS_ASSERT(initialGuessType_!=Unspecified);
   TEUCHOS_ASSERT(precFactory_!=Teuchos::null);

   // build user specified preconditioner
   ModifiableLinearOp & invM = state.getModifiableOp("prec");
   if(invM==Teuchos::null)
      invM = Teko::buildInverse(*precFactory_,lo,state);
   else
      Teko::rebuildInverse(*precFactory_,lo,invM);

   // conditional on initial guess type, build the smoother
   switch(initialGuessType_) {
   case RequestInitialGuess:
      return buildSmootherLinearOp(lo,invM,sweepCount_,initialGuessBlock_);
   case DestAsInitialGuess:
      return buildSmootherLinearOp(lo,invM,sweepCount_,true); // use an initial guess
   case NoInitialGuess:
      return buildSmootherLinearOp(lo,invM,sweepCount_,false); // no initial guess
   case Unspecified:
   default:
      TEUCHOS_ASSERT(false);
   }

   // should never get here
   TEUCHOS_ASSERT(false);
   return Teuchos::null;
}

/** \brief This function builds the internals of the preconditioner factory
  *        from a parameter list.
  */
void SmootherPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList & settings)
{
   // declare all strings used by this initialization routine
   //////////////////////////////////////

   const std::string str_sweepCount = "Sweep Count";
   const std::string str_initialGuessBlock = "Initial Guess Block";
   const std::string str_destAsInitialGuess = "Destination As Initial Guess";
   const std::string str_precType = "Preconditioner Type";

   // default parameters
   //////////////////////////////////////

   initialGuessType_ = Unspecified;
   initialGuessBlock_ = 0;
   sweepCount_ = 0;
   precFactory_ = Teuchos::null;

   // get sweep counts
   //////////////////////////////////////

   if(settings.isParameter(str_sweepCount))
      sweepCount_ = settings.get<int>(str_sweepCount);

   // get initial guess
   //////////////////////////////////////

   if(settings.isParameter(str_initialGuessBlock)) {
      initialGuessBlock_ = settings.get<int>(str_initialGuessBlock);
      initialGuessType_ = RequestInitialGuess;
   }

   if(settings.isParameter(str_destAsInitialGuess)) {
      bool useDest = settings.get<bool>(str_destAsInitialGuess);
      if(useDest) {
         TEST_FOR_EXCEPTION(initialGuessType_!=Unspecified, std::runtime_error,
                            "Cannot set both \"" << str_initialGuessBlock  <<  
                            "\" and \""          << str_destAsInitialGuess << "\"");

         initialGuessType_ = DestAsInitialGuess;
      }
   }

   // safe to assume if the other values are not turned on there is no initial guess
   if(initialGuessType_==Unspecified) 
      initialGuessType_ = NoInitialGuess;

   // get preconditioner factory
   //////////////////////////////////////
 
   TEST_FOR_EXCEPTION(not settings.isParameter(str_precType),std::runtime_error,
                      "Parameter \"" << str_precType << "\" is required by a Teko::SmootherPreconditionerFactory");
      
   // grab library and preconditioner name
   Teuchos::RCP<const InverseLibrary> il = getInverseLibrary();
   std::string precName = settings.get<std::string>(str_precType);

   // build preconditioner factory
   precFactory_ = il->getInverseFactory(precName);
   TEST_FOR_EXCEPTION(precFactory_==Teuchos::null,std::runtime_error,
                      "ERROR: \"" << str_precType << "\" = " << precName 
                   << " could not be found");

   // post conditions required to be satisfied
   //////////////////////////////////////

   TEUCHOS_ASSERT(sweepCount_>0);
   TEUCHOS_ASSERT(initialGuessType_!=Unspecified);
   TEUCHOS_ASSERT(precFactory_!=Teuchos::null);
}

} // end namespace Teko
