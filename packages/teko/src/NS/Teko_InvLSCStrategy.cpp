#include "NS/Teko_InvLSCStrategy.hpp"

#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_EpetraLinearOp.hpp"

#include "Epetra_Vector.h"
#include "Epetra_Map.h"

#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_VectorOut.h"

#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"

// Teko includes
#include "Teko_Utilities.hpp"
#include "NS/Teko_LSCPreconditionerFactory.hpp"
#include "Epetra/Teko_EpetraHelpers.hpp"
#include "Epetra/Teko_EpetraOperatorWrapper.hpp"

using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcp_const_cast;

namespace Teko {
namespace NS {

// helper function for a lid driven cavity-like problem
// This function _WILL_ change the operator
Teko::ModifiableLinearOp reduceCrsOperator(Teko::ModifiableLinearOp & op,const std::vector<int> & zeroIndicies)
{
   // Extract a non-const version of the operator
   RCP<Epetra_Operator> eOp = get_Epetra_Operator(*rcp_const_cast<Thyra::LinearOpBase<double> >(op));
   RCP<Epetra_CrsMatrix> eCrsOp = rcp_dynamic_cast<Epetra_CrsMatrix>(eOp);

   // build vector for scaling
   Epetra_Vector scaling(eCrsOp->OperatorRangeMap());
   scaling.PutScalar(1.0);

   // figure out which zero index (if any) this processor owns
   std::vector<std::pair<int,int> > localIndicies; // local, global indicies
   for(unsigned int i=0;i<zeroIndicies.size();i++) {
      // check if this zero index is owned by this processor
      int local = eCrsOp->LRID(zeroIndicies[i]);
      if(local>=0) 
         localIndicies.push_back(std::make_pair(local,zeroIndicies[i]));
   }

   // set a number of rows to zero
   for(unsigned int i=0;i<localIndicies.size();i++) 
      TEST_FOR_EXCEPT(scaling.ReplaceGlobalValue(localIndicies[i].second,0,0.0));

   // wipe out all but the desired rows and columns
   eCrsOp->RightScale(scaling);
   eCrsOp->LeftScale(scaling);

   #if 1
   // so the matrix is still invertable...set the digaonals of the
   // wiped rows to 1
   for(unsigned int i=0;i<localIndicies.size();i++) {
      double value = 1.0;
      int index = localIndicies[i].second;

      // set diagonal to one
      TEST_FOR_EXCEPT(eCrsOp->ReplaceGlobalValues(index,1,&value,&index));
   } 
   #endif

   // now wrap the crs matrix in a ZeroedEpetraOperator
   // return Thyra::epetraLinearOp(rcp(new Teko::Epetra::ZeroedOperator(zeroIndicies,eCrsOp)));
   return Thyra::nonconstEpetraLinearOp(eCrsOp);
}

/////////////////////////////////////////////////////////////////////////////
// InvLSCStrategy Implementation
/////////////////////////////////////////////////////////////////////////////

// constructors
/////////////////////////////////////////////////////////////////////////////
InvLSCStrategy::InvLSCStrategy()
   : massMatrix_(Teuchos::null), invFactoryF_(Teuchos::null), invFactoryS_(Teuchos::null), eigSolveParam_(5)
   , rowZeroingNeeded_(false), useFullLDU_(false), useMass_(false), useLumping_(false), useWScaling_(false), scaleType_(Diagonal)
   , isSymmetric_(true)
{ }

InvLSCStrategy::InvLSCStrategy(const Teuchos::RCP<InverseFactory> & factory,bool rzn)
   : massMatrix_(Teuchos::null), invFactoryF_(factory), invFactoryS_(factory), eigSolveParam_(5), rowZeroingNeeded_(rzn)
   , useFullLDU_(false), useMass_(false), useLumping_(false), useWScaling_(false), scaleType_(Diagonal)
   , isSymmetric_(true)
{ }

InvLSCStrategy::InvLSCStrategy(const Teuchos::RCP<InverseFactory> & invFactF,
                               const Teuchos::RCP<InverseFactory> & invFactS,
                               bool rzn)
   : massMatrix_(Teuchos::null), invFactoryF_(invFactF), invFactoryS_(invFactS), eigSolveParam_(5), rowZeroingNeeded_(rzn)
   , useFullLDU_(false), useMass_(false), useLumping_(false), useWScaling_(false), scaleType_(Diagonal)
   , isSymmetric_(true)
{ }

InvLSCStrategy::InvLSCStrategy(const Teuchos::RCP<InverseFactory> & factory,LinearOp & mass,bool rzn)
   : massMatrix_(mass), invFactoryF_(factory), invFactoryS_(factory), eigSolveParam_(5), rowZeroingNeeded_(rzn)
   , useFullLDU_(false), useMass_(false), useLumping_(false), useWScaling_(false), scaleType_(Diagonal)
   , isSymmetric_(true)
{ }

InvLSCStrategy::InvLSCStrategy(const Teuchos::RCP<InverseFactory> & invFactF,
                               const Teuchos::RCP<InverseFactory> & invFactS,
                               LinearOp & mass,bool rzn)
   : massMatrix_(mass), invFactoryF_(invFactF), invFactoryS_(invFactS), eigSolveParam_(5), rowZeroingNeeded_(rzn)
   , useFullLDU_(false), useMass_(false), useLumping_(false), useWScaling_(false), scaleType_(Diagonal)
   , isSymmetric_(true)
{ }

/////////////////////////////////////////////////////////////////////////////

void InvLSCStrategy::buildState(BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   Teko_DEBUG_SCOPE("InvLSCStrategy::buildState",10);

   LSCPrecondState * lscState = dynamic_cast<LSCPrecondState*>(&state);
   TEUCHOS_ASSERT(lscState!=0);

   // if neccessary save state information
   if(not lscState->isInitialized()) {
      Teko_DEBUG_EXPR(Teuchos::Time timer(""));

      // construct operators
      {
         Teko_DEBUG_SCOPE("LSC::buildState constructing operators",1);
         Teko_DEBUG_EXPR(timer.start(true));

         initializeState(A,lscState);

         Teko_DEBUG_EXPR(timer.stop());
         Teko_DEBUG_MSG("LSC::buildState BuildOpsTime = " << timer.totalElapsedTime(),1);
      }

      // Build the inverses
      {
         Teko_DEBUG_SCOPE("LSC::buildState calculating inverses",1);
         Teko_DEBUG_EXPR(timer.start(true));

         computeInverses(A,lscState);

         Teko_DEBUG_EXPR(timer.stop());
         Teko_DEBUG_MSG("LSC::buildState BuildInvTime = " << timer.totalElapsedTime(),1);
      }
   }
}

// functions inherited from LSCStrategy
LinearOp InvLSCStrategy::getInvBQBt(const BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   return state.getInverse("invBQBtmC");
}

LinearOp InvLSCStrategy::getInvBHBt(const BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   return state.getInverse("invBHBtmC");
}

LinearOp InvLSCStrategy::getInvF(const BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   return state.getInverse("invF");
}

LinearOp InvLSCStrategy::getInvAlphaD(const BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   LSCPrecondState * lscState = dynamic_cast<LSCPrecondState*>(&state);
   TEUCHOS_ASSERT(lscState!=0);
   TEUCHOS_ASSERT(lscState->isInitialized())

   return lscState->aiD_;
}

LinearOp InvLSCStrategy::getInvMass(const BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   LSCPrecondState * lscState = dynamic_cast<LSCPrecondState*>(&state);
   TEUCHOS_ASSERT(lscState!=0);
   TEUCHOS_ASSERT(lscState->isInitialized())

   return lscState->invMass_;
}

LinearOp InvLSCStrategy::getHScaling(const BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   if(hScaling_!=Teuchos::null) return hScaling_;
   return getInvMass(A,state);
}

//! Initialize the state object using this blocked linear operator
void InvLSCStrategy::initializeState(const BlockedLinearOp & A,LSCPrecondState * state) const
{
   Teko_DEBUG_SCOPE("InvLSCStrategy::initializeState",10);

   const LinearOp F  = getBlock(0,0,A);
   const LinearOp Bt = getBlock(0,1,A);
   const LinearOp B  = getBlock(1,0,A);
   const LinearOp C  = getBlock(1,1,A);

   LinearOp D = B;
   LinearOp G = isSymmetric_ ? Bt : adjoint(D);

   bool isStabilized = (not isZeroOp(C));

   // The logic follows like this
   //    if there is no mass matrix available --> build from F
   //    if there is a mass matrix and the inverse hasn't yet been built
   //       --> build from the mass matrix
   //    otherwise, there is already an invMass_ matrix that is appropriate
   //       --> use that one
   if(massMatrix_==Teuchos::null) {
      Teko_DEBUG_MSG("LSC::initializeState Build Scaling <F> type \"" 
                   << getDiagonalName(scaleType_) << "\"" ,1);
      state->invMass_ = getInvDiagonalOp(F,scaleType_);
   }
   else if(state->invMass_==Teuchos::null) {
      Teko_DEBUG_MSG("LSC::initializeState Build Scaling <mass> type \"" 
                   << getDiagonalName(scaleType_) << "\"" ,1);
      state->invMass_ = getInvDiagonalOp(massMatrix_,scaleType_);
   }
   // else "invMass_" should be set and there is no reason to rebuild it

   // compute BQBt
   state->BQBt_ = explicitMultiply(B,state->invMass_,Bt,state->BQBt_);
   Teko_DEBUG_MSG("Computed BQBt",10);

   // if there is no H-Scaling
   if(wScaling_!=Teuchos::null && hScaling_==Teuchos::null) {
      // from W vector build H operator scaling
      RCP<const Thyra::VectorBase<double> > w = wScaling_->col(0);
      RCP<const Thyra::VectorBase<double> > iQu 
            = rcp_dynamic_cast<const Thyra::DiagonalLinearOpBase<double> >(state->invMass_)->getDiag();
      RCP<Thyra::VectorBase<double> > h = Thyra::createMember(iQu->space());

      Thyra::put_scalar(0.0,h.ptr());
      Thyra::ele_wise_prod(1.0,*w,*iQu,h.ptr());
      hScaling_ = Teuchos::rcp(new Thyra::DefaultDiagonalLinearOp<double>(h));
   } 

   LinearOp H = hScaling_;
   if(H==Teuchos::null && not isSymmetric_)
      H = state->invMass_;

   // setup the scaling operator
   if(H==Teuchos::null)
      state->BHBt_ = state->BQBt_;
   else {
      RCP<Teuchos::Time> time = Teuchos::TimeMonitor::getNewTimer("InvLSCStrategy::initializeState Build BHBt");
      Teuchos::TimeMonitor timer(*time);

      // compute BHBt
      state->BHBt_ = explicitMultiply(D,H,G,state->BHBt_);
   }

   // if this is a stable discretization...we are done!
   if(not isStabilized) {
      state->addInverse("BQBtmC",state->BQBt_);
      state->addInverse("BHBtmC",state->BHBt_);
      state->gamma_ = 0.0;
      state->alpha_ = 0.0;
      state->aiD_ = Teuchos::null;

      state->setInitialized(true);

      return;
   }

   // for Epetra_CrsMatrix...zero out certain rows: this ensures spectral radius is correct
   LinearOp modF = F;
   const RCP<const Epetra_Operator> epF = Thyra::get_Epetra_Operator(*F);
   if(epF!=Teuchos::null && rowZeroingNeeded_) {
      // try to get a CRS matrix
      const RCP<const Epetra_CrsMatrix> crsF = rcp_dynamic_cast<const Epetra_CrsMatrix>(epF);

      // if it is a CRS matrix get rows that need to be zeroed
      if(crsF!=Teuchos::null) {
         std::vector<int> zeroIndices;
          
         // get rows in need of zeroing
         Teko::Epetra::identityRowIndices(crsF->RowMap(), *crsF,zeroIndices);

         // build an operator that zeros those rows
         modF = Thyra::epetraLinearOp(rcp(new Teko::Epetra::ZeroedOperator(zeroIndices,crsF)));
      }
   }

   // compute gamma
   Teko_DEBUG_MSG("Calculating gamma",10);
   LinearOp iQuF = multiply(state->invMass_,modF);

   // do 6 power iterations to compute spectral radius: EHSST2007 Eq. 4.28
   Teko::LinearOp stabMatrix; // this is the pressure stabilization matrix to use
   state->gamma_ = std::fabs(Teko::computeSpectralRad(iQuF,5e-2,false,eigSolveParam_))/3.0; 
   Teko_DEBUG_MSG("Calculated gamma",10);
   if(userPresStabMat_!=Teuchos::null) {
      Teko::LinearOp invDGl = Teko::getInvDiagonalOp(userPresStabMat_);
      Teko::LinearOp gammaOp = multiply(invDGl,C);
      state->gamma_ *= std::fabs(Teko::computeSpectralRad(gammaOp,5e-2,false,eigSolveParam_));
      stabMatrix = userPresStabMat_;
   } else 
      stabMatrix = C;

   // compute alpha scaled inv(D): EHSST2007 Eq. 4.29
   // construct B_idF_Bt and save it for refilling later: This could reuse BQBt graph
   LinearOp invDiagF = getInvDiagonalOp(F);
   Teko::ModifiableLinearOp modB_idF_Bt = state->getInverse("BidFBt");
   modB_idF_Bt = explicitMultiply(B,invDiagF,Bt,modB_idF_Bt);
   state->addInverse("BidFBt",modB_idF_Bt);
   const LinearOp B_idF_Bt = modB_idF_Bt;

   MultiVector vec_D = getDiagonal(B_idF_Bt); // this memory could be reused
   update(-1.0,getDiagonal(C),1.0,vec_D); // vec_D = diag(B*inv(diag(F))*Bt)-diag(C)
   const LinearOp invD = buildInvDiagonal(vec_D,"inv(D)");

   Teko_DEBUG_MSG("Calculating alpha",10);
   const LinearOp BidFBtidD = multiply<double>(B_idF_Bt,invD);
   double num = std::fabs(Teko::computeSpectralRad(BidFBtidD,5e-2,false,eigSolveParam_));
   Teko_DEBUG_MSG("Calculated alpha",10);
   state->alpha_ = 1.0/num;
   state->aiD_ = Thyra::scale(state->alpha_,invD);

   // now build B*Q*Bt-gamma*C
   Teko::ModifiableLinearOp BQBtmC = state->getInverse("BQBtmC");
   BQBtmC = explicitAdd(state->BQBt_,scale(-state->gamma_,stabMatrix),BQBtmC);
   state->addInverse("BQBtmC",BQBtmC);

   // now build B*H*Bt-gamma*C
   Teko::ModifiableLinearOp BHBtmC = state->getInverse("BHBtmC");
   if(H==Teuchos::null)
      BHBtmC = BQBtmC;
   else {
      BHBtmC = explicitAdd(state->BHBt_,scale(-state->gamma_,stabMatrix),BHBtmC);
   }
   state->addInverse("BHBtmC",BHBtmC);

   Teko_DEBUG_MSG_BEGIN(5)
      DEBUG_STREAM << "LSC Gamma Parameter = " << state->gamma_ << std::endl;
      DEBUG_STREAM << "LSC Alpha Parameter = " << state->alpha_ << std::endl;
   Teko_DEBUG_MSG_END()

   state->setInitialized(true);
}

/** Compute the inverses required for the LSC Schur complement
  *
  * \note This method assumes that the BQBt and BHBt operators have
  *       been constructed.
  */
void InvLSCStrategy::computeInverses(const BlockedLinearOp & A,LSCPrecondState * state) const
{
   Teko_DEBUG_SCOPE("InvLSCStrategy::computeInverses",10);
   Teko_DEBUG_EXPR(Teuchos::Time invTimer(""));

   const LinearOp F  = getBlock(0,0,A);

   /////////////////////////////////////////////////////////

   // (re)build the inverse of F
   Teko_DEBUG_MSG("LSC::computeInverses Building inv(F)",1);
   Teko_DEBUG_EXPR(invTimer.start(true));
   InverseLinearOp invF = state->getInverse("invF");
   if(invF==Teuchos::null) {
      invF = buildInverse(*invFactoryF_,F);
      state->addInverse("invF",invF); 
   } else {
      rebuildInverse(*invFactoryF_,F,invF);
   }
   Teko_DEBUG_EXPR(invTimer.stop());
   Teko_DEBUG_MSG("LSC::computeInverses GetInvF = " << invTimer.totalElapsedTime(),1);

   /////////////////////////////////////////////////////////

   // (re)build the inverse of BQBt 
   Teko_DEBUG_MSG("LSC::computeInverses Building inv(BQBtmC)",1);
   Teko_DEBUG_EXPR(invTimer.start(true));
   const LinearOp BQBt = state->getInverse("BQBtmC");
   InverseLinearOp invBQBt = state->getInverse("invBQBtmC");
   if(invBQBt==Teuchos::null) {
      invBQBt = buildInverse(*invFactoryS_,BQBt);
      state->addInverse("invBQBtmC",invBQBt); 
   } else {
      rebuildInverse(*invFactoryS_,BQBt,invBQBt);
   }
   Teko_DEBUG_EXPR(invTimer.stop());
   Teko_DEBUG_MSG("LSC::computeInverses GetInvBQBt = " << invTimer.totalElapsedTime(),1);

   /////////////////////////////////////////////////////////

   // Compute the inverse of BHBt or just use BQBt
   ModifiableLinearOp invBHBt = state->getInverse("invBHBtmC");
   if(hScaling_!=Teuchos::null || not isSymmetric_) {
      // (re)build the inverse of BHBt 
      Teko_DEBUG_MSG("LSC::computeInverses Building inv(BHBtmC)",1);
      Teko_DEBUG_EXPR(invTimer.start(true));
      const LinearOp BHBt = state->getInverse("BHBtmC");
      if(invBHBt==Teuchos::null) {
         invBHBt = buildInverse(*invFactoryS_,BHBt);
         state->addInverse("invBHBtmC",invBHBt); 
      } else {
         rebuildInverse(*invFactoryS_,BHBt,invBHBt);
      }
      Teko_DEBUG_EXPR(invTimer.stop());
      Teko_DEBUG_MSG("LSC::computeInverses GetInvBHBt = " << invTimer.totalElapsedTime(),1);
   } 
   else if(invBHBt==Teuchos::null) {
      // just use the Q version
      state->addInverse("invBHBtmC",invBQBt); 
   }
}

//! Initialize from a parameter list
void InvLSCStrategy::initializeFromParameterList(const Teuchos::ParameterList & pl,const InverseLibrary & invLib) 
{
   // get string specifying inverse
   std::string invStr="", invVStr="", invPStr="";
   bool rowZeroing = true;
   bool useLDU = false;
   scaleType_ = Diagonal;

   // "parse" the parameter list
   if(pl.isParameter("Inverse Type"))
      invStr = pl.get<std::string>("Inverse Type");
   if(pl.isParameter("Inverse Velocity Type"))
      invVStr = pl.get<std::string>("Inverse Velocity Type");
   if(pl.isParameter("Inverse Pressure Type")) 
      invPStr = pl.get<std::string>("Inverse Pressure Type");
   if(pl.isParameter("Ignore Boundary Rows"))
      rowZeroing = pl.get<bool>("Ignore Boundary Rows");
   if(pl.isParameter("Use LDU"))
      useLDU = pl.get<bool>("Use LDU");
   if(pl.isParameter("Use Mass Scaling"))
      useMass_ = pl.get<bool>("Use Mass Scaling");
   // if(pl.isParameter("Use Lumping"))
   //    useLumping_ = pl.get<bool>("Use Lumping");
   if(pl.isParameter("Use W-Scaling"))
      useWScaling_ = pl.get<bool>("Use W-Scaling");
   if(pl.isParameter("Eigen Solver Iterations"))
      eigSolveParam_ = pl.get<int>("Eigen Solver Iterations");
   if(pl.isParameter("Scaling Type")) {
      scaleType_ = getDiagonalType(pl.get<std::string>("Scaling Type"));
      TEST_FOR_EXCEPT(scaleType_==NotDiag);
   }

   Teko_DEBUG_MSG_BEGIN(5)
      DEBUG_STREAM << "LSC Inverse Strategy Parameters: " << std::endl;
      DEBUG_STREAM << "   inv type   = \"" << invStr  << "\"" << std::endl;
      DEBUG_STREAM << "   inv v type = \"" << invVStr << "\"" << std::endl;
      DEBUG_STREAM << "   inv p type = \"" << invPStr << "\"" << std::endl;
      DEBUG_STREAM << "   bndry rows = " << rowZeroing << std::endl;
      DEBUG_STREAM << "   use ldu    = " << useLDU << std::endl;
      DEBUG_STREAM << "   use mass    = " << useMass_ << std::endl;
      // DEBUG_STREAM << "   use lumping    = " << useLumping_ << std::endl;
      DEBUG_STREAM << "   use w-scaling    = " << useWScaling_ << std::endl;
      DEBUG_STREAM << "   scale type    = " << getDiagonalName(scaleType_) << std::endl;
      DEBUG_STREAM << "LSC  Inverse Strategy Parameter list: " << std::endl;
      pl.print(DEBUG_STREAM);
   Teko_DEBUG_MSG_END()

   // set defaults as needed
   if(invStr=="") invStr = "Amesos";
   if(invVStr=="") invVStr = invStr;
   if(invPStr=="") invPStr = invStr;

   // build velocity inverse factory
   invFactoryF_ = invLib.getInverseFactory(invVStr);
   invFactoryS_ = invFactoryF_; // by default these are the same
   if(invVStr!=invPStr) // if different, build pressure inverse factory
      invFactoryS_ = invLib.getInverseFactory(invPStr);

   // set other parameters
   setUseFullLDU(useLDU);
   setRowZeroing(rowZeroing);
}

//! For assiting in construction of the preconditioner
Teuchos::RCP<Teuchos::ParameterList> InvLSCStrategy::getRequestedParameters() const 
{
   Teuchos::RCP<Teuchos::ParameterList> result;
   Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList());

   // grab parameters from F solver
   RCP<Teuchos::ParameterList> fList = invFactoryF_->getRequestedParameters();
   if(fList!=Teuchos::null) {
      Teuchos::ParameterList::ConstIterator itr;
      for(itr=fList->begin();itr!=fList->end();++itr)
         pl->setEntry(itr->first,itr->second);
      result = pl;
   }

   // grab parameters from S solver
   RCP<Teuchos::ParameterList> sList = invFactoryS_->getRequestedParameters();
   if(sList!=Teuchos::null) {
      Teuchos::ParameterList::ConstIterator itr;
      for(itr=sList->begin();itr!=sList->end();++itr)
         pl->setEntry(itr->first,itr->second);
      result = pl;
   }

   // use the mass matrix
   if(useMass_) {
      pl->set<Teko::LinearOp>("Velocity Mass Matrix", Teuchos::null,"Velocity mass matrix");
      result = pl;
   }

   // use the mass matrix
   if(useWScaling_) {
      pl->set<Teko::LinearOp>("W-Scaling Vector", Teuchos::null,"W-Scaling Vector");
      result = pl;
   }

   return result;
}

//! For assiting in construction of the preconditioner
bool InvLSCStrategy::updateRequestedParameters(const Teuchos::ParameterList & pl) 
{
   Teko_DEBUG_SCOPE("InvLSCStrategy::updateRequestedParameters",10);
   bool result = true;
 
   // update requested parameters in solvers
   result &= invFactoryF_->updateRequestedParameters(pl);
   result &= invFactoryS_->updateRequestedParameters(pl);

   // set the mass matrix: throw if the strategy is not the right type
   if(useMass_) {
      Teko::LinearOp mass = pl.get<Teko::LinearOp>("Velocity Mass Matrix");

      // we must have a mass matrix
      if(mass==Teuchos::null)
         result &= false;
      else
         setMassMatrix(mass);
   }

   // use W scaling matrix
   if(useWScaling_) {
      Teko::MultiVector wScale = pl.get<Teko::MultiVector>("W-Scaling Vector");

      if(wScale==Teuchos::null)
         result &= false;
      else
         setWScaling(wScale);
   }

   return result;
}

} // end namespace NS
} // end namespace Teko