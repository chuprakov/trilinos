// //////////////////////////////////////////////////////
// TSFCoreNonlinEpetraNPFO.cpp

#include "TSFCoreNonlinEpetraNPFO.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFCoreEpetraVector.hpp"
#include "TSFCoreEpetraMultiVector.hpp"
#include "TSFCoreEpetraLinearOp.hpp"
#include "AbstractFactoryStd.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Ifpack_CrsRiluk.h"

#ifdef _DEBUG
#define EpetraNPFO_VALIDATE_L_IN_RANGE(l) TEST_FOR_EXCEPTION( l > epetra_np_->Nu() || l < 0, std::invalid_argument, "Error!  l == " << l << " out of range")
#else
#define EpetraNPFO_VALIDATE_L_IN_RANGE(l)
#endif

namespace TSFCore {
namespace Nonlin {

//
// EpetraNPFO
//

// Constructors / Initializers / accessors

EpetraNPFO::EpetraNPFO(
	const int      maxLinSolveIter
	,const double  relLinSolveTol
	)
	:maxLinSolveIter_(maxLinSolveIter)
	,relLinSolveTol_(relLinSolveTol)
	,isInitialized_(false)
{
	// Turn off output from Aztec by default!
	aztecOO_.SetAztecOption(AZ_output,AZ_none);
}

EpetraNPFO::EpetraNPFO(
  const Teuchos::RefCountPtr<Epetra::NonlinearProblemFirstOrder>   &epetra_np
  )
{
  initialize(epetra_np);
}

void EpetraNPFO::initialize(
  const Teuchos::RefCountPtr<Epetra::NonlinearProblemFirstOrder>   &epetra_np
  )
{
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;

  epetra_np_ = epetra_np;

  // Make sure that the Epetra Nonlinear Problem is initialized
  epetra_np_->initialize(false);

  // The the number of auxiliary variables and responce functions
  const int
    Nu  = epetra_np_->Nu(),
    nrf = epetra_np_->numResponseFunctions();
  
  // Set up vector spaces
  space_y_ = rcp(new EpetraVectorSpace(epetra_np_->map_y()));
  space_c_ = rcp(new EpetraVectorSpace(epetra_np_->map_c()));
  space_u_.resize(Nu); for(int l=1;l<=Nu;++l) space_u_[l-1] = rcp(new EpetraVectorSpace(epetra_np_->map_u(l)));
  space_g_ = rcp(new EpetraVectorSpace(epetra_np_->map_g()));

  // Set up variable bounds and initial guesses

  // y
  yL_ = rcp(new EpetraVector(rcp(&const_cast<Epetra_Vector&>(epetra_np_->yL()),false),space_y_));
  yU_ = rcp(new EpetraVector(rcp(&const_cast<Epetra_Vector&>(epetra_np_->yU()),false),space_y_));
  y0_ = rcp(new EpetraVector(rcp(&const_cast<Epetra_Vector&>(epetra_np_->y0()),false),space_y_));
  // u
  uL_.resize(Nu);
  uU_.resize(Nu);
  u0_.resize(Nu);
  u_in_.resize(Nu);
  for(int l=1;l<=Nu;++l) {
    uL_[l-1] = rcp(new EpetraVector(rcp(&const_cast<Epetra_Vector&>(epetra_np_->uL(l)),false),space_u_[l-1]));
    uU_[l-1] = rcp(new EpetraVector(rcp(&const_cast<Epetra_Vector&>(epetra_np_->uU(l)),false),space_u_[l-1]));
    u0_[l-1] = rcp(new EpetraVector(rcp(&const_cast<Epetra_Vector&>(epetra_np_->u0(l)),false),space_u_[l-1]));
  }
  // g
  gL_ = rcp(new EpetraVector(rcp(&const_cast<Epetra_Vector&>(epetra_np_->gL()),false),space_g_));
  gU_ = rcp(new EpetraVector(rcp(&const_cast<Epetra_Vector&>(epetra_np_->gU()),false),space_g_));

  // Setup factories
	factory_DcDy_ = Teuchos::rcp(new MemMngPack::AbstractFactoryStd<LinearOpWithSolve<Scalar>,LinearOpWithSolveAztecOO>());
  factory_DcDu_.resize(Nu);
  for(int l=1;l<=Nu;++l) {
    typedef LinearOp<Scalar> T_itfc;
    typedef MemMngPack::PostModNothing<T_itfc> pm_t;
    typedef MemMngPack::AbstractFactoryStd<T_itfc,T_itfc,pm_t,DcDu_Allocator>  af_t;
    factory_DcDu_[l-1] = Teuchos::rcp( new af_t( pm_t(), DcDu_Allocator( epetra_np_->use_DcDu_op(l) ) ) );
  }

  // Setup cache arrays for dealing with DcDu
  DcDu_op_.resize(Nu);
  DcDu_mv_.resize(Nu);
  DcDu_updated_.resize(Nu);
  epetra_DcDu_op_.resize(Nu);
  epetra_DcDu_mv_.resize(Nu);
  epetra_DcDu_args_.resize(Nu);

  // Setup cache arrays for dealing with D\gDu
  DgDu_.resize(Nu);
  DgDu_updated_.resize(Nu);
  epetra_DgDu_.resize(Nu);
  epetra_DgDu_args_.resize(Nu);

  // ToDo: Initialize everything else!

  unsetQuantities();

  isInitialized_ = true;

}

// Overridden from NonlinearProblem

void EpetraNPFO::initialize( bool testSetup )
{
  epetra_np_->initialize(testSetup);
  // Already initialized!
}

bool EpetraNPFO::isInitialized() const
{
  return isInitialized_;
}

int EpetraNPFO::Nu() const
{
  return epetra_np_->Nu();
}

Teuchos::RefCountPtr<const VectorSpace<Scalar> >
EpetraNPFO::space_y() const
{
	return space_y_;
}

Teuchos::RefCountPtr<const VectorSpace<Scalar> >
EpetraNPFO::space_u(int l) const
{
  EpetraNPFO_VALIDATE_L_IN_RANGE(l);
  return space_u_[l-1];
}

Teuchos::RefCountPtr<const VectorSpace<Scalar> >
EpetraNPFO::space_c() const
{
  return space_c_;
}

Teuchos::RefCountPtr<const VectorSpace<Scalar> >
EpetraNPFO::space_g() const
{
	return space_g_;
}

const Vector<Scalar>&
EpetraNPFO::yL() const
{
	return *yL_;
}

const Vector<Scalar>&
EpetraNPFO::yU() const
{
	return *yU_;
}

const Vector<Scalar>&
EpetraNPFO::uL(int l) const
{
  EpetraNPFO_VALIDATE_L_IN_RANGE(l);
  return *uL_[l-1];
}

const Vector<Scalar>&
EpetraNPFO::uU(int l) const
{
  EpetraNPFO_VALIDATE_L_IN_RANGE(l);
  return *uU_[l-1];
}

const Vector<Scalar>&
EpetraNPFO::gL() const
{
	return *gL_;
}

const Vector<Scalar>&
EpetraNPFO::gU() const
{
	return *gU_;
}

const Vector<Scalar>&
EpetraNPFO::y0() const
{
  return *y0_;
}

const Vector<Scalar>&
EpetraNPFO::u0(int l) const
{
  EpetraNPFO_VALIDATE_L_IN_RANGE(l);
  return *u0_[l-1];
}

void EpetraNPFO::set_c(Vector<Scalar>* c)
{
  using DynamicCastHelperPack::dyn_cast;
  if(c) {
    c_ = &dyn_cast<TSFCore::EpetraVector>(*c);
    c_updated_ = false;
  }
  else {
    c_ = NULL;
  }
}

Vector<Scalar>* EpetraNPFO::get_c()
{
  return c_;
}

void EpetraNPFO::set_g(Vector<Scalar>* g)
{
  using DynamicCastHelperPack::dyn_cast;
  if(g) {
    g_ = &dyn_cast<TSFCore::EpetraVector>(*g);
    g_updated_ = false;
  }
  else {
    g_ = NULL;
  }
}

Vector<Scalar>* EpetraNPFO::get_g()
{
  return g_;
}

void EpetraNPFO::unsetQuantities()
{
  c_ = NULL;
  g_ = NULL;
  DcDy_ = NULL;
  std::fill( DcDu_op_.begin(), DcDu_op_.end(), (EpetraLinearOp*)NULL );
  std::fill( DcDu_mv_.begin(), DcDu_mv_.end(), (EpetraMultiVector*)NULL );
  DgDy_ = NULL;
  std::fill( DgDu_.begin(), DgDu_.end(), (EpetraMultiVector*)NULL );
}

void EpetraNPFO::calc_c(
	const Vector<Scalar>     &y
	,const Vector<Scalar>*   u[]
	,bool                    newPoint
	) const
{
  calc_Dc(y,u,newPoint,false);
}

void EpetraNPFO::calc_g(
	const Vector<Scalar>     &y
	,const Vector<Scalar>*   u[]
	,bool                    newPoint
	) const
{
  calc_Dg(y,u,newPoint,false);
}

// Overridden from NonlinearProblemFirstOrder

Teuchos::RefCountPtr< const MemMngPack::AbstractFactory<LinearOpWithSolve<Scalar> > >
EpetraNPFO::factory_DcDy() const
{
  return factory_DcDy_;
}

Teuchos::RefCountPtr< const MemMngPack::AbstractFactory<LinearOp<Scalar > > >
EpetraNPFO::factory_DcDu(int l) const
{
  EpetraNPFO_VALIDATE_L_IN_RANGE(l);
  return factory_DcDu_[l-1];
}

ETransp EpetraNPFO::opDcDy() const
{
  return NOTRANS;
}

ETransp EpetraNPFO::opDcDu(int l) const
{
  return NOTRANS;
}

void EpetraNPFO::set_DcDy(LinearOpWithSolve<Scalar>* DcDy)
{
  using DynamicCastHelperPack::dyn_cast;
  if(DcDy) {
    DcDy_ = &dyn_cast<TSFCore::Nonlin::LinearOpWithSolveAztecOO>(*DcDy);
    DcDy_updated_ = false;
  }
  else {
    DcDy_ = NULL;
  }
}

LinearOpWithSolve<Scalar>* EpetraNPFO::get_DcDy()
{
  return DcDy_;
}

void EpetraNPFO::set_DcDu(int l, LinearOp<Scalar>* DcDu_l)
{
  EpetraNPFO_VALIDATE_L_IN_RANGE(l);
  using DynamicCastHelperPack::dyn_cast;
  if(DcDu_l) {
    if(epetra_np_->use_DcDu_op(l)) {
      DcDu_op_[l-1] = &dyn_cast<TSFCore::EpetraLinearOp>(*DcDu_l);
    }
    else {
      DcDu_mv_[l-1] = &dyn_cast<TSFCore::EpetraMultiVector>(*DcDu_l);
    }
    DcDu_updated_[l-1] = false;
  }
  else {
    DcDu_op_[l-1] = NULL;
    DcDu_mv_[l-1] = NULL;
  }
}

LinearOp<Scalar>* EpetraNPFO::get_DcDu(int l)
{
  EpetraNPFO_VALIDATE_L_IN_RANGE(l);
  if(epetra_np_->use_DcDu_op(l)) {
    return DcDu_op_[l-1];
  }
  return DcDu_mv_[l-1];
}

void EpetraNPFO::set_DgDy(MultiVector<Scalar>* DgDy)
{
  using DynamicCastHelperPack::dyn_cast;
  if(DgDy) {
    DgDy_ = &dyn_cast<TSFCore::EpetraMultiVector>(*DgDy);
    DgDy_updated_ = false;
  }
  else {
    DgDy_ = NULL;
  }
}

MultiVector<Scalar>* EpetraNPFO::get_DgDy()
{
  return DgDy_;
}

void EpetraNPFO::set_DgDu(int l, MultiVector<Scalar>* DgDu_l)
{
  EpetraNPFO_VALIDATE_L_IN_RANGE(l);
  using DynamicCastHelperPack::dyn_cast;
  if(DgDu_l) {
    DgDu_[l-1] = &dyn_cast<TSFCore::EpetraMultiVector>(*DgDu_l);
    DgDu_updated_[l-1] = false;
  }
  else {
    DgDu_[l-1] = NULL;
  }
}

MultiVector<Scalar>* EpetraNPFO::get_DgDu(int l)
{
  return DgDu_[l-1];
}

void EpetraNPFO::calc_DcDy(
	const Vector<Scalar>     &y
	,const Vector<Scalar>*   u[]
	,bool                    newPoint
	) const
{
  calc_Dc(y,u,newPoint,true);
}

void EpetraNPFO::calc_DcDu(
	int                      l
	,const Vector<Scalar>    &y
	,const Vector<Scalar>*   u[]
	,bool                    newPoint
	) const
{
  calc_Dc(y,u,newPoint,true);
}

void EpetraNPFO::calc_DgDy(
	const Vector<Scalar>     &y
	,const Vector<Scalar>*   u[]
	,bool                    newPoint
	) const
{
  calc_Dg(y,u,newPoint,true);
}

void EpetraNPFO::calc_DgDu(
	int                      l
	,const Vector<Scalar>    &y
	,const Vector<Scalar>*   u[]
	,bool                    newPoint
	) const
{
  calc_Dg(y,u,newPoint,true);
}

// private

const Epetra_Vector**
EpetraNPFO::set_u( const Vector<Scalar>* u_in[], bool newPoint ) const
{
	using DynamicCastHelperPack::dyn_cast;
  updateNewPoint(newPoint);
  if( u_in && newPoint ) {
    for( int l = 1; l <= u_in_.size(); ++l ) {
      if(u_in[l-1]) {
        u_in_[l-1] = &*(dyn_cast<const EpetraVector>(*u_in[l-1]).epetra_vec());
      }
      else {
        u_in_[l-1] = NULL;  // Tell epetra_np_ to use epetra_np_->u0(l)
      }
    }
  }
  return &u_in_[0];
}

void EpetraNPFO::updateNewPoint( bool newPoint ) const
{
	if(newPoint) {
		const int Nu = epetra_np_->Nu();
		c_updated_ = g_updated_ = DcDy_updated_ = DgDy_updated_ = false;
		for(int l=1;l<=Nu;++l) DcDu_updated_[l-1] = DgDu_updated_[l-1] = false;
	}
}

void EpetraNPFO::calc_Dc(
  const Vector<Scalar>     &y_in
  ,const Vector<Scalar>*   u_in[]
  ,bool                    newPoint
  ,bool                    computeGradients
  ) const
{
  const int Nu = epetra_np_->Nu();
  //
  // Get Epetra objects for y and u
  //
  const Epetra_Vector &y  = get_epetra_vec(y_in);
  const Epetra_Vector **u = set_u( u_in, newPoint );
  //
  // Pick out the raw Epetra objects and in the process set their
  // adpater objects to uninitialized.
	//
	bool computeSomething = false;
  // c
  Teuchos::RefCountPtr<Epetra_Vector>  epetra_c;
  if(c_ && !c_updated_) {
		c_->setUninitialized( &epetra_c );
		computeSomething = true;
	}
  // DcDy
  Teuchos::RefCountPtr<Epetra_Operator>  epetra_DcDy_op;
  Teuchos::RefCountPtr<Epetra_Operator>  epetra_DcDy_prec;
  const bool allow_specialized_DcDy_prec = true; // ToDo: Make this an external parameter!
  const bool specialized_DcDy_prec = ( epetra_np_->specialized_DcDy_prec() && allow_specialized_DcDy_prec );
  if(computeGradients && DcDy_ && !DcDy_updated_) {
    DcDy_->setUninitialized( &epetra_DcDy_op, NULL, NULL, &epetra_DcDy_prec, NULL );
    if( !epetra_DcDy_op.get() ) {
      epetra_DcDy_op = epetra_np_->create_DcDy_op();
			DcDy_->set_trace_out(
				Teuchos::rcp(new std::ofstream("LinearOpWithSolveAztecOO.out") )
				); // ToDo: Make this more flexible!
    }
    if( !epetra_DcDy_prec.get() && specialized_DcDy_prec ) {
      epetra_DcDy_op = epetra_np_->create_DcDy_prec();
    }
		computeSomething = true;
  }
  // DcDu
  for(int l=1;l<=Nu;++l) {
    if(computeGradients && DcDu_op_[l-1] && !DcDu_updated_[l-1]) {
      DcDu_op_[l-1]->setUninitialized( &epetra_DcDu_op_[l-1] ); // Grab the RCP<Epetra_Operator>
      if(!epetra_DcDu_op_[l-1].get())
        epetra_DcDu_op_[l-1] = epetra_np_->create_DcDu_op(l);
      epetra_DcDu_args_[l-1] = Epetra::EpetraOp_or_EpetraMV(&*epetra_DcDu_op_[l-1]);
			computeSomething = true;
    }
    else if(computeGradients && DcDu_mv_[l-1] && !DcDu_updated_[l-1]) {
      DcDu_mv_[l-1]->setUninitialized( &epetra_DcDu_mv_[l-1] ); // Grap the RCP<Epetra_MultiVector>
      if(!epetra_DcDu_mv_[l-1].get())
        epetra_DcDu_mv_[l-1] = epetra_np_->create_DcDu_mv(l);
      epetra_DcDu_args_[l-1] = Epetra::EpetraOp_or_EpetraMV(&*epetra_DcDu_mv_[l-1]);
			computeSomething = true;
    }
    else {
      epetra_DcDu_args_[l-1] = Epetra::EpetraOp_or_EpetraMV(); // No DcDu(l) to compute!
    }
  }

	if(computeSomething) {
		
		//
		// Compute the set Epetra objects
		//
		epetra_np_->calc_Dc(
			y
			,u
			,epetra_c.get()
			,epetra_DcDy_op.get()
			,specialized_DcDy_prec ? &*epetra_DcDy_prec : NULL
			,&epetra_DcDu_args_[0]
			);
		
		//
		// Put the raw Epetra objects back into the adpater objects
		//
		// c
		if( epetra_c.get() ) {
			c_->initialize( epetra_c, space_c_ );
			c_updated_ = true;
		}
		// DcDy
		if( epetra_DcDy_op.get() ) {
			// Setup the externally defined preconditioner
			const bool usePrec = true; // ToDo: Make an external option
			if( usePrec && !specialized_DcDy_prec ) {
				precGenerator().setupPrec(epetra_DcDy_op,&epetra_DcDy_prec);
			}
			// Set the options for AztecOO::Iterate(...)
			DcDy_->maxIter(maxLinSolveIter());
			DcDy_->relTol(relLinSolveTol());
			// Finally initialize the aggregate LinearOpWithSolve object
			DcDy_->initialize(
				epetra_DcDy_op                                               // Op
				,epetra_np_->opDcDy() == Epetra::NOTRANS ? NOTRANS : TRANS   // Op_trans
				,Teuchos::rcp(&aztecOO_,false)                               // solver
				,epetra_DcDy_prec                                            // Prec
				,epetra_np_->opDcDy() == Epetra::NOTRANS ? NOTRANS : TRANS   // Prec_trans
				,epetra_np_->adjointSupported()                              // adjointSupported
				);
			DcDy_updated_ = true;
		}
		// DcDu
		for(int l=1;l<=Nu;++l) {
			if(epetra_DcDu_op_[l-1].get()) {
				DcDu_op_[l-1]->initialize(epetra_DcDu_op_[l-1]);
				epetra_DcDu_op_[l-1] = Teuchos::null;
				DcDu_updated_[l-1] = true;
			}
			else if(epetra_DcDu_mv_[l-1].get()) {
				DcDu_mv_[l-1]->initialize(epetra_DcDu_mv_[l-1],space_c_);
				epetra_DcDu_mv_[l-1] = Teuchos::null;
				DcDu_updated_[l-1] = true;
			}
		}
	}

}

void EpetraNPFO::calc_Dg(
  const Vector<Scalar>     &y_in
  ,const Vector<Scalar>*   u_in[]
  ,bool                    newPoint
  ,bool                    computeGradients
  ) const
{
  const int Nu = epetra_np_->Nu();
  //
  // Get Epetra objects for y and u
  //
  const Epetra_Vector &y  = get_epetra_vec(y_in);
  const Epetra_Vector **u = set_u( u_in, newPoint );
  //
  // Pick out the raw Epetra objects and in the process set their
  // adpater objects to uninitialized.
  //
	bool computeSomething = false;
  // g
  Teuchos::RefCountPtr<Epetra_Vector>  epetra_g;
  if(g_ && !g_updated_ ) {
		g_->setUninitialized( &epetra_g );
		computeSomething = true;
	}
  // DgDy
  Teuchos::RefCountPtr<Epetra_MultiVector>  epetra_DgDy;
  if(computeGradients && DgDy_ && !DgDy_updated_) {
    DgDy_->setUninitialized( &epetra_DgDy );
    if( !epetra_DgDy.get() ) {
      epetra_DgDy = Teuchos::rcp(new Epetra_MultiVector(*epetra_np_->map_y(),space_g_->dim()));
    }
		computeSomething = true;
  }
  // DgDu
  for(int l=1;l<=Nu;++l) {
    if(computeGradients && DgDu_[l-1] && !DgDu_updated_[l-1]) {
      DgDu_[l-1]->setUninitialized( &epetra_DgDu_[l-1] ); // Grap the RCP<Epetra_MultiVector>
      if(!epetra_DgDu_[l-1].get())
        epetra_DgDu_[l-1] = Teuchos::rcp(new Epetra_MultiVector(*epetra_np_->map_u(l),space_g_->dim()));
      epetra_DgDu_args_[l-1] = &*epetra_DgDu_[l-1];
			computeSomething = true;
    }
    else {
      epetra_DgDu_args_[l-1] = NULL;
    }
  }
	
	if(computeSomething) {

		//
		// Compute the set Epetra objects
		//
		epetra_np_->calc_Dg(
			y
			,u
			,epetra_g.get()
			,epetra_DgDy.get()
			,&epetra_DgDu_args_[0]
			);
		
		//
		// Put the raw Epetra objects back into the adpater objects
		//
		if( epetra_g.get() ) {
			g_->initialize( epetra_g, space_g_ );
			g_updated_ = true;
		}
		// DgDy
		if(epetra_DgDy.get()) {
			DgDy_->initialize( epetra_DgDy, space_y_ );
			epetra_DgDy = Teuchos::null;
			DgDy_updated_ = true;
		}
		// DgDu
		for(int l=1;l<=Nu;++l) {
			if(epetra_DgDu_[l-1].get()) {
				DgDu_[l-1]->initialize( epetra_DgDu_[l-1], space_u_[l-1] );
				epetra_DgDu_[l-1] = Teuchos::null;
				DgDu_updated_[l-1] = true;
			}
		}

	}

}

//
// EpetraNPFO::DcDu_Allocator
//

EpetraNPFO::DcDu_Allocator::DcDu_Allocator(
  const bool  useEO
  )
  :useEO_(useEO)
{}

const EpetraNPFO::DcDu_Allocator::ptr_t
EpetraNPFO::DcDu_Allocator::allocate() const
{
  if(useEO_) {
    // Wrap an Epetra_Operator
    return Teuchos::rcp(new EpetraLinearOp());
  }
  // Wrap an Epetra_MultiVector
  return Teuchos::rcp(new EpetraMultiVector());
}

} // namespace Nonlin
} // namespace TSFCore
