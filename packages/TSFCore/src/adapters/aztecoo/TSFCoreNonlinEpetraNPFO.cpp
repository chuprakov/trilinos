// //////////////////////////////////////////////////////
// TSFCoreNonlinEpetraNPFO.cpp

#include "TSFCoreNonlinEpetraNPFO.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFCoreEpetraVector.hpp"
#include "TSFCoreEpetraMultiVector.hpp"
#include "TSFCoreEpetraLinearOp.hpp"


#ifdef _DEBUG
#define EpetraNPFO_VALIDATE_L_IN_RANGE(l) TEST_FOR_EXCEPTION( l > epetra_np_->Nu() || l < 0, std::invalid_argument, "Error!  l == " << l << " out of range")
#else
#define EpetraNPFO_VALIDATE_L_IN_RANGE(l)
#endif

namespace TSFCore {
namespace Nonlin {

// Constructors / Initializers / accessors

EpetraNPFO::EpetraNPFO()
  :isInitialized_(false)
{}

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
  for(int l=1;l<=Nu;++l) {
    uL_[l-1] = rcp(new EpetraVector(rcp(&const_cast<Epetra_Vector&>(epetra_np_->uL(l)),false),space_u_[l-1]));
    uU_[l-1] = rcp(new EpetraVector(rcp(&const_cast<Epetra_Vector&>(epetra_np_->uU(l)),false),space_u_[l-1]));
    u0_[l-1] = rcp(new EpetraVector(rcp(&const_cast<Epetra_Vector&>(epetra_np_->u0(l)),false),space_u_[l-1]));
  }
  // g
  gL_ = rcp(new EpetraVector(rcp(&const_cast<Epetra_Vector&>(epetra_np_->gL()),false),space_g_));
  gU_ = rcp(new EpetraVector(rcp(&const_cast<Epetra_Vector&>(epetra_np_->gU()),false),space_g_));

  // ToDo: Initialize everything else!

  unsetQuantities();

  isInitialized_ = true;

}

// Overridden from NonlinearProblem

void EpetraNPFO::initialize( bool testSetup )
{
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
  // ToDo: Unset the rest ...
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
  assert(0);
  return NOTRANS;
}

ETransp EpetraNPFO::opDcDu(int l) const
{
  assert(0);
  return NOTRANS;
}

void EpetraNPFO::set_DcDy(LinearOpWithSolve<Scalar>* DcDy)
{
  assert(0);
}

LinearOpWithSolve<Scalar>* EpetraNPFO::get_DcDy()
{
  assert(0);
  return NULL;
}

void EpetraNPFO::set_DcDu(int l, LinearOp<Scalar>* DcDu_l)
{
  assert(0);
}

LinearOp<Scalar>* EpetraNPFO::get_DcDu(int l)
{
  assert(0);
  return NULL;
}

void EpetraNPFO::set_DgDy(MultiVector<Scalar>* DgDy)
{
  assert(0);
}

MultiVector<Scalar>* EpetraNPFO::get_DgDy()
{
  assert(0);
  return NULL;
}

void EpetraNPFO::set_DgDu(int l, MultiVector<Scalar>* DgDu_l)
{
  assert(0);
}

MultiVector<Scalar>* EpetraNPFO::get_DgDu(int l)
{
  assert(0);
  return NULL;
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

void EpetraNPFO::set_u( const Vector<Scalar>* u_in[], bool newPoint ) const
{
  // ToDo: Fill this in!
}

void EpetraNPFO::calc_Dc(
  const Vector<Scalar>     &y_in
  ,const Vector<Scalar>*   u_in[]
  ,bool                    newPoint
  ,bool                    computeGradients
  ) const
{
  const Epetra_Vector &y = get_epetra_vec(y_in);
  set_u( u_in, newPoint );

  //
  // Pick out the raw Epetra objects and in the process set their
  // adpater objects to uninitialized.
  //
  Teuchos::RefCountPtr<Epetra_Vector>  epetra_c;
  if(c_ && !c_updated_) c_->setUninitialized( &epetra_c );

  //
  // Compute the set Epetra objects
  //
  epetra_np_->calc_Dc(
    y
    ,NULL            // ToDo: Put in array for u
    ,epetra_c.get()
    ,NULL            // ToDo: Put in DcDy
    ,NULL            // ToDo: Put in array for DcDy
    );
  // ToDo: Finish for DcDy, DcDu ...

  //
  // Put the raw Epetra objects back into the adpater objects
  //
  if( epetra_c.get() ) {
    c_->initialize( epetra_c, space_c_ );
    c_updated_ = true;
  }
  // ToDo: Finish for DcDy, DcDu ...

}

void EpetraNPFO::calc_Dg(
  const Vector<Scalar>     &y
  ,const Vector<Scalar>*   u[]
  ,bool                    newPoint
  ,bool                    computeGradients
  ) const
{
  assert(0);
}

} // namespace Nonlin
} // namespace TSFCore

