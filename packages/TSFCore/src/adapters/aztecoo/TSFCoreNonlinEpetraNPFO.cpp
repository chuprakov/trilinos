// //////////////////////////////////////////////////////
// TSFCoreNonlinEpetraNPFO.cpp
//
// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "TSFCoreNonlinEpetraNPFO.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFCoreEpetraVector.hpp"
#include "TSFCoreEpetraMultiVector.hpp"
#include "TSFCoreEpetraLinearOp.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Comm.h"
#include "Ifpack_CrsRiluk.h"
#include "Teuchos_Time.hpp"
#include "Teuchos_dyn_cast.hpp"

#ifdef _DEBUG
#define EpetraNPFO_VALIDATE_L_IN_RANGE(l) TEST_FOR_EXCEPTION( l > epetra_np_->Nu() || l < 0, std::invalid_argument, "Error!  l == " << l << " out of range")
#else
#define EpetraNPFO_VALIDATE_L_IN_RANGE(l)
#endif


namespace {

void sqrt( Epetra_Vector *v )
{
	const int localDim = v->Map().NumMyPoints();
	double *v_ptr = &(*v)[0];
	for( int i = 0; i < localDim; ++i, ++v_ptr )
		*v_ptr = std::sqrt(*v_ptr);
}

} // namespace

namespace TSFCore {
namespace Nonlin {

//
// EpetraNPFO
//

// Constructors / Initializers / accessors

EpetraNPFO::EpetraNPFO(
	const bool           autoScaleStateConstraints
	,const bool          autoScaleStateVariables
	,const int           maxLinSolveIter
	,const double        relLinSolveTol
	,const bool          usePrec
	,const bool          testOperators
	,const std::string   &yGuessFileNameBase
	,const std::string   &yFinalFileNameBase
	)
	:autoScaleStateConstraints_(autoScaleStateConstraints)
	,autoScaleStateVariables_(autoScaleStateVariables)
	,maxLinSolveIter_(maxLinSolveIter)
	,relLinSolveTol_(relLinSolveTol)
	,usePrec_(usePrec)
	,testOperators_(testOperators)
	,yGuessFileNameBase_(yGuessFileNameBase)
	,yFinalFileNameBase_(yFinalFileNameBase)
	,isInitialized_(false)
	,numProc_(1)
	,procRank_(0)
{
	// Turn off output from Aztec by default!
	aztecOO_.SetAztecOption(AZ_output,AZ_none);
}

void EpetraNPFO::initialize(
  const Teuchos::RefCountPtr<Epetra::NonlinearProblemFirstOrder>   &epetra_np
  )
{
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
	using Teuchos::dyn_cast;

	Teuchos::Time timer("");

	TEST_FOR_EXCEPTION( !epetra_np.get(), std::invalid_argument, "Error!" );

	const Epetra_Comm &epetra_comm = epetra_np->map_y()->Comm();
	numProc_  = epetra_comm.NumProc();
	procRank_ = epetra_comm.MyPID();
  epetra_np_ = epetra_np;

  // Make sure that the Epetra Nonlinear Problem is initialized
  epetra_np_->initialize(false);

  // Get the the number of auxiliary variables
  const int
    Nu  = epetra_np_->Nu();

  // Setup cache arrays for dealing with DcDu
  DcDu_op_.resize(Nu);
  DcDu_mv_.resize(Nu);
  DcDu_updated_.resize(Nu);
  epetra_DcDu_op_.resize(Nu);
  epetra_DcDu_mv_.resize(Nu);
  epetra_DcDu_args_.resize(Nu);

  // Setup cache arrays for dealing with DgDu
  DgDu_.resize(Nu);
  DgDu_updated_.resize(Nu);
  epetra_DgDu_.resize(Nu);
  epetra_DgDu_args_.resize(Nu);

  unsetQuantities();
  
  // Set up vector spaces
  space_y_ = rcp(new EpetraVectorSpace(epetra_np_->map_y()));
  space_c_ = rcp(new EpetraVectorSpace(epetra_np_->map_c()));
  space_u_.resize(Nu); for(int l=1;l<=Nu;++l) space_u_[l-1] = rcp(new EpetraVectorSpace(epetra_np_->map_u(l)));
  space_g_ = rcp(new EpetraVectorSpace(epetra_np_->map_g()));

	//
  // Set up variable bounds and initial guesses (unscaled)
	//

  // y
  yL_ = rcp_dynamic_cast<EpetraVector>(space_y_->createMember()); *yL_->epetra_vec() = epetra_np_->yL();
  yU_ = rcp_dynamic_cast<EpetraVector>(space_y_->createMember()); *yU_->epetra_vec() = epetra_np_->yU();
  y0_ = rcp_dynamic_cast<EpetraVector>(space_y_->createMember()); *y0_->epetra_vec() = epetra_np_->y0();
	read_y_guess( &*y0_ );
  // u
  uL_.resize(Nu);
  uU_.resize(Nu);
  u0_.resize(Nu);
  u_unscaled_.resize(Nu);
  for(int l=1;l<=Nu;++l) {
    uL_[l-1] = rcp_dynamic_cast<EpetraVector>(space_u_[l-1]->createMember()); *uL_[l-1]->epetra_vec() = epetra_np_->uL(l);
    uU_[l-1] = rcp_dynamic_cast<EpetraVector>(space_u_[l-1]->createMember()); *uU_[l-1]->epetra_vec() = epetra_np_->uU(l);
    u0_[l-1] = rcp_dynamic_cast<EpetraVector>(space_u_[l-1]->createMember()); *u0_[l-1]->epetra_vec() = epetra_np_->u0(l);
		// ToDo: Optionally read in u0[l] from file!
  }
  // g
  gL_ = rcp(new EpetraVector(rcp(&const_cast<Epetra_Vector&>(epetra_np_->gL()),false),space_g_));
  gU_ = rcp(new EpetraVector(rcp(&const_cast<Epetra_Vector&>(epetra_np_->gU()),false),space_g_));
	
	//
	// Compute automatic scaling vectors for c and y
	//
	computeScaling();

	//
	// Rescale the variables and bounds
	//
	scale_y( *y0_->epetra_vec(), &*y0_->epetra_vec() );
	// ToDo: Most modify yL and yU for scaling and maintain constants for unbounded values!

  // Setup factories
	factory_DcDy_ = Teuchos::rcp(new Teuchos::AbstractFactoryStd<LinearOpWithSolve<Scalar>,LinearOpWithSolveAztecOO>());
  factory_DcDu_.resize(Nu);
  for(int l=1;l<=Nu;++l) {
    typedef LinearOp<Scalar> T_itfc;
    typedef Teuchos::PostModNothing<T_itfc> pm_t;
    typedef Teuchos::AbstractFactoryStd<T_itfc,T_itfc,pm_t,DcDu_Allocator>  af_t;
    factory_DcDu_[l-1] = Teuchos::rcp( new af_t( pm_t(), DcDu_Allocator( epetra_np_->use_DcDu_op(l) ) ) );
  }

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

Teuchos::RefCountPtr<const VectorSpace<EpetraNPFO::Scalar> >
EpetraNPFO::space_y() const
{
	return space_y_;
}

Teuchos::RefCountPtr<const VectorSpace<EpetraNPFO::Scalar> >
EpetraNPFO::space_u(int l) const
{
  EpetraNPFO_VALIDATE_L_IN_RANGE(l);
  return space_u_[l-1];
}

Teuchos::RefCountPtr<const VectorSpace<EpetraNPFO::Scalar> >
EpetraNPFO::space_c() const
{
  return space_c_;
}

Teuchos::RefCountPtr<const VectorSpace<EpetraNPFO::Scalar> >
EpetraNPFO::space_g() const
{
	return space_g_;
}

const Vector<EpetraNPFO::Scalar>&
EpetraNPFO::yL() const
{
	return *yL_;
}

const Vector<EpetraNPFO::Scalar>&
EpetraNPFO::yU() const
{
	return *yU_;
}

const Vector<EpetraNPFO::Scalar>&
EpetraNPFO::uL(int l) const
{
  EpetraNPFO_VALIDATE_L_IN_RANGE(l);
  return *uL_[l-1];
}

const Vector<EpetraNPFO::Scalar>&
EpetraNPFO::uU(int l) const
{
  EpetraNPFO_VALIDATE_L_IN_RANGE(l);
  return *uU_[l-1];
}

const Vector<EpetraNPFO::Scalar>&
EpetraNPFO::gL() const
{
	return *gL_;
}

const Vector<EpetraNPFO::Scalar>&
EpetraNPFO::gU() const
{
	return *gU_;
}

const Vector<EpetraNPFO::Scalar>&
EpetraNPFO::y0() const
{
  return *y0_;
}

const Vector<EpetraNPFO::Scalar>&
EpetraNPFO::u0(int l) const
{
  EpetraNPFO_VALIDATE_L_IN_RANGE(l);
  return *u0_[l-1];
}

void EpetraNPFO::set_c(Vector<Scalar>* c)
{
  using Teuchos::dyn_cast;
  if(c) {
    c_ = &dyn_cast<TSFCore::EpetraVector>(*c);
    c_updated_ = false;
  }
  else {
    c_ = NULL;
  }
}

Vector<EpetraNPFO::Scalar>* EpetraNPFO::get_c()
{
  return c_;
}

void EpetraNPFO::set_g(Vector<Scalar>* g)
{
  using Teuchos::dyn_cast;
  if(g) {
    g_ = &dyn_cast<TSFCore::EpetraVector>(*g);
    g_updated_ = false;
  }
  else {
    g_ = NULL;
  }
}

Vector<EpetraNPFO::Scalar>* EpetraNPFO::get_g()
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

void EpetraNPFO::reportFinalSolution(
		const Vector<Scalar>     &y_in
		,const Vector<Scalar>*   u_in[]
		,bool                    solved
	)
{
  // Get Epetra objects for y and u
  const Epetra_Vector &y  = set_y(y_in);
  const Epetra_Vector **u = set_u( u_in, true );
	// Write solution files
	write_y_final(y);
	// ToDo: Write the final solution for u[l] to file!
	// Report the final solution
	epetra_np_->reportFinalSolution(y,u,solved);
}

// Overridden from NonlinearProblemFirstOrder

Teuchos::RefCountPtr< const Teuchos::AbstractFactory<LinearOpWithSolve<EpetraNPFO::Scalar> > >
EpetraNPFO::factory_DcDy() const
{
  return factory_DcDy_;
}

Teuchos::RefCountPtr< const Teuchos::AbstractFactory<LinearOp<EpetraNPFO::Scalar > > >
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
  using Teuchos::dyn_cast;
  if(DcDy) {
    DcDy_ = &dyn_cast<TSFCore::Nonlin::LinearOpWithSolveAztecOO>(*DcDy);
    DcDy_updated_ = false;
  }
  else {
    DcDy_ = NULL;
  }
}

LinearOpWithSolve<EpetraNPFO::Scalar>* EpetraNPFO::get_DcDy()
{
  return DcDy_;
}

void EpetraNPFO::set_DcDu(int l, LinearOp<Scalar>* DcDu_l)
{
  EpetraNPFO_VALIDATE_L_IN_RANGE(l);
  using Teuchos::dyn_cast;
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

LinearOp<EpetraNPFO::Scalar>* EpetraNPFO::get_DcDu(int l)
{
  EpetraNPFO_VALIDATE_L_IN_RANGE(l);
  if(epetra_np_->use_DcDu_op(l)) {
    return DcDu_op_[l-1];
  }
  return DcDu_mv_[l-1];
}

void EpetraNPFO::set_DgDy(MultiVector<Scalar>* DgDy)
{
  using Teuchos::dyn_cast;
  if(DgDy) {
    DgDy_ = &dyn_cast<TSFCore::EpetraMultiVector>(*DgDy);
    DgDy_updated_ = false;
  }
  else {
    DgDy_ = NULL;
  }
}

MultiVector<EpetraNPFO::Scalar>* EpetraNPFO::get_DgDy()
{
  return DgDy_;
}

void EpetraNPFO::set_DgDu(int l, MultiVector<Scalar>* DgDu_l)
{
  EpetraNPFO_VALIDATE_L_IN_RANGE(l);
  using Teuchos::dyn_cast;
  if(DgDu_l) {
    DgDu_[l-1] = &dyn_cast<TSFCore::EpetraMultiVector>(*DgDu_l);
    DgDu_updated_[l-1] = false;
  }
  else {
    DgDu_[l-1] = NULL;
  }
}

MultiVector<EpetraNPFO::Scalar>* EpetraNPFO::get_DgDu(int l)
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

const Epetra_Vector&
EpetraNPFO::set_y( const Vector<Scalar> &y_scaled_in ) const
{
	using Teuchos::dyn_cast;
	const Epetra_Vector &y_scaled = *dyn_cast<const EpetraVector>(y_scaled_in).epetra_vec();
	if(!y_scaling_.get())
		return y_scaled; // There is no scaling
	// Unscale y_scaled into y_unscaled_ and return
	unscale_y( y_scaled, &*y_unscaled_ );
	return *y_unscaled_;
}

const Epetra_Vector**
EpetraNPFO::set_u( const Vector<Scalar>* u_in[], bool newPoint ) const
{
	// ToDo: Must unscaled u_in into u_unscaled_ if we are scaling u!
	using Teuchos::dyn_cast;
  updateNewPoint(newPoint);
  if( u_in && newPoint ) {
    for( int l = 1; l <= static_cast<int>(u_unscaled_.size()); ++l ) {
      if(u_in[l-1]) {
        u_unscaled_[l-1] = &*(dyn_cast<const EpetraVector>(*u_in[l-1]).epetra_vec());
      }
      else {
        u_unscaled_[l-1] = &*u0_[l-1]->epetra_vec();
      }
    }
  }
  return &u_unscaled_[0];
}

void EpetraNPFO::read_y_guess( EpetraVector *y )
{
	if(yGuessFileNameBase().length()) {
		std::ostringstream y_guess_file_name;
		y_guess_file_name << yGuessFileNameBase() << "." << numProc_ << "." << procRank_;
		std::ifstream y_guess_file(y_guess_file_name.str().c_str());
		TEST_FOR_EXCEPTION(
			!y_guess_file, std::runtime_error
			,"EpetraNPFO::read_y_guess(...): Error, the file \""<<y_guess_file_name.str()<<"\" could not be opended for input!"
			);
		Epetra_Vector &epetra_y = *y->epetra_vec();
		const int local_dim = epetra_np_->map_y()->NumMyElements();
		int read_local_dim;
#ifdef _DEBUG
		TEST_FOR_EXCEPTION( y_guess_file.eof(), std::runtime_error, "Error!" );
#endif
		y_guess_file >> read_local_dim;
		for( int i = 0; i < local_dim; ++i ) {
#ifdef _DEBUG
			TEST_FOR_EXCEPTION( y_guess_file.eof(), std::runtime_error, "Error!" );
#endif
			int i_read = -1;
			y_guess_file >> i_read;
#ifdef _DEBUG
			TEST_FOR_EXCEPTION( i != i_read || y_guess_file.eof(), std::runtime_error, "Error!" );
#endif
			y_guess_file >> epetra_y[i];
		}
	}
}

void EpetraNPFO::write_y_final( const Epetra_Vector &epetra_y )
{
	if(yFinalFileNameBase().length()) {
		std::ostringstream y_final_file_name;
		y_final_file_name << yFinalFileNameBase() << "." << numProc_ << "." << procRank_;
		std::ofstream y_final_file(y_final_file_name.str().c_str());
		TEST_FOR_EXCEPTION(
			!y_final_file, std::runtime_error
			,"EpetraNPFO::write_y_final(...): Error, the file \""<<y_final_file_name.str()<<"\" could not be opended for output!"
			);
		y_final_file.precision(23); // Should print enough digits to read back in same binary number!
		const int local_dim = epetra_np_->map_y()->NumMyElements();
		y_final_file << local_dim << std::endl;
		for( int i = 0; i < local_dim; ++i ) {
			y_final_file << "  " << i << "  " << epetra_y[i] << std::endl;
		}
	}
}

void EpetraNPFO::computeScaling()
{
	using Teuchos::dyn_cast;
	using Teuchos::rcp;
	using Teuchos::rcp_dynamic_cast;
	Teuchos::Time timer("");

	y_unscaled_ = Teuchos::rcp(new Epetra_Vector(*space_y_->epetra_map()));
	// Setup scaling based on initial matrix values
	if( autoScaleStateConstraints() || autoScaleStateVariables() ) {
		if(get_trace_out().get())
			trace_out()
				<< "\nEpetraNPFO::initialize(epetra_np): Computing automatic scaling vectors from DcDy evaluated at y0, {u0(l)} ...\n";
		
		// Get unscaled initial guess
		const Epetra_Vector &y  = *y0_->epetra_vec();
		const Epetra_Vector **u = set_u( NULL, true );

		Teuchos::RefCountPtr<Epetra_Operator>  epetra_DcDy_op = epetra_np_->create_DcDy_op();
		Epetra_RowMatrix &epetra_DcDy_rm = dyn_cast<Epetra_RowMatrix>(*epetra_DcDy_op); // Must support this interface!

		// ToDo: Save this jacobian object for the first Jacobian returned
		// from the abstract factory!  We can then check to see if the
		// first point is the base point (y0,u0) in which case the
		// Jacobian DcDy is already computed.  That way the cost of
		// computing these scaling factors will only include the extra
		// InvRowSum() and InvColSum() calls.

		if(get_trace_out().get())
			trace_out()
				<< "\nEpetraNPFO::initialize(epetra_np): Calling " << typeid(*epetra_np_).name()
				<< ".calc_Dc(...) to compute inital DcDy matrix ...\n";
		timer.start(true);
		epetra_np_->calc_Dc(
			y
			,u
			,NULL                       // c
			,epetra_DcDy_op.get()       // DcDy_op
			,NULL                       // DcDy_prec
			,NULL                       // DcDu
			);
		timer.stop();
		if(get_trace_out().get())
			trace_out()
				<< "\n  => time = " << timer.totalElapsedTime() << " sec\n";
		
		if(get_trace_out().get())
			trace_out()
				<< "\nEpetraNPFO::initialize(epetra_np): "
				<< "Computing row scaling (c_scaling) and column scaling (y_scaling) ...";
		Teuchos::RefCountPtr<EpetraVector>
			rowScaling = rcp_dynamic_cast<EpetraVector>(space_c_->createMember()),
			colScaling = rcp_dynamic_cast<EpetraVector>(space_y_->createMember());
		timer.start(true);
		TEST_FOR_EXCEPTION(0!=epetra_DcDy_rm.InvRowSums(*rowScaling->epetra_vec()),std::runtime_error,"Error!");
		TEST_FOR_EXCEPTION(0!=epetra_DcDy_rm.InvColSums(*colScaling->epetra_vec()),std::runtime_error,"Error!");
		timer.stop();
		if(get_trace_out().get()) {
			double min_row[1], min_col[1];
			TEST_FOR_EXCEPTION(0!=rowScaling->epetra_vec()->MinValue(min_row),std::runtime_error,"Error!");
			TEST_FOR_EXCEPTION(0!=colScaling->epetra_vec()->MinValue(min_col),std::runtime_error,"Error!");
			trace_out()
				<< "\n  => time = " << timer.totalElapsedTime() << " sec\n"
				<< "\n  max{rowScaling} = " << norm_inf(*rowScaling)
				<< "\n  min{rowScaling} = " << min_row[0]
				<< "\n  max{colScaling} = " << norm_inf(*colScaling)
				<< "\n  min{colScaling} = " << min_col[0]
				<< std::endl;
		}

		if( autoScaleStateConstraints() ) {
			if(get_trace_out().get())
				trace_out() << "\nEpetraNPFO::initialize(epetra_np): autoScaleStateConstraints==true: setting c_scaling = rowScaling ...\n";
			c_scaling_ = rowScaling;
/*
			if(get_trace_out().get())
				trace_out() << "\nEpetraNPFO::initialize(epetra_np): autoScaleStateConstraints==true: setting c_scaling = sqrt(rowScaling) ...\n";
			c_scaling_ = rowScaling;
			sqrt(&*c_scaling_->epetra_vec());
			if(get_trace_out().get()) {
				double min[1];
				TEST_FOR_EXCEPTION(0!=c_scaling_->epetra_vec()->MinValue(min),std::runtime_error,"Error!");
				trace_out()
					<< "\n  => time = " << timer.totalElapsedTime() << " sec\n"
					<< "\n  max{c_scaling} = " << norm_inf(*c_scaling_)
					<< "\n  min{c_scaling} = " << min[0]
					<< std::endl;
			}
*/
		}

		if( autoScaleStateVariables() ) {
			if(get_trace_out().get())
				trace_out() << "\nEpetraNPFO::initialize(epetra_np): autoScaleStateVariables==true: setting y_scaling = colScaling ...\n";
			y_scaling_ = colScaling;
/*
			if(get_trace_out().get())
				trace_out() << "\nEpetraNPFO::initialize(epetra_np): autoScaleStateVariables==true: setting y_scaling = sqrt(colScaling) ...\n";
			y_scaling_ = colScaling;
			sqrt(&*y_scaling_->epetra_vec());
			if(get_trace_out().get()) {
				double min[1];
				TEST_FOR_EXCEPTION(0!=y_scaling_->epetra_vec()->MinValue(min),std::runtime_error,"Error!");
				trace_out()
					<< "\n  => time = " << timer.totalElapsedTime() << " sec\n"
					<< "\n  max{y_scaling} = " << norm_inf(*y_scaling_)
					<< "\n  min{y_scaling} = " << min[0]
					<< std::endl;
			}
*/
		}

	}

	if( y_scaling_.get() ) {
		y_scaling_inv_ = rcp_dynamic_cast<EpetraVector>(space_y_->createMember());
		y_scaling_inv_->epetra_vec()->Reciprocal(*y_scaling_->epetra_vec());
	}

}

void EpetraNPFO::scale_y( const Epetra_MultiVector &y_unscaled, Epetra_MultiVector *y_scaled ) const
{
	if( y_scaling_.get() ) {
		y_scaled->Multiply( 1.0, *y_scaling_->epetra_vec(), y_unscaled, 0.0 );
	}
	else {
		if( y_scaled != &y_unscaled  )
			*y_scaled = y_unscaled;
	}
}

void EpetraNPFO::unscale_y( const Epetra_MultiVector &y_scaled, Epetra_MultiVector *y_unscaled ) const
{
	if( y_scaling_.get() ) {
		y_unscaled->ReciprocalMultiply( 1.0, *y_scaling_->epetra_vec(), y_scaled, 0.0 );
	}
	else {
		if( y_unscaled != &y_scaled  )
			*y_unscaled = y_scaled;
	}
}

void EpetraNPFO::scale_c( const Epetra_MultiVector &c_unscaled, Epetra_MultiVector *c_scaled ) const
{
	if( c_scaling_.get() ) {
		c_scaled->Multiply( 1.0, *c_scaling_->epetra_vec(), c_unscaled, 0.0 );
	}
	else {
		if( c_scaled != &c_unscaled )
			*c_scaled = c_unscaled;
	}
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
	using Teuchos::dyn_cast;
  const int Nu = epetra_np_->Nu();
	Teuchos::Time timer("");
  //
  // Get Epetra objects for y and u
  //
  const Epetra_Vector &y  = set_y(y_in);
  const Epetra_Vector **u = set_u( u_in, newPoint );
  //
  // Pick out the raw Epetra objects and in the process set their
  // adpater objects to uninitialized.
	//
	bool computeSomething = false;
  // c
  Teuchos::RefCountPtr<Epetra_Vector>  epetra_c;
  if(c_ && !c_updated_) {
		if(get_trace_out().get())
			trace_out()
				<< "\nEpetraNPFO::calc_Dc(...): Computing a Epetra_Vector object for c ... \n";
		c_->setUninitialized( &epetra_c );
		computeSomething = true;
	}
  // DcDy
  Teuchos::RefCountPtr<Epetra_Operator>  epetra_DcDy_op;
  Teuchos::RefCountPtr<Epetra_Operator>  epetra_DcDy_prec;
  const bool allow_specialized_DcDy_prec = true; // ToDo: Make this an external parameter!
  const bool specialized_DcDy_prec = ( epetra_np_->specialized_DcDy_prec() && allow_specialized_DcDy_prec );
  if(computeGradients && DcDy_ && !DcDy_updated_) {
		if( DcDy_->Op().get() && epetra_np_->DcDy_op_is_const() ) {
			if(get_trace_out().get())
				trace_out()
					<< "\nEpetraNPFO::calc_Dc(...): The DcDy is currently initialized and is constant so we will not recompute it ... \n";
		}
		else {
			if(get_trace_out().get())
				trace_out()
					<< "\nEpetraNPFO::calc_Dc(...): Computing a new Epetra_Operator object for DcDy ... \n";
			DcDy_->setUninitialized( &epetra_DcDy_op, NULL, NULL, &epetra_DcDy_prec, NULL );
			if( !epetra_DcDy_op.get() ) {
				epetra_DcDy_op = epetra_np_->create_DcDy_op();
				DcDy_->set_trace_out(get_trace_out());
			}
			if( !epetra_DcDy_prec.get() && specialized_DcDy_prec ) {
				epetra_DcDy_op = epetra_np_->create_DcDy_prec();
			}
			computeSomething = true;
		}
  }
  // DcDu
  for(int l=1;l<=Nu;++l) {
    if(computeGradients && DcDu_op_[l-1] && !DcDu_updated_[l-1]) {
			if( DcDu_op_[l-1]->epetra_op().get() && epetra_np_->DcDu_is_const(l) ) {
				if(get_trace_out().get())
					trace_out()
						<< "\nEpetraNPFO::calc_Dc(...): The DcDu_op("<<l<<") is currently initialized and is constant so we will not recompute it ... \n";
			}
			else {
				if(get_trace_out().get())
					trace_out()
						<< "\nEpetraNPFO::calc_Dc(...): Computing a Epetra_Operator object for DcDu("<<l<<") ... \n";
				DcDu_op_[l-1]->setUninitialized( &epetra_DcDu_op_[l-1] ); // Grab the RCP<Epetra_Operator>
				if(!epetra_DcDu_op_[l-1].get())
					epetra_DcDu_op_[l-1] = epetra_np_->create_DcDu_op(l);
				epetra_DcDu_args_[l-1] = Epetra::EpetraOp_or_EpetraMV(&*epetra_DcDu_op_[l-1]);
				computeSomething = true;
			}
    }
    else if(computeGradients && DcDu_mv_[l-1] && !DcDu_updated_[l-1]) {
			if( DcDu_mv_[l-1]->epetra_multi_vec().get() && epetra_np_->DcDu_is_const(l) ) {
				if(get_trace_out().get())
					trace_out()
						<< "\nEpetraNPFO::calc_Dc(...): The DcDu_mv("<<l<<") is currently initialized and is constant so we will not recompute it ... \n";
			}
			else {
				if(get_trace_out().get())
					trace_out()
						<< "\nEpetraNPFO::calc_Dc(...): Computing a DcDu("<<l<<") Epetra_MultiVector object ... \n";
				DcDu_mv_[l-1]->setUninitialized( &epetra_DcDu_mv_[l-1] ); // Grap the RCP<Epetra_MultiVector>
				if(!epetra_DcDu_mv_[l-1].get())
					epetra_DcDu_mv_[l-1] = epetra_np_->create_DcDu_mv(l);
				epetra_DcDu_args_[l-1] = Epetra::EpetraOp_or_EpetraMV(&*epetra_DcDu_mv_[l-1]);
				computeSomething = true;
			}
		}
		else {
      epetra_DcDu_args_[l-1] = Epetra::EpetraOp_or_EpetraMV(); // No DcDu(l) to compute!
    }
  }

	if(computeSomething) {
		//
		// Compute the set Epetra objects
		//
		if(get_trace_out().get())
			trace_out()
				<< "\nEpetraNPFO::calc_Dc(...): Calling " << typeid(*epetra_np_).name()
				<< ".calc_Dc(...) ...\n";
		timer.start(true);
		epetra_np_->calc_Dc(
			y
			,u
			,epetra_c.get()
			,epetra_DcDy_op.get()
			,specialized_DcDy_prec ? &*epetra_DcDy_prec : NULL
			,&epetra_DcDu_args_[0]
			);
		timer.stop();
		if(get_trace_out().get())
			trace_out()
				<< "\n  => time = " << timer.totalElapsedTime() << " sec\n";
		//
		// Put the raw Epetra objects back into the adpater objects
		//
		// c
		if( epetra_c.get() ) {
			if( c_scaling_.get() )
				scale_c( *epetra_c, &*epetra_c );
			c_->initialize( epetra_c, space_c_ );
			c_updated_ = true;
		}
		// DcDy
		if( epetra_DcDy_op.get() ) {
			// Test DcDy
			if(testOperators()) {
				if(get_trace_out().get())
					trace_out()
						<< "\nEpetraNPFO::calc_Dc(...): Testing the \'" << typeid(*epetra_DcDy_op).name()
						<< "\' object for DcDy ...\n";
				bool result = linearOpTester_.check(EpetraLinearOp(epetra_DcDy_op),get_trace_out().get());
				TEST_FOR_EXCEPTION(
					!result, std::runtime_error
					,"EpetraNPFO::calc_Dc(...): Error, test of DcDy operator object failed!"
					);
			}
			// Scale DcDy (before creating the preconditioner)
			if( c_scaling_.get() || y_scaling_.get() ) {
				Epetra_RowMatrix &epetra_DcDy_rm = dyn_cast<Epetra_RowMatrix>(*epetra_DcDy_op);
				if(c_scaling_.get())
					TEST_FOR_EXCEPTION(0!=epetra_DcDy_rm.LeftScale(*c_scaling_->epetra_vec()),std::logic_error,"Error!");
				if(y_scaling_.get())
					TEST_FOR_EXCEPTION(0!=epetra_DcDy_rm.RightScale(*y_scaling_inv_->epetra_vec()),std::logic_error,"Error!");
			}
			// Setup the matrix-based preconditioner
			if( usePrec() && !specialized_DcDy_prec ) {
				if(get_trace_out().get())
					trace_out()
						<< "\nEpetraNPFO::calc_Dc(...): Generating a preconditioner matrix for DcDy ... \n";
				precGenerator().setupPrec(epetra_DcDy_op,&epetra_DcDy_prec,get_trace_out().get());
				// Test preconditioner for DcDy
				if(testOperators()) {
					if(get_trace_out().get())
						trace_out()
							<< "\nEpetraNPFO::calc_Dc(...): Testing the \'" << typeid(*epetra_DcDy_prec).name()
							<< "\' object for the preconditioner for DcDy ...\n";
					bool result = linearOpTester_.check(EpetraLinearOp(epetra_DcDy_prec),get_trace_out().get());
					TEST_FOR_EXCEPTION(
						!result, std::runtime_error
						,"EpetraNPFO::calc_Dc(...): Error, test of DcDy preconditioner object failed!"
						);
				}
			}
			// Set the options for AztecOO::Iterate(...)
			DcDy_->maxIter(maxLinSolveIter());
			DcDy_->relTol(relLinSolveTol());
			DcDy_->minRelTol(1.0); // ToDo: Make this an external parameter
			// Set the linear system scaling options
			DcDy_->linearSystemScaler() = linearSystemScaler();
			// Finally initialize the aggregate LinearOpWithSolve object
			DcDy_->initialize(
				epetra_DcDy_op                                               // Op
				,epetra_np_->opDcDy() == Epetra::NOTRANS ? NOTRANS : TRANS   // Op_trans
				,Teuchos::rcp(&aztecOO_,false)                               // solver
				,epetra_DcDy_prec                                            // Prec
				,epetra_np_->opDcDy() == Epetra::NOTRANS ? NOTRANS : TRANS   // Prec_trans
				,Epetra::ProductOperator::APPLY_MODE_APPLY_INVERSE           // Prec_inverse (This is only for Ifpack!)
				,epetra_np_->adjointSupported()                              // adjointSupported
				);
			DcDy_updated_ = true;
		}
		// DcDu
		for(int l=1;l<=Nu;++l) {
			if(epetra_DcDu_op_[l-1].get()) {
				Teuchos::RefCountPtr<Epetra_Operator> &epetra_DcDu_op_l = epetra_DcDu_op_[l-1];
				if( c_scaling_.get() ) {
					Epetra_RowMatrix &epetra_DcDu_rm_l = dyn_cast<Epetra_RowMatrix>(*epetra_DcDu_op_l);
					TEST_FOR_EXCEPTION(0!=epetra_DcDu_rm_l.LeftScale(*c_scaling_->epetra_vec()),std::logic_error,"Error!");
				}
				DcDu_op_[l-1]->initialize(epetra_DcDu_op_l);
				if(testOperators()) {
					if(get_trace_out().get())
						trace_out()
							<< "\nEpetraNPFO::calc_Dc(...): Testing the \'" << typeid(*epetra_DcDu_op_l).name()
							<< "\' object for DcDu("<<l<<") ...\n";
					bool result = linearOpTester_.check(*DcDu_op_[l-1],get_trace_out().get());
					TEST_FOR_EXCEPTION(
						!result, std::runtime_error
						,"EpetraNPFO::calc_Dc(...): Error, test of DcDu("<<l<<") operator object failed!"
						);
				}
				epetra_DcDu_op_l = Teuchos::null;
				DcDu_updated_[l-1] = true;
			}
			else if(epetra_DcDu_mv_[l-1].get()) {
				Teuchos::RefCountPtr<Epetra_MultiVector> &epetra_DcDu_mv_l = epetra_DcDu_mv_[l-1];
				if( c_scaling_.get() ) {
					scale_c( *epetra_DcDu_mv_l, &*epetra_DcDu_mv_l );
				}
				DcDu_mv_[l-1]->initialize(epetra_DcDu_mv_l,space_c_);
				if(testOperators()) {
					if(get_trace_out().get())
						trace_out()
							<< "\nEpetraNPFO::calc_Dc(...): Testing the \'" << typeid(*epetra_DcDu_mv_l).name()
							<< "\' object for DcDu("<<l<<") ...\n";
					bool result = linearOpTester_.check(*DcDu_mv_[l-1],get_trace_out().get());
					TEST_FOR_EXCEPTION(
						!result, std::runtime_error
						,"EpetraNPFO::calc_Dc(...): Error, test of DcDu("<<l<<") operator object failed!"
						);
				}
				epetra_DcDu_mv_l = Teuchos::null;
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
	Teuchos::Time timer("");
  //
  // Get Epetra objects for y and u (unscaled)
  //
  const Epetra_Vector &y  = set_y(y_in);
  const Epetra_Vector **u = set_u( u_in, newPoint );
  //
  // Pick out the raw Epetra objects and in the process set their
  // adpater objects to uninitialized.
  //
	bool computeSomething = false;
  // g
  Teuchos::RefCountPtr<Epetra_Vector>  epetra_g;
  if(g_ && !g_updated_ ) {
		if(get_trace_out().get())
			trace_out()
				<< "\nEpetraNPFO::calc_Dg(...): Computing a Epetra_Vector object for g ... \n";
		g_->setUninitialized( &epetra_g );
		computeSomething = true;
	}
  // DgDy
  Teuchos::RefCountPtr<Epetra_MultiVector>  epetra_DgDy;
  if(computeGradients && DgDy_ && !DgDy_updated_) {
		if(get_trace_out().get())
			trace_out()
				<< "\nEpetraNPFO::calc_Dg(...): Computing a Epetra_MultiVector object for DgDy ... \n";
    DgDy_->setUninitialized( &epetra_DgDy );
    if( !epetra_DgDy.get() ) {
      epetra_DgDy = Teuchos::rcp(new Epetra_MultiVector(*epetra_np_->map_y(),space_g_->dim()));
    }
		computeSomething = true;
  }
  // DgDu
  for(int l=1;l<=Nu;++l) {
    if(computeGradients && DgDu_[l-1] && !DgDu_updated_[l-1]) {
			if(get_trace_out().get())
				trace_out()
					<< "\nEpetraNPFO::calc_Dg(...): Computing a Epetra_MultiVector object for DgDu("<<l<<") ... \n";
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
		if(get_trace_out().get())
			trace_out()
				<< "\nEpetraNPFO::calc_Dg(...): Calling " << typeid(*epetra_np_).name()
				<< ".calc_Dg(...) ...\n";
		timer.start(true);
		epetra_np_->calc_Dg(
			y
			,u
			,epetra_g.get()
			,epetra_DgDy.get()
			,&epetra_DgDu_args_[0]
			);
		timer.stop();
		if(get_trace_out().get())
			trace_out()
				<< "\n  => time = " << timer.totalElapsedTime() << " sec\n";
		
		//
		// Put the raw Epetra objects back into the adpater objects
		//
		if( epetra_g.get() ) {
			g_->initialize( epetra_g, space_g_ );
			g_updated_ = true;
		}
		// DgDy
		if(epetra_DgDy.get()) {
			if( y_scaling_.get() )
				unscale_y( *epetra_DgDy, &*epetra_DgDy );
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
