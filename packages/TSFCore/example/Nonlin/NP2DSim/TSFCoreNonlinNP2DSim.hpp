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

// ////////////////////////////////////////////////////////////
// TSFCoreNonlinNP2DSim.hpp

//#define TSFCORE_NONLIN_NP_2D_SIM_USE_SIMPLE_GMRES_SOLVER

#ifndef TSFCORE_NONLIN_NP_2D_SIM_HPP
#define TSFCORE_NONLIN_NP_2D_SIM_HPP

#include "TSFCoreNonlinNP2DSimDecl.hpp"
#include "TSFCoreSerialVectorSpace.hpp"
#include "TSFCoreMultiVector.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "TSFCoreExplicitVectorView.hpp"
#include "TSFCoreExplicitMultiVectorView.hpp"
#include "TSFCoreSolversBiCGSolver.hpp"
#include "TSFCoreSolversNormedConvergenceTester.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_TestForException.hpp"
#include "AbstractFactoryStd.hpp"

#ifdef TSFCORE_NONLIN_NP_2D_SIM_USE_SIMPLE_GMRES_SOLVER
#include "TSFCoreSolversSimpleGMRESSolver.hpp"
#else
#include "TSFCoreSolversGMRESSolver.hpp"
#endif


namespace TSFCore {
namespace Nonlin {

template<class Scalar>
NP2DSim<Scalar>::NP2DSim(
	const Scalar                                             a
	,const Scalar                                            b
	,const Scalar                                            d
	,const Scalar                                            lin_sol_tol
	,const ELinearSolverType                                 linearSolverType
	,const bool                                              dumpToFile
	,const Teuchos::RefCountPtr<const VectorSpace<Scalar> >  &space_y_c
	)
	:linearSolverType_(linearSolverType),dumpToFile_(dumpToFile),isInitialized_(false)
	,a_(a),b_(b),d_(d)
	,lin_sol_tol_(lin_sol_tol),c_(NULL),DcDy_(NULL)
{
	if(space_y_c.get()) {
		TEST_FOR_EXCEPTION( !(space_y_c->dim()==2), std::logic_error, "NP2DSim<Scalar>::NP2DSim(...): Error!" );
		space_y_c_ = space_y_c;
	}
	else {
		space_y_c_ = Teuchos::rcp(new SerialVectorSpace<Scalar>(2));
	}
	factory_DcDy_ = Teuchos::rcp(new MemMngPack::AbstractFactoryStd<LinearOpWithSolve<Scalar>,LinearOpWithSolveIter<Scalar> >());
}

template<class Scalar>
void NP2DSim<Scalar>::set_y0( const Scalar y01, const Scalar y02 )
{
	set_ele( 1, y01, &*y0_ );
	set_ele( 2, y02, &*y0_ );
}

// Overridden from NonlinearProblem

template<class Scalar>
void NP2DSim<Scalar>::initialize( bool testSetup )
{
	if(isInitialized_) return;
	const Scalar inf_bnd = NonlinearProblem<Scalar>::infiniteBound();
	yL_ = space_y_c_->createMember(); assign(yL_.get(),-inf_bnd);
	yU_ = space_y_c_->createMember(); assign(yU_.get(),+inf_bnd);
	y0_ = space_y_c_->createMember(); assign(y0_.get(),0.0);
	isInitialized_ = true;
}

template<class Scalar>
bool NP2DSim<Scalar>::isInitialized() const
{
	return isInitialized_;
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpace<Scalar> >
NP2DSim<Scalar>::space_y() const
{
	return space_y_c_;
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpace<Scalar> >
NP2DSim<Scalar>::space_c() const
{
	return space_y_c_;
}

template<class Scalar>
const Vector<Scalar>& NP2DSim<Scalar>::yL() const
{
	return *yL_;
}

template<class Scalar>
const Vector<Scalar>& NP2DSim<Scalar>::yU() const
{
	return *yU_;
}

template<class Scalar>
const Vector<Scalar>& NP2DSim<Scalar>::y0() const
{
	return *y0_;
}

template<class Scalar>
void NP2DSim<Scalar>::set_c(Vector<Scalar>* c)
{
#ifdef _DEBUG
	if(c) {
		TEST_FOR_EXCEPTION(
			!c->space()->isCompatible(*this->space_c()), Exceptions::IncompatibleVectorSpaces
		, "NP2DSim<Scalar>::set_c(...): Error!" );
	}
#endif
	c_ = c;
}

template<class Scalar>
Vector<Scalar>* NP2DSim<Scalar>::get_c()
{
	return c_;
}

template<class Scalar>
void NP2DSim<Scalar>::unsetQuantities()
{
	c_    = NULL;
	DcDy_ = NULL;
}

template<class Scalar>
void NP2DSim<Scalar>::calc_c(
	const Vector<Scalar>     &y_in
	,const Vector<Scalar>*   u_in[]
	,bool                    newPoint
	) const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		!y_in.space()->isCompatible(*this->space_y()), Exceptions::IncompatibleVectorSpaces
		, "NP2DSim<Scalar>::calc_c(...): Error!" );
	TEST_FOR_EXCEPTION(
		u_in != NULL, std::invalid_argument
		, "NP2DSim<Scalar>::calc_c(...): Error!" );
	TEST_FOR_EXCEPTION(
		c_ == NULL, std::logic_error
		, "NP2DSim<Scalar>::calc_c(...): Error!" );
#endif
	// Get at the data
	ExplicitVectorView<Scalar>        y(y_in);
	ExplicitMutableVectorView<Scalar> c(*c_);
	// Compute c(y)
	c(1) =        y(1)      + y(2)*y(2) - a_  ;
	c(2) = d_ * ( y(1)*y(1) - y(2)      - b_ );
}

// Overridden from NonlinearProblemFirstOrder

template<class Scalar>
Teuchos::RefCountPtr< const MemMngPack::AbstractFactory<LinearOpWithSolve<Scalar> > >
NP2DSim<Scalar>::factory_DcDy() const
{
	return factory_DcDy_;
}

template<class Scalar>
ETransp NP2DSim<Scalar>::opDcDy() const
{
	return NOTRANS;
}

template<class Scalar>
void NP2DSim<Scalar>::set_DcDy(LinearOpWithSolve<Scalar>* DcDy)
{
	using Teuchos::dyn_cast;
	if(DcDy)  DcDy_ = &dyn_cast<LinearOpWithSolveIter<Scalar> >(*DcDy);
	else      DcDy_ = NULL;
}

template<class Scalar>
LinearOpWithSolve<Scalar>* NP2DSim<Scalar>::get_DcDy()
{
	return DcDy_;
}

template<class Scalar>
void NP2DSim<Scalar>::calc_DcDy(
	const Vector<Scalar>     &y_in
	,const Vector<Scalar>*   u_in[]
	,bool                    newPoint
	) const
{
	// Validate preconditions
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		!y_in.space()->isCompatible(*this->space_y()), Exceptions::IncompatibleVectorSpaces
		, "NP2DSim<Scalar>::calc_DcDy(...): Error!" );
	TEST_FOR_EXCEPTION(
		u_in != NULL, std::invalid_argument
		, "NP2DSim<Scalar>::calc_DcDy(...): Error!" );
	TEST_FOR_EXCEPTION(
		DcDy_ == NULL, std::logic_error
		, "NP2DSim<Scalar>::calc_DcDy(...): Error!" );
#endif
	// Get the underlying multi-vector operator, iterative linear
	// solver, convergence tester and preconditioner
	LinearOpWithSolveIterState<Scalar>
		DcDy_state = DcDy_->setUninitialized();
	Teuchos::RefCountPtr<MultiVector<Scalar> >
		DcDy_mv = Teuchos::rcp_dynamic_cast<MultiVector<Scalar> >(
			Teuchos::rcp_const_cast<LinearOp<Scalar> >(DcDy_state.M)
			);
	Teuchos::RefCountPtr<const Solvers::IterativeLinearSolver<Scalar> >
		DcDy_solver = DcDy_state.solver;
	Teuchos::RefCountPtr<Solvers::ConvergenceTester<Scalar> >
		convTester = DcDy_state.convTester;
	Teuchos::RefCountPtr<MultiVector<Scalar> >
		DcDy_prec = Teuchos::rcp_dynamic_cast<MultiVector<Scalar> >(
			Teuchos::rcp_const_cast<LinearOp<Scalar> >(DcDy_state.M_tilde_left_inv)
			);
	// Create these objects if they have not been created already
	if( !DcDy_mv.get() ) DcDy_mv = this->space_y()->createMembers(2);
	if( !DcDy_solver.get() ) {
		// Create output stream
		if( dumpToFile() && !linear_solver_out_.get() )
			linear_solver_out_ = Teuchos::rcp(new std::ofstream("NP2DSim::linearSolver.out"));
		// Create the linear solver
		switch(linearSolverType()) {
			case LINEAR_SOLVER_BICG: {
				Teuchos::RefCountPtr<Solvers::BiCGSolver<Scalar> >
					bicg_solver = Teuchos::rcp(new Solvers::BiCGSolver<Scalar>());
				if(dumpToFile()) {
					bicg_solver->set_out(linear_solver_out_);
					bicg_solver->dump_all(true);
				}
				DcDy_solver = bicg_solver;
				break;
			}
			case LINEAR_SOLVER_GMRES: {
#ifdef TSFCORE_NONLIN_NP_2D_SIM_USE_SIMPLE_GMRES_SOLVER
				Teuchos::RefCountPtr<Solvers::SimpleGMRESSolver<Scalar> >
					gmres_solver = Teuchos::rcp(new Solvers::SimpleGMRESSolver<Scalar>());
#else
				Teuchos::RefCountPtr<Solvers::GMRESSolver<Scalar> >
					gmres_solver = Teuchos::rcp(new Solvers::GMRESSolver<Scalar>());
#endif
				// ToDo: Put in this type of printing?
				//if(dumpToFile()) {
				//	gmres_solver->set_out(linear_solver_out);
				//	gmres_solver->dump_all(true);
				//}
				DcDy_solver = gmres_solver;
				break;
			}
		}
  }
	if( !convTester.get() ) convTester = Teuchos::rcp(new Solvers::NormedConvergenceTester<Scalar>(lin_sol_tol_));
	// Get at the state vector and Jacobian data
	if(true) {
		ExplicitVectorView<Scalar> y(y_in);
		ExplicitMutableMultiVectorView<Scalar> DcDy(*DcDy_mv);
		// Fill the Jacobian
		DcDy(1,1) = 1.0;           DcDy(1,2) = 2.0*y(2);
		DcDy(2,1) = 2.0*d_*y(1);   DcDy(2,2) = -d_;
	}
	// Initialize the precondtioner to a Diagonal matrix only the first time!
	if( !DcDy_prec.get() ) {
		if(linearSolverType()==LINEAR_SOLVER_BICG) {
			DcDy_prec  = this->space_y()->createMembers(2);
			ExplicitMutableMultiVectorView<Scalar> M_tilde(*DcDy_prec);
			// Fill the preconditioner \tilde{M}
			M_tilde(1,1) = 1.0;       M_tilde(1,2) = 0.0;
			M_tilde(2,1) = 0.0;       M_tilde(2,2) = -1.0/d_;
		}
		// GMRES can not handle preconditioners yet!
	}
	// Initialize the DcDy object given its initialized parts
	DcDy_->initialize(
		DcDy_mv, NOTRANS      // M, M_trans
		,DcDy_solver          // solver
		,convTester           // default convTester
		,DcDy_prec, NOTRANS   // M_tilde_left_inv, M_tilde_left_inv_trans
		);
}

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLIN_NP_2D_SIM_HPP
