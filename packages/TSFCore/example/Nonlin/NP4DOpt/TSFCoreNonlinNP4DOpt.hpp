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

// ///////////////////////////////////////////////////////////////////
// TSFCoreNonlinNP4DOpt.hpp

#ifndef TSFCORE_NONLIN_NP_4D_OPT_HPP
#define TSFCORE_NONLIN_NP_4D_OPT_HPP

#include "TSFCoreNonlinNP4DOptDecl.hpp"
#include "TSFCoreMultiVectorAllocator.hpp"
#include "TSFCoreExplicitVectorView.hpp"
#include "TSFCoreExplicitMultiVectorView.hpp"
#include "TSFCoreNonlinLinearOpWithSolve.hpp"
#include "Teuchos_TestForException.hpp"

// Debugging ony
//#include <iostream>
//#include "TSFExtended/src/Core/TSFCoreTestingTools.hpp"

#ifdef _DEBUG
#define NP4DOPT_VALIDATE_L_IN_RANGE(l) TEST_FOR_EXCEPTION( (l) != 1, std::invalid_argument, "Error!  l == " << l << " != 1!")
#else
#define NP4DOPT_VALIDATE_L_IN_RANGE(l)
#endif

namespace TSFCore {
namespace Nonlin {

// Constructors / Initializers / accessors

template<class Scalar>
NP4DOpt<Scalar>::NP4DOpt(
	const Scalar                                             yt1
	,const Scalar                                            yt2
	,const Scalar                                            ut1
	,const Scalar                                            ut2
	,const Scalar                                            d
	,const Scalar                                            lin_sol_tol
	,const typename NP2DSim<Scalar>::ELinearSolverType       linearSolverType
	,const bool                                              dumpToFile
	,const Teuchos::RefCountPtr<const VectorSpace<Scalar> >  &space_y_c
	)
	:np2dsim_(ut1,ut2,d,lin_sol_tol,linearSolverType,dumpToFile,space_y_c)
	,yt1_(yt1),yt2_(yt2),ut1_(ut1),ut2_(ut2)
{
	typedef LinearOp<Scalar>                LO;
	typedef MultiVector<Scalar>             MV;
	typedef MultiVectorAllocator<Scalar>    MVA;
	typedef MemMngPack::PostModNothing<MV>  PMN;
	space_g_      = np2dsim_.space_y()->smallVecSpcFcty()->createVecSpc(NUM_RESPONSE_FUNCTIONS);
	factory_DcDu_ = Teuchos::rcp(new MemMngPack::AbstractFactoryStd<LO,MV,PMN,MVA>(PMN(),MVA(np2dsim_.space_c(),2)));
}

template<class Scalar>
void NP4DOpt<Scalar>::set_u0( const Scalar u01, const Scalar u02 )
{
	TSFCore::set_ele( 1, u01, &*u0_ );
	TSFCore::set_ele( 2, u02, &*u0_ );
}

// Overridden from NonlinearProblem

template<class Scalar>
void NP4DOpt<Scalar>::initialize( bool testSetup )
{
	if(np2dsim_.isInitialized()) return;
	np2dsim_.initialize(testSetup);
	const VectorSpace<Scalar> &space_y_c = *np2dsim_.space_y();
	const Scalar inf_bnd = NonlinearProblem<Scalar>::infiniteBound();
	uL_ = space_y_c.createMember(); assign(uL_.get(),-inf_bnd);
	uU_ = space_y_c.createMember(); assign(uU_.get(),+inf_bnd);
	u0_ = space_y_c.createMember(); { ExplicitMutableVectorView<Scalar> u0(*u0_); u0(1)=ut1_; u0(2)=ut2_; }
	gL_ = space_g_->createMember(); assign(gL_.get(),-10.0);
	gU_ = space_g_->createMember(); assign(gU_.get(),+10.0);
}

template<class Scalar>
bool NP4DOpt<Scalar>::isInitialized() const
{
	return np2dsim_.isInitialized();
}

template<class Scalar>
int NP4DOpt<Scalar>::Nu() const
{
	return 1;
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpace<Scalar> >
NP4DOpt<Scalar>::space_y() const
{
	return np2dsim_.space_y();
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpace<Scalar> >
NP4DOpt<Scalar>::space_u(int l) const
{
	NP4DOPT_VALIDATE_L_IN_RANGE(l);
	return np2dsim_.space_y();  // This works!
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpace<Scalar> >
NP4DOpt<Scalar>::space_c() const
{
	return np2dsim_.space_c();
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpace<Scalar> >
NP4DOpt<Scalar>::space_g() const
{
	return space_g_;
}

template<class Scalar>
const Vector<Scalar>&
NP4DOpt<Scalar>::yL() const
{
	return np2dsim_.yL();
}

template<class Scalar>
const Vector<Scalar>&
NP4DOpt<Scalar>::yU() const
{
	return np2dsim_.yU();
}

template<class Scalar>
const Vector<Scalar>&
NP4DOpt<Scalar>::uL(int l) const
{
	NP4DOPT_VALIDATE_L_IN_RANGE(l);
	return *uL_;
}

template<class Scalar>
const Vector<Scalar>&
NP4DOpt<Scalar>::uU(int l) const
{
	NP4DOPT_VALIDATE_L_IN_RANGE(l);
	return *uU_;
}

template<class Scalar>
const Vector<Scalar>&
NP4DOpt<Scalar>::gL() const
{
	return *gL_;
}

template<class Scalar>
const Vector<Scalar>&
NP4DOpt<Scalar>::gU() const
{
	return *gU_;
}

template<class Scalar>
const Vector<Scalar>&
NP4DOpt<Scalar>::y0() const
{
	return np2dsim_.y0();
}

template<class Scalar>
const Vector<Scalar>&
NP4DOpt<Scalar>::u0(int l) const
{
	NP4DOPT_VALIDATE_L_IN_RANGE(l);
	return *u0_;
}

template<class Scalar>
void NP4DOpt<Scalar>::set_c(Vector<Scalar>* c)
{
	np2dsim_.set_c(c);
}

template<class Scalar>
Vector<Scalar>* NP4DOpt<Scalar>::get_c()
{
	return np2dsim_.get_c();
}

template<class Scalar>
void NP4DOpt<Scalar>::set_g(Vector<Scalar>* g)
{
#ifdef _DEBUG
	if(g) {
		TEST_FOR_EXCEPTION(
			!g->space()->isCompatible(*this->space_g()), Exceptions::IncompatibleVectorSpaces
		, "NP4DOpt<Scalar>::set_g(...): Error!" );
	}
#endif
	g_ = g;
}

template<class Scalar>
Vector<Scalar>* NP4DOpt<Scalar>::get_g()
{
	return g_;	
}

template<class Scalar>
void NP4DOpt<Scalar>::unsetQuantities()
{
	np2dsim_.unsetQuantities();
	g_ = NULL; DcDu_ = NULL; DgDy_ = NULL; DgDu_ = NULL;
}

template<class Scalar>
void NP4DOpt<Scalar>::calc_c(
	const Vector<Scalar>     &y
	,const Vector<Scalar>*   u[]
	,bool                    newPoint
	) const
{
	set_u(u,newPoint);
	np2dsim_.calc_c(y,NULL,newPoint);
}

template<class Scalar>
void NP4DOpt<Scalar>::calc_g(
	const Vector<Scalar>     &y_in
	,const Vector<Scalar>*   u_in[]
	,bool                    newPoint
	) const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( u_in==NULL, std::invalid_argument, "Error!" );
	TEST_FOR_EXCEPTION( !u_in[0]->space()->isCompatible(*this->space_u(1)), Exceptions::IncompatibleVectorSpaces, "Error!" );
	TEST_FOR_EXCEPTION( g_ == NULL, std::logic_error, "Error!" );
#endif
	// Get at the data
	ExplicitVectorView<Scalar>        y(y_in);
	ExplicitVectorView<Scalar>        u(*u_in[0]);
	ExplicitMutableVectorView<Scalar> g(*g_);
	// Compute g(u)
	g(1) = y(1) - yt1_;
	g(2) = y(2) - yt2_;
	g(3) = u(1) - ut1_;
	g(4) = u(2) - ut2_;
}

// Overridden from NonlinearProblemFirstOrder

template<class Scalar>
Teuchos::RefCountPtr< const MemMngPack::AbstractFactory<LinearOpWithSolve<Scalar> > >
NP4DOpt<Scalar>::factory_DcDy() const
{
	return np2dsim_.factory_DcDy();
}

template<class Scalar>
Teuchos::RefCountPtr< const MemMngPack::AbstractFactory<LinearOp<Scalar > > >
NP4DOpt<Scalar>::factory_DcDu(int l) const
{
	NP4DOPT_VALIDATE_L_IN_RANGE(l);
	return factory_DcDu_;
}

template<class Scalar>
ETransp NP4DOpt<Scalar>::opDcDy() const
{
	return np2dsim_.opDcDy();
}

template<class Scalar>
ETransp NP4DOpt<Scalar>::opDcDu(int l) const
{
	NP4DOPT_VALIDATE_L_IN_RANGE(l);
	return TRANS;
}

template<class Scalar>
void NP4DOpt<Scalar>::set_DcDy(LinearOpWithSolve<Scalar>* DcDy)
{
	np2dsim_.set_DcDy(DcDy);
}

template<class Scalar>
LinearOpWithSolve<Scalar>* NP4DOpt<Scalar>::get_DcDy()
{
	return np2dsim_.get_DcDy();
}

template<class Scalar>
void NP4DOpt<Scalar>::set_DcDu(int l, LinearOp<Scalar>* DcDu_l)
{
	NP4DOPT_VALIDATE_L_IN_RANGE(l);
	using Teuchos::dyn_cast;
	if(DcDu_l)  DcDu_ = &dyn_cast<MultiVector<Scalar> >(*DcDu_l);
	else        DcDu_ = NULL;
}

template<class Scalar>
LinearOp<Scalar>* NP4DOpt<Scalar>::get_DcDu(int l)
{
	NP4DOPT_VALIDATE_L_IN_RANGE(l);
	return DcDu_;
}

template<class Scalar>
void NP4DOpt<Scalar>::set_DgDy(MultiVector<Scalar>* DgDy)
{
#ifdef _DEBUG
	if(DgDy) {
		TEST_FOR_EXCEPTION(
			!DgDy->range()->isCompatible(*this->space_y())
			,Exceptions::IncompatibleVectorSpaces
			, "NP4DOpt<Scalar>::set_DgDy(...): Error!" );
		TEST_FOR_EXCEPTION(
			!(DgDy->domain()->dim() == NUM_RESPONSE_FUNCTIONS)
			,Exceptions::IncompatibleVectorSpaces
			, "NP4DOpt<Scalar>::set_DgDy(...): Error!" );
	}
#endif
	DgDy_ = DgDy;
}

template<class Scalar>
MultiVector<Scalar>* NP4DOpt<Scalar>::get_DgDy()
{
	return DgDy_;
}

template<class Scalar>
void NP4DOpt<Scalar>::set_DgDu(int l, MultiVector<Scalar>* DgDu_l)
{
	NP4DOPT_VALIDATE_L_IN_RANGE(l);
#ifdef _DEBUG
	if(DgDu_l) {
		TEST_FOR_EXCEPTION(
			!DgDu_l->range()->isCompatible(*this->space_u(1))
			,Exceptions::IncompatibleVectorSpaces
			, "NP4DOpt<Scalar>::set_DgDu(...): Error!" );
		TEST_FOR_EXCEPTION(
			!(DgDu_l->domain()->dim() == NUM_RESPONSE_FUNCTIONS)
			,Exceptions::IncompatibleVectorSpaces
			, "NP4DOpt<Scalar>::set_DgDu(...): Error!" );
	}
#endif
	DgDu_ = DgDu_l;
}

template<class Scalar>
MultiVector<Scalar>* NP4DOpt<Scalar>::get_DgDu(int l)
{
	NP4DOPT_VALIDATE_L_IN_RANGE(l);
	return DgDu_;
}

template<class Scalar>
void NP4DOpt<Scalar>::calc_DcDy(
	const Vector<Scalar>     &y
	,const Vector<Scalar>*   u[]
	,bool                    newPoint
	) const
{
	set_u(u,newPoint);
	np2dsim_.calc_DcDy(y,NULL,newPoint);
}

template<class Scalar>
void NP4DOpt<Scalar>::calc_DcDu(
	int                      l
	,const Vector<Scalar>    &y_in
	,const Vector<Scalar>*   u_in[]
	,bool                    newPoint
	) const
{
	NP4DOPT_VALIDATE_L_IN_RANGE(l);
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( u_in==NULL, std::invalid_argument, "Error!" );
	TEST_FOR_EXCEPTION( !u_in[0]->space()->isCompatible(*this->space_u(1)), Exceptions::IncompatibleVectorSpaces, "Error!" );
	TEST_FOR_EXCEPTION( DcDu_ == NULL, std::logic_error, "Error!" );
#endif
	ExplicitVectorView<Scalar> u(*u_in[0]);
	ExplicitMutableMultiVectorView<Scalar> DcDu_t(*DcDu_);
	// Fill the Jacobian DcDu'
	DcDu_t(1,1) = -1.0;   DcDu_t(1,2) = 0.0;
	DcDu_t(2,1) = 0.0;    DcDu_t(2,2) = -np2dsim_.get_d();
}

template<class Scalar>
void NP4DOpt<Scalar>::calc_DgDy(
	const Vector<Scalar>     &y_in
	,const Vector<Scalar>*   u_in[]
	,bool                    newPoint
	) const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( DgDy_ == NULL, std::logic_error, "Error!" );
#endif
	assign( DgDy_, 0.0 );
	ExplicitVectorView<Scalar> y(y_in);
	ExplicitMutableMultiVectorView<Scalar> DgDy_t(*DgDy_);
	// Fill nonzeros in DgDy'
	DgDy_t(1,1) = 1.0; DgDy_t(2,2) = 1.0;
}

template<class Scalar>
void NP4DOpt<Scalar>::calc_DgDu(
	int                      l
	,const Vector<Scalar>    &y_in
	,const Vector<Scalar>*   u_in[]
	,bool                    newPoint
	) const
{
	NP4DOPT_VALIDATE_L_IN_RANGE(l);
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( u_in==NULL, std::invalid_argument, "Error!" );
	TEST_FOR_EXCEPTION( !u_in[0]->space()->isCompatible(*this->space_u(1)), Exceptions::IncompatibleVectorSpaces, "Error!" );
	TEST_FOR_EXCEPTION( DgDu_ == NULL, std::logic_error, "Error!" );
#endif
	assign( DgDu_, 0.0 );
	ExplicitVectorView<Scalar> u(*u_in[0]);
	ExplicitMutableMultiVectorView<Scalar> DgDu_t(*DgDu_);
	// Fill nonzeros in DgDu'
	DgDu_t(1,3) = 1.0;  DgDu_t(2,4) = 1.0;
}

// private

template<class Scalar>
void NP4DOpt<Scalar>::set_u(const Vector<Scalar>* u_in[], bool newPoint) const
{
#ifdef _DEBUG
	if(u_in) TEST_FOR_EXCEPTION( !u_in[0]->space()->isCompatible(*this->space_u(1)), Exceptions::IncompatibleVectorSpaces, "Error!" );
#endif
	if(!newPoint) return;
	if(u_in) {
		ExplicitVectorView<Scalar> u(*u_in[0]);
		const_cast<NP2DSim<Scalar>&>(np2dsim_).set_a(u(1));
		const_cast<NP2DSim<Scalar>&>(np2dsim_).set_b(u(2));
	}
	else {
		ExplicitVectorView<Scalar> u(*u0_);
		const_cast<NP2DSim<Scalar>&>(np2dsim_).set_a(u(1));
		const_cast<NP2DSim<Scalar>&>(np2dsim_).set_b(u(2));
	}
}

} // namespace Nonlin
} // namespace TSFCore

#undef NP4DOPT_VALIDATE_L_IN_RANGE

#endif // TSFCORE_NONLIN_NP_4D_OPT_HPP
