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

// //////////////////////////////////////////////////////////////////
// TSFCoreNonlinNonlinearProblemFirstOrderTester.hpp

#ifndef TSFCORE_NONLIN_NONLINEAR_PROBLEM_FIRST_ORDER_TESTER_HPP
#define TSFCORE_NONLIN_NONLINEAR_PROBLEM_FIRST_ORDER_TESTER_HPP

#include "TSFCoreNonlinNonlinearProblemFirstOrderTesterDecl.hpp"
#include "TSFCoreTestingTools.hpp"
#include "TSFCoreNonlinNonlinearProblem.hpp"
#include "TSFCoreNonlinNonlinearProblemFirstOrder.hpp"
#include "TSFCoreNonlinLinearOpWithSolve.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVector.hpp"
#include "Teuchos_AbstractFactory.hpp"

namespace TSFCore {
namespace Nonlin {

template<class Scalar>
bool
NonlinearProblemFirstOrderTester<Scalar>::doTest(
	NonlinearProblemFirstOrder<Scalar>* prob, std::ostream* o, EOutputLevel olevel
	) const
{
	namespace mmp = MemMngPack;
	using std::endl;

	if(o) *o << "Entering NonlinearProblemFirstOrderTester<Scalar>::do_test(prob,...)\n\n";
	// Removed std::boolalpha from output stream because it doesn't exist for the GNU 2.96 compiler (HKT, 09/22/03)
	// if(o) *o << "Entering NonlinearProblemFirstOrderTester<Scalar>::do_test(prob,...)\n\n" << std::boolalpha;

	bool result, success = true;
	int l;

	// Basic info

	const int Nu = prob->Nu();
	const int nrf = prob->numResponseFunctions();

	if(o) *o << "prob->adjointSupported() -> "; const bool adjointSupported = prob->adjointSupported(); if(o) *o << adjointSupported << endl;

	// Factories

	typedef Teuchos::RefCountPtr< const Teuchos::AbstractFactory<LinearOpWithSolve<Scalar> > >   DcDy_fcty_ptr_t;
	typedef Teuchos::RefCountPtr< const Teuchos::AbstractFactory<LinearOp<Scalar > > >           DcDu_fcty_ptr_t;

	if(o) *o << "prob->factory_DcDy().get() -> ";       const DcDy_fcty_ptr_t DcDy_fcty = prob->factory_DcDy(); if(o) *o << DcDy_fcty.get() << endl;
	if(o) *o << "prob->factory_DcDy().get()!=NULL -> "; if(!(result=(DcDy_fcty.get()!=NULL))) success=false;    if(o) *o << passfail(result) << endl;

	for( l = 1; l <= Nu; ++l ) {
		if(o) *o << "prob->factory_DcDu("<<l<<").get() -> ";       const DcDu_fcty_ptr_t DcDy_l_fcty = prob->factory_DcDu(l); if(o) *o << DcDy_l_fcty.get() << endl;
		if(o) *o << "prob->factory_DcDu("<<l<<").get()!=NULL -> "; if(!(result=(DcDy_l_fcty.get()!=NULL))) success = false;   if(o) *o << passfail(result) << endl;
	}

	if(o) *o << "prob->opDcDy() -> "; const ETransp opDcDy = prob->opDcDy(); if(o) *o << toString(opDcDy) << endl;
	for( l = 1; l <= Nu; ++l ) {
		if(o) *o << "prob->opDcDu("<<l<<") -> "; const ETransp opDcDu_l = prob->opDcDu(l); if(o) *o << toString(opDcDu_l) << endl;
	}

	// Calulations

	if(o) *o << "DcDy = prob->factory_DcDy()->create();\n";
	Teuchos::RefCountPtr<LinearOpWithSolve<Scalar> > DcDy = DcDy_fcty->create();
	if(o) *o << "prob->set_DcDy(DcDy)\n"; prob->set_DcDy(DcDy.get());
	if(o) *o << "prob->get_DcDy() -> ";
	LinearOpWithSolve<Scalar> *get_DcDy = prob->get_DcDy();
	if(o) *o << get_DcDy << endl;

	std::vector<Teuchos::RefCountPtr<LinearOp<Scalar> > >  DcDu(Nu);
	for( l = 1; l <= Nu; ++l ) {
		if(o) *o << "DcDu("<<l<<") = prob->factory_DcDu(l)->create();\n";
		DcDu[l-l] = prob->factory_DcDu(l)->create();
		if(o) *o << "prob->set_DcDu("<<l<<",DcDu("<<l<<"))\n";  prob->set_DcDu(l,DcDu[l-1].get());
		if(o) *o << "prob->get_DcDu("<<l<<") -> "; LinearOp<Scalar> *get_DcDu_l = prob->get_DcDu(l); if(o) *o << get_DcDu_l << endl;
	}

	Teuchos::RefCountPtr<MultiVector<Scalar> >                DgDy;
	std::vector<Teuchos::RefCountPtr<MultiVector<Scalar> > >  DgDu(Nu);
	if(nrf) {
		if(o) *o << "DgDy = prob->space_y()->createMembers(prob->space_g()->dim());\n";
		DgDy = prob->space_y()->createMembers(nrf);
		if(o) *o << "prob->set_DgDy(DgDy)\n"; prob->set_DgDy(DgDy.get());
		if(o) *o << "prob->get_DgDy() -> "; MultiVector<Scalar> *get_DgDy = prob->get_DgDy(); if(o) *o << get_DgDy << endl;
		
		for( l = 1; l <= Nu; ++l ) {
			if(o) *o << "DgDu("<<l<<") = prob->space_u("<<l<<")->createMembers(prob->space_g()->dim());\n";
			DgDu[l-l] = prob->space_u(l)->createMembers(nrf);
			if(o) *o << "prob->set_DgDu("<<l<<",DgDu("<<l<<"))\n";  prob->set_DgDu(l,DgDu[l-1].get());
			if(o) *o << "prob->get_DgDu("<<l<<") -> "; MultiVector<Scalar> *get_DgDu_l = prob->get_DgDu(l); if(o) *o << get_DgDu_l << endl;
		}
	}

	const Vector<Scalar>                &y0 = prob->y0();
	std::vector<const Vector<Scalar>*>  u0_store(Nu); for(l=1;l<=Nu;++l) u0_store[l-1] = &prob->u0(l);
	const Vector<Scalar>                **u0 = NULL; if(Nu) u0 = &u0_store[0];
	
	if(o) *o << "prob->calc_DcDy(y0,u0,true);\n"; prob->calc_DcDy(y0,u0,true);
	if(o && olevel>=VERBOSE_OUTPUT) *o << "*prob->get_DcDy() = \n" << *prob->get_DcDy();

	for( l = 1; l <= Nu; ++l ) {
		if(o) *o << "prob->calc_DcDu("<<l<<",y0,u0,true);\n"; prob->calc_DcDu(l,y0,u0,true);
		if(o && olevel>=VERBOSE_OUTPUT) *o << "*prob->get_DcDu("<<l<<") = \n" << *prob->get_DcDu(l);
	}
	
	if(nrf) {
		if(o) *o << "prob->calc_DgDy(y0,u0,true);\n"; prob->calc_DgDy(y0,u0,true);
		if(o && olevel>=VERBOSE_OUTPUT) *o << "*prob->get_DgDy() = \n" << *prob->get_DgDy();

		for( l = 1; l <= Nu; ++l ) {
			if(o) *o << "prob->calc_DgDu("<<l<<",y0,u0,true);\n"; prob->calc_DgDu(l,y0,u0,true);
			if(o && olevel>=VERBOSE_OUTPUT) *o << "*prob->get_DgDu("<<l<<") = \n" << *prob->get_DgDu(l);
		}
	}

	// Unset quantities
	
	if(o) *o << "prob->unsetQuantities();\n"; prob->unsetQuantities();
	if(o) *o << "prob->get_c() -> "; Vector<Scalar> *get_c = prob->get_c(); if(o) *o << get_c << endl;
	if(nrf) {
		if(o) *o << "prob->get_g() -> "; Vector<Scalar> *get_g = prob->get_g(); if(o) *o << get_g << endl;
	}
	// ToDo: add calls to get_DcDy() ... etc.

	// Test DcDy

	if(true) {
		// Test NONTRANS
		Teuchos::RefCountPtr<Vector<Scalar> >
			t1 = DcDy->range()->createMember(),
			t2 = DcDy->domain()->createMember(),
			t3 = DcDy->range()->createMember(),
			t4 = DcDy->range()->createMember();
		if(o) *o << "assign(t1.get(),0.5);\n";                    TSFCore::assign(t1.get(),0.5);
		if(o) *o << "assign(t2.get(),0.0);\n";                    TSFCore::assign(t2.get(),0.0);
		if(o) *o << "DcDy->solve(NOTRANS,*t1,t2.get());\n";
		DcDy->solve(NOTRANS,*t1,t2.get());
		if(o && olevel>=VERBOSE_OUTPUT) *o << "*t2 =\n" << *t2;
		if(o) *o << "DcDy->apply(NOTRANS,*t2,t3.get());\n";
		DcDy->apply(NOTRANS,*t2,t3.get());
		if(o && olevel>=VERBOSE_OUTPUT) *o << "*t3 =\n" << *t3;
		if(o) *o << "assign(t4.get(),*t3);\n";                    TSFCore::assign(t4.get(),*t3);
		if(o) *o << "Vp_StV(t4.get(),-1.0,*t1);\n";               TSFCore::Vp_StV(t4.get(),-1.0,*t1);
		if(o) *o << "norm_inf(*t4) -> " << TSFCore::norm_inf(*t4) << std::endl;
		if(o && olevel>=VERBOSE_OUTPUT) *o << "*t4 =\n" << *t4;
	}
	if(prob->adjointSupported()) {
		// Test TRANS
		Teuchos::RefCountPtr<Vector<Scalar> >
			t1 = DcDy->domain()->createMember(),
			t2 = DcDy->range()->createMember(),
			t3 = DcDy->domain()->createMember(),
			t4 = DcDy->domain()->createMember();
		if(o) *o << "assign(t1.get(),0.5);\n";                    TSFCore::assign(t1.get(),0.5);
		if(o) *o << "assign(t2.get(),0.0);\n";                    TSFCore::assign(t2.get(),0.0);
		if(o) *o << "DcDy->solve(TRANS,*t1,t2.get());\n";
		DcDy->solve(TRANS,*t1,t2.get());
		if(o && olevel>=VERBOSE_OUTPUT) *o << "*t2 =\n" << *t2;
		if(o) *o << "DcDy->apply(TRANS,*t2,t3.get());\n";
		DcDy->apply(TRANS,*t2,t3.get());
		if(o && olevel>=VERBOSE_OUTPUT) *o << "*t3 =\n" << *t3;
		if(o) *o << "assign(t4.get(),*t3);\n";                    TSFCore::assign(t4.get(),*t3);
		if(o) *o << "Vp_StV(t4.get(),-1.0,*t1);\n";               TSFCore::Vp_StV(t4.get(),-1.0,*t1);
		if(o) *o << "norm_inf(*t4) -> " << TSFCore::norm_inf(*t4) << std::endl;
		if(o && olevel>=VERBOSE_OUTPUT) *o << "*t4 =\n" << *t4;
	}

	if(o) *o << "\nLeaving NonlinearProblemFirstOrderTester<Scalar>::do_test(prob,...)\n\n";

	return success;
}

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLIN_NONLINEAR_PROBLEM_FIRST_ORDER_TESTER_HPP
