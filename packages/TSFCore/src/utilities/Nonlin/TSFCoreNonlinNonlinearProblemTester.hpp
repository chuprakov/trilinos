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
// TSFCoreNonlinNonlinearProblemTester.hpp

#ifndef TSFCORE_NONLIN_NONLINEAR_PROBLEM_TESTER_HPP
#define TSFCORE_NONLIN_NONLINEAR_PROBLEM_TESTER_HPP

#include "TSFCoreNonlinNonlinearProblemTesterDecl.hpp"
#include "TSFCoreTestingTools.hpp"
#include "TSFCoreNonlinNonlinearProblem.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVector.hpp"

namespace TSFCore {
namespace Nonlin {

template<class Scalar>
bool
NonlinearProblemTester<Scalar>::doTest(
	NonlinearProblem<Scalar>* prob, std::ostream* o, EOutputLevel olevel
	) const
{
	namespace mmp = MemMngPack;
	using std::endl;
	
	if(o) *o << "Entering NonlinearProblemTester<Scalar>::do_test(prob,...)\n\n";
	//if(o) *o << "Entering NonlinearProblemTester<Scalar>::do_test(prob,...)\n\n" << std::boolalpha;
	// Removed boolalpha because it doesn't exist for GNU 2.96 compilers (HKT 09/22/03)
	// To Do (HKT):  Create boolalpha capabilities 

	bool result, success = true;
	int l;

	// Basic info

	if(o) *o << "prob->initialize(true);\n";         prob->initialize(true);
	if(o) *o << "prob->isInitialized() -> ";         const bool isInitialized = prob->isInitialized();   if(o) *o << isInitialized << endl;
	if(o) *o << "prob->Nu() -> ";                    const int Nu = prob->Nu();                          if(o) *o << Nu << endl;
	if(o) *o << "prob->numResponseFunctions() -> ";  const int nrf = prob->numResponseFunctions();       if(o) *o << nrf << endl;

	// Vector spaces

	typedef Teuchos::RefCountPtr<const VectorSpace<Scalar> >  vs_ptr_t;

	if(o) *o << "prob->space_y().get() -> ";           const vs_ptr_t space_y = prob->space_y();               if(o) *o << space_y.get() << endl;
	if(o) *o << "prob->space_y().get()!=NULL -> ";     if(!(result=(space_y.get()!=NULL))) success = false;    if(o) *o << passfail(result) << endl;
	if(o) *o << "ny = prob->space_y()->dim() -> ";     const Index ny = space_y->dim();                        if(o) *o << ny << endl;

	for( l = 1; l <= Nu; ++l ) {
		if(o) *o << "prob->space_u("<<l<<").get() -> ";                const vs_ptr_t space_u_l = prob->space_u(l);           if(o) *o << space_u_l.get() << endl;
		if(o) *o << "prob->space_u("<<l<<").get()!=NULL -> ";          if(!(result=(space_u_l.get()!=NULL))) success = false; if(o) *o << passfail(result) << endl;
		if(o) *o << "nu("<<l<<") = prob->space_u("<<l<<")->dim() -> "; const Index nu_l = space_u_l->dim();                   if(o) *o << nu_l << endl;
	}

	if(o) *o << "prob->space_c().get() -> ";           const vs_ptr_t space_c = prob->space_c();               if(o) *o << space_c.get() << endl;
	if(o) *o << "prob->space_c().get()!=NULL -> ";     if(!(result=(space_c.get()!=NULL))) success = false;    if(o) *o << passfail(result) << endl;
	if(o) *o << "nc = prob->space_c()->dim() -> ";     const Index nc = space_c->dim();                        if(o) *o << nc << endl;
	if(o) *o << "ny-nc = " << (ny-nc) << " == 0 -> ";  if(!(result=(ny-nc==0))) success = false;               if(o) *o << passfail(result) << endl;

	if(o) *o << "prob->space_g().get() -> ";           const vs_ptr_t space_g = prob->space_g();               if(o) *o << space_g.get() << endl;
	if(nrf) {
		if(o) *o << "prob->space_g().get()!=NULL -> ";                        if(!(result=(space_g.get()!=NULL))) success = false;  if(o) *o << passfail(result) << endl;
		if(o) *o << "ng = prob->space_g()->dim() -> ";                        const Index ng = space_g->dim();                      if(o) *o << nc << endl;
		if(o) *o << "numResponseFunctions-ng -> " << (nrf-ng) << " == 0 -> "; if(!(result=(nrf-ng==0))) success = false;            if(o) *o << passfail(result) << endl;
	}
	else {
		if(o) *o << "prob->space_g().get()==NULL -> ";    if(!(result=(space_g.get()==NULL))) success = false;    if(o) *o << passfail(result) << endl;
	}

	// Bounds and initial guesses

	if(o && olevel>=VERBOSE_OUTPUT) *o << "prob->yL() = \n" << prob->yL();
	if(o && olevel>=VERBOSE_OUTPUT) *o << "prob->yU() = \n" << prob->yU();

	for( l = 1; l <= Nu; ++l ) {
		if(o && olevel>=VERBOSE_OUTPUT) *o << "prob->uL("<<l<<") = \n" << prob->uL(l);
		if(o && olevel>=VERBOSE_OUTPUT) *o << "prob->uU("<<l<<") = \n" << prob->uU(l);
	}

	if(nrf) {
		if(o && olevel>=VERBOSE_OUTPUT) *o << "prob->gL() = \n" << prob->gL();
		if(o && olevel>=VERBOSE_OUTPUT) *o << "prob->gU() = \n" << prob->gU();
	}
	
	if(o && olevel>=VERBOSE_OUTPUT) *o << "prob->y0() = \n" << prob->y0();

	for( l = 1; l <= Nu; ++l ) {
		if(o && olevel>=VERBOSE_OUTPUT) *o << "prob->u0("<<l<<") = \n" << prob->u0(l);
	}

	// Calculation members

	if(o) *o << "c = prob->space_y()->createMember();\n";
	Teuchos::RefCountPtr<Vector<Scalar> > c = prob->space_y()->createMember(), g;
	if(o) *o << "prob->set_c(c)\n";  prob->set_c(c.get());
	if(o) *o << "prob->get_c() -> "; Vector<Scalar> *get_c = prob->get_c(); if(o) *o << get_c << endl;

	if(nrf) {
		if(o) *o << "g = prob->space_y()->createMember();\n";
		g = prob->space_g()->createMember();
		if(o) *o << "prob->set_g(g)\n";  prob->set_g(g.get());
		if(o) *o << "prob->get_g() -> "; Vector<Scalar> *get_g = prob->get_g(); if(o) *o << get_g << endl;
	}

	const Vector<Scalar>                &y0 = prob->y0();
	std::vector<const Vector<Scalar>*>  u0_store(Nu); for(l=1;l<=Nu;++l) u0_store[l-1] = &prob->u0(l);
	const Vector<Scalar>                **u0 = NULL; if(Nu) u0 = &u0_store[0];
	
	if(o) *o << "prob->calc_c(y0,u0,true);\n"; prob->calc_c(y0,u0,true);
	if(o && olevel>=VERBOSE_OUTPUT) *o << "*prob->get_c() = \n" << *prob->get_c();

	if(nrf) {
		if(o) *o << "prob->calc_g(y0,u0,true);\n"; prob->calc_g(y0,u0,true);
		if(o && olevel>=VERBOSE_OUTPUT) *o << "*prob->get_g() = \n" << *prob->get_g();
	}

	// Unset quantities
	
	if(o) *o << "prob->unsetQuantities();\n"; prob->unsetQuantities();
	if(o) *o << "prob->get_c() -> "; get_c = prob->get_c(); if(o) *o << get_c << endl;
	if(nrf) {
		if(o) *o << "prob->get_g() -> "; Vector<Scalar> *get_g = prob->get_g(); if(o) *o << get_g << endl;
	}

	if(o) *o << "\nLeaving NonlinearProblemTester<Scalar>::do_test(prob,...)\n\n";

	return success;
}

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLIN_NONLINEAR_PROBLEM_TESTER_HPP
