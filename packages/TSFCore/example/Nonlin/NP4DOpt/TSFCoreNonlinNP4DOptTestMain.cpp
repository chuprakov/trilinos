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

// //////////////////////////////////////////////////////////////////////////////////////
// TSFCoreNonlinNP4DOptDeclTestMain.cpp

#include "TSFCoreNonlinNP4DOpt.hpp"
#include "TSFCoreNonlinNonlinearProblemTester.hpp"
#include "TSFCoreNonlinNonlinearProblemFirstOrderTester.hpp"
#include "TSFCoreNonlinSimpleNewtonSolver.hpp"

int main()
{
	namespace mmp = MemMngPack;
	using TSFCore::assign;

	bool result, success = true;

	try {

		std::cout
			<< std::endl
			<< "************************************************\n"
			<< "*** Testing TSFCore::Nonlin::NP4DOpt<> class ***\n"
			<< "************************************************\n";
	
		typedef double Scalar;
		std::cout << "\n*** Construct np4dopt\n";
		TSFCore::Nonlin::NP4DOpt<Scalar> np4dopt;

		std::cout << "\n*** Testing the NonlinearProblem<> interface ...\n\n";
	
		typedef TSFCore::Nonlin::NonlinearProblemTester<Scalar> np_tester_t;
		np_tester_t  np_tester;
	
		result = np_tester.doTest(&np4dopt,&std::cout,np_tester_t::VERBOSE_OUTPUT);
		if(!result) success = 0;

		std::cout << "\n*** Testing the NonlinearProblemFirstOrder<> interface ...\n\n";
	
		typedef TSFCore::Nonlin::NonlinearProblemFirstOrderTester<Scalar> npfo_tester_t;
		npfo_tester_t  npfo_tester;
	
		result = npfo_tester.doTest(&np4dopt,&std::cout,npfo_tester_t::VERBOSE_OUTPUT);
		if(!result) success = 0;

		std::cout << "\n*** Solve the equations c(y)==0 ...\n\n";

		TSFCore::Nonlin::SimpleNewtonSolver<Scalar> newtonSolver;
		newtonSolver.set_summaryOut(Teuchos::rcp(new std::ofstream("NewtonSummary.out")));
		Teuchos::RefCountPtr<TSFCore::Vector<Scalar> > y = np4dopt.space_y()->createMember();
		assign( y.get(), np4dopt.y0() );
		newtonSolver.solve( &np4dopt, y.get(), &std::cout, true );

	}
	catch( const std::exception& excpt) {
		std::cerr << "\nCaught a standard exception of type \'" << typeid(excpt).name() << "\': what() = \"" << excpt.what() << "\"\n";
		success = false;
	}
	catch( ... ) {
		std::cerr << "\nCaught an unknown exception\n";
		success = false;
	}

	if(success)
		std::cout << "\nCongratulations! All of the tests performed passed!\n";
	else
		std::cout << "\nOohs! At least one of the tests failed!\n";
	
	return success ? 0 : 1;

}
