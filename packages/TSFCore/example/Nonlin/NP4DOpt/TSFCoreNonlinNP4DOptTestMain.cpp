// //////////////////////////////////////////////////////////////////////////////////////
// TSFCoreNonlinNP4DOptDeclTestMain.cpp

#include <iostream>
#include <typeinfo>

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
		mmp::ref_count_ptr<TSFCore::Vector<Scalar> > y = np4dopt.space_y()->createMember();
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
		std::cout << "\nCongradulations! All of the tests performed passed!\n";
	else
		std::cout << "\nOohs! At least one of the tests failed!\n";
	
	return success ? 0 : -1;

}
