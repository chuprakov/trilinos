// //////////////////////////////////////////////////////////////////////////////////////
// TSFCoreNonlinNP2DSimDeclTestMain.cpp

#include <iostream>
#include <typeinfo>

#include "TSFCoreNonlinNP2DSim.hpp"
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
			<< "*** Testing TSFCore::Nonlin::NP2DSim<> class ***\n"
			<< "************************************************\n";
	
		typedef double Scalar;
		std::cout << "\n*** Construct np2dsim\n";
		TSFCore::Nonlin::NP2DSim<Scalar>  np2dsim;

		std::cout << "\n*** Testing the NonlinearProblem<> interface ...\n\n";
	
		typedef TSFCore::Nonlin::NonlinearProblemTester<Scalar> np_tester_t;
		np_tester_t  np_tester;
	
		result = np_tester.doTest(&np2dsim,&std::cout,np_tester_t::VERBOSE_OUTPUT);
		if(!result) success = 0;

		std::cout << "\n*** Testing the NonlinearProblemFirstOrder<> interface ...\n\n";
	
		typedef TSFCore::Nonlin::NonlinearProblemFirstOrderTester<Scalar> npfo_tester_t;
		npfo_tester_t  npfo_tester;
	
		result = npfo_tester.doTest(&np2dsim,&std::cout,npfo_tester_t::VERBOSE_OUTPUT);
		if(!result) success = 0;

		std::cout << "\n*** Solve the equations ...\n\n";

		TSFCore::Nonlin::SimpleNewtonSolver<Scalar> newtonSolver;
		mmp::ref_count_ptr<TSFCore::Vector<Scalar> > y = np2dsim.space_y()->createMember();
		assign( y.get(), np2dsim.y0() );
		newtonSolver.solve( &np2dsim, y.get(), &std::cout, true );

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
