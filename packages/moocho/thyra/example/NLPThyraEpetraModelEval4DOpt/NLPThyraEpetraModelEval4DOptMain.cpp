// //////////////////////////////////////////////////////////////////
// NLPThyraEpetraModelEval4DOptMain.cpp

#include <iostream>

#include "NLPInterfacePack_NLPThyraModelEvaluator.hpp"
#include "EpetraModelEval4DOpt.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#include "MoochoPack_MoochoSolver.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

/* * \mainpage Example optimization problem using <tt>Thyra::ModelEvaluator</tt> and <tt>EpetraExt::ModelEvaluator</tt>.

ToDo: Finish documentation!

*/

int main( int argc, char* argv[] )
{
	using MoochoPack::MoochoSolver;
	using NLPInterfacePack::NLPThyraModelEvaluator;
	using Teuchos::CommandLineProcessor;
	typedef AbstractLinAlgPack::value_type  Scalar;

#ifdef HAVE_MPI
	MPI_Init(&argc,&argv);
#endif

	try {
	
		//
		// Get options from the command line
		//
		
		Scalar       xt0         = 1.0;
		Scalar       xt1         = 1.0;
		Scalar       pt0         = 2.0;
		Scalar       pt1         = 0.0;
		Scalar       d           = 10.0;
    Scalar       x00         = 1.0;
    Scalar       x01         = 1.0;
    Scalar       p00         = 2.0;
    Scalar       p01         = 0.0;
		bool         do_sim      = false;

		CommandLineProcessor  command_line_processor(false); // Don't throw exceptions

		command_line_processor.setOption( "xt0", &xt0 );
		command_line_processor.setOption( "xt1", &xt1 );
		command_line_processor.setOption( "pt0", &pt0 );
		command_line_processor.setOption( "pt1", &pt1 );
		command_line_processor.setOption( "d", &d );
		command_line_processor.setOption( "x00", &x00 );
		command_line_processor.setOption( "x01", &x01 );
		command_line_processor.setOption( "p00", &p00 );
		command_line_processor.setOption( "p01", &p01 );
		command_line_processor.setOption( "do-sim", "do-opt",  &do_sim, "Flag for if only the square constraints are solved" );
	
		CommandLineProcessor::EParseCommandLineReturn
			parse_return = command_line_processor.parse(argc,argv,&std::cerr);

		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
			return parse_return;

		//
		// Create the NLP
		//
		
    // Create the EpetraExt::ModelEvaluator object
    
    EpetraModelEval4DOpt epetraModel(xt0,xt1,pt0,pt1,d,x00,x01,p00,p01);

    // Create the Thyra::EpetraModelEvaluator object
    
    Thyra::EpetraModelEvaluator thyraModel; // Sets default options!
    thyraModel.initialize(
      Teuchos::rcp(&epetraModel,false)
      ,Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory())
      );

    // ToDo: Specify an initial guess not from the Thyra::EpetraModelEvaluator object
    
		NLPThyraModelEvaluator nlp(
			Teuchos::rcp(&thyraModel,false)
			,do_sim ? -1 : 1
			,do_sim ? -1 : 1
			);
    
    //
    // Solve the NLP
    //

		// Create the solver object
		MoochoSolver  solver;

		// Set the NLP
		solver.set_nlp( Teuchos::rcp(&nlp,false) );

		// Solve the NLP
		const MoochoSolver::ESolutionStatus	solution_status = solver.solve_nlp();
		
		//
		// Return the solution status (0 if sucessfull)
		//

		return solution_status;

	}
	catch(const std::exception& excpt) {
		std::cerr << "\nCaught a std::exception " << excpt.what() << std::endl;
	}
	catch(...) {
		std::cerr << "\nCaught an unknown exception\n";
	}

	return MoochoSolver::SOLVE_RETURN_EXCEPTION;
}
