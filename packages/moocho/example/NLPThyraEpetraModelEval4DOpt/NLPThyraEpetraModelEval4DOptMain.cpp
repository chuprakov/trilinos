// //////////////////////////////////////////////////////////////////
// NLPThyraEpetraModelEval4DOptMain.cpp

#include <iostream>

#include "NLPInterfacePack/src/abstract/thyra/NLPThyraModelEvaluator.hpp"
#include "EpetraModelEval4DOpt.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#include "MoochoPack/configurations/MoochoSolver.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

/** \mainpage Example optimization problem using <tt>Thyra::ModelEvaluator</tt> and <tt>EpetraExt::ModelEvaluator</tt>.

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
		
		Scalar       yt1         = 1.0;
		Scalar       yt2         = 1.0;
		Scalar       ut1         = 2.0;
		Scalar       ut2         = 0.0;
		Scalar       d           = 10.0;
    Scalar       y01         = 1.0;
    Scalar       y02         = 1.0;
    Scalar       u01         = 2.0;
    Scalar       u02         = 0.0;
		bool         do_sim      = false;

		CommandLineProcessor  command_line_processor(false); // Don't throw exceptions

		command_line_processor.setOption( "yt1", &yt1 );
		command_line_processor.setOption( "yt2", &yt2 );
		command_line_processor.setOption( "ut1", &ut1 );
		command_line_processor.setOption( "ut2", &ut2 );
		command_line_processor.setOption( "d", &d );
		command_line_processor.setOption( "y01", &y01 );
		command_line_processor.setOption( "y02", &y02 );
		command_line_processor.setOption( "u01", &u01 );
		command_line_processor.setOption( "u02", &u02 );
		command_line_processor.setOption( "do-sim", "do-opt",  &do_sim, "Flag for if only the square constraints are solved" );
	
		CommandLineProcessor::EParseCommandLineReturn
			parse_return = command_line_processor.parse(argc,argv,&std::cerr);

		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
			return parse_return;

		//
		// Create the NLP
		//
		
    // Create the Epetra::NonlinearProblemFirstOrder object
    
    EpetraModelEval4DOpt epetra_np(yt1,yt2,ut1,ut2,d,y01,y02,u01,u02);

    // Create the TSFCore::Nonlin::NonlinearProblemFirstOrder object
    
    Thyra::EpetraModelEvaluator np; // Sets default options!
    np.initialize(
      Teuchos::rcp(&epetra_np,false)
      ,Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory())
      );
    
    // Setup NLPTSFCoreNP for just one objective function

		const int                              u_indep_ind[1]     = { 1 };
		const int                              num_obj            = 1;
		const int                              obj_ind[num_obj]   = { 1 };
		const Scalar                           obj_wgt[num_obj]   = { 1.0 };
		const NLPThyraModelEvaluator::EObjPow  obj_pow[num_obj]   = { NLPThyraModelEvaluator::OBJ_LINEAR };

    // ToDo: Specify an initial guess not from the TSFCore::Nonlin::NonlinearProblem object

		NLPThyraModelEvaluator nlp(
			Teuchos::rcp(&np,false)
			,do_sim ? 0    : 1
			,do_sim ? 0    : u_indep_ind 
			,do_sim ? 0    : num_obj
			,do_sim ? NULL : obj_ind
			,do_sim ? NULL : obj_wgt
			,do_sim ? NULL : obj_pow
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
