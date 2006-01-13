// //////////////////////////////////////////////////////////////////
// NLPThyraEpetraModelEval4DOptMain.cpp

#include "GLpApp_AdvDiffReactOptModel.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "NLPInterfacePack_NLPThyraModelEvaluator.hpp"
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#include "MoochoPack_MoochoSolver.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

/* * \mainpage Example optimization problem using <tt>Thyra::ModelEvaluator</tt> and <tt>EpetraExt::ModelEvaluator</tt>.

ToDo: Finish documentation!

*/

int main( int argc, char* argv[] )
{
	using MoochoPack::MoochoSolver;
	using NLPInterfacePack::NLPThyraModelEvaluator;
	using Teuchos::CommandLineProcessor;
	typedef AbstractLinAlgPack::value_type  Scalar;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

	try {
	
		//
		// Get options from the command line
		//

    std::string  geomFile    = "";
    double       beta        = 1.0;
    double       x0          = 0.0;
    double       p0          = 1.0;
		bool         do_sim      = false;

		CommandLineProcessor  command_line_processor(false); // Don't throw exceptions

		command_line_processor.setOption( "geom-file", &geomFile, "Base name of geometry file." );
		command_line_processor.setOption( "beta", &beta, "Regularization." );
		command_line_processor.setOption( "x0", &x0, "Initial guess for the state." );
		command_line_processor.setOption( "p0", &p0, "Initial guess or nonminal value for control." );
		command_line_processor.setOption( "do-sim", "do-opt",  &do_sim, "Flag for if only the square constraints are solved" );
	
		CommandLineProcessor::EParseCommandLineReturn
			parse_return = command_line_processor.parse(argc,argv,&std::cerr);

		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
			return parse_return;

    TEST_FOR_EXCEPTION(geomFile=="",std::logic_error,"Error, you must specify a geometry file as --geom-file=???, see --help");

		//
		// Create the NLP
		//

    // Create the GenSQP data pool object

#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    GLpApp::GLpYUEpetraDataPool dat ( Teuchos::rcp(&comm,false), beta, geomFile.c_str() );
		
    // Create the Epetra-centric model

    GLpApp::AdvDiffReactOptModel epetraModel(Teuchos::rcp(&dat,false),x0,p0);

    // Create the Thyra-wrapped model with a linear solver set
    
    Thyra::EpetraModelEvaluator thyraModel; // Sets default options!
    thyraModel.initialize(
      Teuchos::rcp(&epetraModel,false)
      ,Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory())
      );
    
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
