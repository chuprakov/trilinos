#include "GLpApp_AdvDiffReactOptModel.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "NLPInterfacePack_NLPThyraModelEvaluator.hpp"
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#include "MoochoPack_MoochoSolver.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_oblackholestream.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

int main( int argc, char* argv[] )
{
	using MoochoPack::MoochoSolver;
	using NLPInterfacePack::NLPThyraModelEvaluator;
	using Teuchos::CommandLineProcessor;
	typedef AbstractLinAlgPack::value_type  Scalar;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  const int procRank = Teuchos::GlobalMPISession::getRank();
  //const int numProcs = Teuchos::GlobalMPISession::getNProc();

  bool dummySuccess = true;

  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

	try {
	
		//
		// Get options from the command line
		//

    std::string  geomFileBase    = "";
    double       beta            = 1.0;
    double       x0              = 0.0;
    double       p0              = 1.0;
    double       reactionRate    = 1.0;
		bool         do_sim          = false;
    bool         printOnAllProcs = true;
    bool         dump_all        = false;

		CommandLineProcessor  command_line_processor(false); // Don't throw exceptions

		command_line_processor.setOption( "geom-file-base", &geomFileBase, "Base name of geometry file." );
		command_line_processor.setOption( "beta", &beta, "Regularization." );
		command_line_processor.setOption( "x0", &x0, "Initial guess for the state." );
		command_line_processor.setOption( "p0", &p0, "Initial guess or nonminal value for control." );
		command_line_processor.setOption( "reaction-rate", &reactionRate, "The rate of the reaction" );
		command_line_processor.setOption( "do-sim", "do-opt",  &do_sim, "Flag for if only the square constraints are solved" );
    command_line_processor.setOption( "print-on-all-procs", "print-on-root-proc", &printOnAllProcs, "Print on all processors or just the root processor?" );
		command_line_processor.setOption( "dump-all", "no-dump-all",  &dump_all, "Flag for if we dump everything to STDOUT" );

		CommandLineProcessor::EParseCommandLineReturn
			parse_return = command_line_processor.parse(argc,argv,&std::cerr);

		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
			return parse_return;

    TEST_FOR_EXCEPTION(geomFileBase=="",std::logic_error,"Error, you must specify a geometry file as --geom-file-base=???, see --help");

    //
    // Setup the output streams
    //

    Teuchos::oblackholestream black_hole_out;
    std::ostream &this_proc_out = ( procRank==0 || printOnAllProcs ? std::cout : black_hole_out );
    Teuchos::VerboseObjectBase::setDefaultOStream(
      Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcp(&this_proc_out,false),"  ")));
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

		//
		// Create the NLP
		//

    *out << "\nCreate the GLpApp::GLpYUEpetraDataPool object ...\n";

#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    GLpApp::GLpYUEpetraDataPool dat ( Teuchos::rcp(&comm,false), beta, geomFileBase.c_str(), false );

    *out << "\nCreate the GLpApp::AdvDiffReactOptModel wrapper object ...\n";

    GLpApp::AdvDiffReactOptModel epetraModel(Teuchos::rcp(&dat,false),x0,p0,reactionRate,dump_all);

    *out << "\nCreate the Thyra::EpetraModelEvaluator wrapper object ...\n";
    
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
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,dummySuccess)

	return MoochoSolver::SOLVE_RETURN_EXCEPTION;
}
