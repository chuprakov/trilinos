#include "GLpApp_AdvDiffReactOptModel.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "NLPInterfacePack_NLPDirectThyraModelEvaluator.hpp"
#include "NLPInterfacePack_NLPFirstOrderThyraModelEvaluator.hpp"
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#include "MoochoPack_MoochoSolver.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_oblackholestream.hpp"
#ifdef HAVE_AZTECOO_THYRA
#  include "Thyra_AztecOOLinearOpWithSolveFactory.hpp"
#endif
#ifdef HAVE_BELOS_THYRA
#  include "Thyra_BelosLinearOpWithSolveFactory.hpp"
#endif
#ifdef HAVE_IFPACK_THYRA
#  include "Thyra_IfpackPreconditionerFactory.hpp"
#endif
#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif
#if defined(HAVE_TEUCHOS_EXTENDED) && defined(HAVE_TEUCHOS_EXPAT)
#  include "Teuchos_FileInputSource.hpp"
#  include "Teuchos_StringInputSource.hpp"
#  include "Teuchos_XMLParameterListReader.hpp"
#  include "Teuchos_XMLParameterListWriter.hpp"
#endif

enum ELOWSFactoryType {
  AMESOS_LOWSF
#ifdef HAVE_AZTECOO_THYRA
  ,AZTECOO_LOWSF
#endif
#ifdef HAVE_BELOS_THYRA
  ,BELOS_LOWSF
#endif
};

const int numLOWSFactoryTypes =
1
#ifdef HAVE_AZTECOO_THYRA
+1
#endif
#ifdef HAVE_BELOS_THYRA
+1
#endif
;

int main( int argc, char* argv[] )
{
	using Teuchos::rcp;
	using Teuchos::OSTab;
	using MoochoPack::MoochoSolver;
	using NLPInterfacePack::NLP;
	using NLPInterfacePack::NLPDirectThyraModelEvaluator;
	using NLPInterfacePack::NLPFirstOrderThyraModelEvaluator;
	using Teuchos::CommandLineProcessor;
	typedef AbstractLinAlgPack::value_type  Scalar;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  bool dummySuccess = true;

  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

	try {
	
		// Create the solver object
		MoochoSolver  solver;

		//
		// Get options from the command line
		//

    const ELOWSFactoryType LOWSFactoryTypeValues[numLOWSFactoryTypes] = {
      AMESOS_LOWSF
#ifdef HAVE_AZTECOO_THYRA
      ,AZTECOO_LOWSF
#endif
#ifdef HAVE_BELOS_THYRA
      ,BELOS_LOWSF
#endif
    };
    const char* LOWSFactoryTypeNames[numLOWSFactoryTypes] = {
      "amesos"
#ifdef HAVE_AZTECOO_THYRA
      ,"aztecoo"
#endif
#ifdef HAVE_BELOS_THYRA
      ,"belos"
#endif
    };

    std::string         geomFileBase    = "";
    double              beta            = 1.0;
    double              x0              = 0.0;
    double              p0              = 1.0;
    double              reactionRate    = 1.0;
    bool                use_direct      = false;
		bool                do_sim          = false;
    ELOWSFactoryType    lowsFactoryType = AMESOS_LOWSF;
#if defined(HAVE_TEUCHOS_EXTENDED) && defined(HAVE_TEUCHOS_EXPAT)
    std::string         lowsfParamsFile = "";
    std::string         lowsfExtraParams = "";
#endif
    bool                printOnAllProcs = true;
    bool                dump_all        = false;

		CommandLineProcessor  clp(false); // Don't throw exceptions

		clp.setOption( "geom-file-base", &geomFileBase, "Base name of geometry file." );
		clp.setOption( "beta", &beta, "Regularization." );
		clp.setOption( "x0", &x0, "Initial guess for the state." );
		clp.setOption( "p0", &p0, "Initial guess or nonminal value for control." );
		clp.setOption( "reaction-rate", &reactionRate, "The rate of the reaction" );
		clp.setOption( "use-direct", "use-first-order",  &use_direct, "Flag for if we use the NLPDirect or NLPFirstOrderInfo implementation." );
		clp.setOption( "do-sim", "do-opt",  &do_sim, "Flag for if only the square constraints are solved" );
		clp.setOption( "lowsf", &lowsFactoryType
                   ,numLOWSFactoryTypes,LOWSFactoryTypeValues,LOWSFactoryTypeNames
                   ,"The implementation for the LinearOpWithSolveFactory object used to solve the state linear systems" );
#if defined(HAVE_TEUCHOS_EXTENDED) && defined(HAVE_TEUCHOS_EXPAT)
		clp.setOption( "lowsf-params-file", &lowsfParamsFile, "LOWSF parameters XML file (must be compatible with --lowsf=???" );
		clp.setOption( "lowsf-extra-params", &lowsfExtraParams, "Extra LOWSF parameters specified as a string in XML format (must be compatible with --lowsf=???" );
#endif
    clp.setOption( "print-on-all-procs", "print-on-root-proc", &printOnAllProcs, "Print on all processors or just the root processor?" );
		clp.setOption( "dump-all", "no-dump-all",  &dump_all, "Flag for if we dump everything to STDOUT" );
    solver.setup_commandline_processor(&clp);

		CommandLineProcessor::EParseCommandLineReturn
			parse_return = clp.parse(argc,argv,&std::cerr);

		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
			return parse_return;

    TEST_FOR_EXCEPTION(geomFileBase=="",std::logic_error,"Error, you must specify a geometry file as --geom-file-base=???, see --help");

    //
    // Setup the output streams
    //

    Teuchos::RefCountPtr<Teuchos::FancyOStream>
      journalOut = Teuchos::rcp(
        new Teuchos::FancyOStream(
          Teuchos::rcp(new std::ofstream("MoochoJournal.out"))
          ,"  "
          )
        );

		//
		// Create the NLP
		//

    *out << "\nCreate the GLpApp::GLpYUEpetraDataPool object ...\n";

#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    GLpApp::GLpYUEpetraDataPool dat(Teuchos::rcp(&comm,false),beta,geomFileBase.c_str(),false);

    *out << "\nCreate the GLpApp::AdvDiffReactOptModel wrapper object ...\n";

    GLpApp::AdvDiffReactOptModel epetraModel(Teuchos::rcp(&dat,false),x0,p0,reactionRate);
    epetraModel.setOStream(journalOut);
    if(dump_all) epetraModel.setVerbLevel(Teuchos::VERB_EXTREME);
    
    Teuchos::RefCountPtr<Thyra::LinearOpWithSolveFactoryBase<Scalar> > lowsFactory;
    switch(lowsFactoryType) {
      case AMESOS_LOWSF:
        *out << "\nCreating a Thyra::AmesosLinearOpWithSolveFactory object ...\n";
        lowsFactory = Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory());
        break;
#ifdef HAVE_AZTECOO_THYRA
      case AZTECOO_LOWSF:
        *out << "\nCreating a Thyra::AztecOOLinearOpWithSolveFactory object ...\n";
        lowsFactory = Teuchos::rcp(new Thyra::AztecOOLinearOpWithSolveFactory());
        break;
#endif // HAVE_AZTECOO_THYRA
#ifdef HAVE_BELOS_THYRA
      case BELOS_LOWSF:
        *out << "\nCreating a Thyra::BelosLinearOpWithSolveFactory object ...\n";
        lowsFactory = Teuchos::rcp(new Thyra::BelosLinearOpWithSolveFactory<Scalar>());
#ifdef HAVE_IFPACK_THYRA
        *out << "\nCreating a Thyra::IfpackPreconditionerFactory object ...\n";
        lowsFactory->setPreconditionerFactory(Teuchos::rcp(new Thyra::IfpackPreconditionerFactory()),"");
#endif // HAVE_IFPACK_THYRA
        break;
#endif // HAVE_BELOS_THYRA
      default:
        TEST_FOR_EXCEPT(true); // should never get here!
    }
    Teuchos::RefCountPtr<Teuchos::ParameterList>
      lowsfPL = Teuchos::rcp(new Teuchos::ParameterList("LOWSF"));
    if(1) {
#if defined(HAVE_TEUCHOS_EXTENDED) && defined(HAVE_TEUCHOS_EXPAT)
      Teuchos::XMLParameterListReader xmlPLReader;
      if(lowsfParamsFile.length()) {
        Teuchos::FileInputSource xmlFile(lowsfParamsFile);
        Teuchos::XMLObject xmlParams = xmlFile.getObject();
        lowsfPL->setParameters(xmlPLReader.toParameterList(xmlParams));
        *out << "\nLOWSF parameters read from the file \""<<lowsfParamsFile<<"\":\n";
        lowsfPL->print(*OSTab(out).getOStream(),0,true);
      }
      if(lowsfExtraParams.length()) {
        Teuchos::StringInputSource xmlStr(lowsfExtraParams);
        Teuchos::XMLObject xmlParams = xmlStr.getObject();
        lowsfPL->setParameters(xmlPLReader.toParameterList(xmlParams));
        *out << "\nExtra LOWSF parameters taken from the command-line:\n";
        lowsfPL->print(*OSTab(out).getOStream(),0,true);
      }
#endif // defined(HAVE_TEUCHOS_EXTENDED) && defined(HAVE_TEUCHOS_EXPAT)
      lowsFactory->setParameterList(lowsfPL);
      *out << "\nList of all valid LOWSF parameters:\n";
      OSTab tab(out);
      *out << lowsFactory->getValidParameters()->name() << " ->\n";
      tab.incrTab();
      lowsFactory->getValidParameters()->print(*out,0,true);
    }
    
    *out << "\nCreate the Thyra::EpetraModelEvaluator wrapper object ...\n";
    
    Thyra::EpetraModelEvaluator thyraModel; // Sets default options!
    thyraModel.setOStream(journalOut);
    thyraModel.initialize(Teuchos::rcp(&epetraModel,false),lowsFactory);
    
    Teuchos::RefCountPtr<NLP> nlp;
    if(use_direct) {
      nlp = rcp(
        new NLPDirectThyraModelEvaluator(
          Teuchos::rcp(&thyraModel,false)
          ,do_sim ? -1 : 0
          ,do_sim ? -1 : 0
          )
        );
    }
    else {
      nlp = rcp(
        new NLPFirstOrderThyraModelEvaluator(
          Teuchos::rcp(&thyraModel,false)
          ,do_sim ? -1 : 0
          ,do_sim ? -1 : 0
          )
        );
    }
    
    //
    // Solve the NLP
    //

    // Set the journal file
    solver.set_journal_out(journalOut);

/*    
    // Create the initial options group that will be overridden
    if(1) {
      Teuchos::RefCountPtr<OptionsFromStreamPack::OptionsFromStream>
        moochoOptions = Teuchos::rcp(new OptionsFromStreamPack::OptionsFromStream());
      std::istringstream iss("DecompositionSystemStateStepBuilderStd{range_space_matrix=ORTHOGONAL}");
      moochoOptions->read_options(iss);
      solver.set_options(moochoOptions);
    }
*/
    
		// Set the NLP
		solver.set_nlp(nlp);

		// Solve the NLP
		const MoochoSolver::ESolutionStatus	solution_status = solver.solve_nlp();

    // Write the LOWSF parameters that were used:
#if defined(HAVE_TEUCHOS_EXTENDED) && defined(HAVE_TEUCHOS_EXPAT)
    if(1) {
      Teuchos::XMLParameterListWriter plWriter;
      Teuchos::XMLObject xml = plWriter.toXML(*lowsfPL);
      std::ofstream of("lowsfParams.used.xml");
      of << xml << std::endl;
    }
#endif // defined(HAVE_TEUCHOS_EXTENDED) && defined(HAVE_TEUCHOS_EXPAT)
		
		//
		// Return the solution status (0 if sucessfull)
		//

		return solution_status;

	}
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,dummySuccess)

	return MoochoSolver::SOLVE_RETURN_EXCEPTION;
}
