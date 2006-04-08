#include "GLpApp_AdvDiffReactOptModel.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "NLPInterfacePack_NLPDirectThyraModelEvaluator.hpp"
#include "NLPInterfacePack_NLPFirstOrderThyraModelEvaluator.hpp"
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#include "Thyra_MultiVectorSerialization.hpp"
#include "Thyra_VectorStdOps.hpp"
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

namespace {

const int maxProcOrder = 4; // Good for 9,999 processors!

std::string getParallelFileName( const std::string &fileNameBase )
{
  std::ostringstream parallelFileName;
  parallelFileName
    << fileNameBase
    << "."
    << std::setfill('0')
    << std::right << std::setw(maxProcOrder)
    << Teuchos::GlobalMPISession::getNProc()
    << "."
    << std::setfill('0')
    << std::right << std::setw(maxProcOrder)
    << Teuchos::GlobalMPISession::getRank()
    ;
  return parallelFileName.str();
}

Teuchos::RefCountPtr<Thyra::VectorBase<double> >
readVectorFromFile(
  const std::string                                                   &fileNameBase
  ,const Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<double> >  &vs
  ,const double                                                       scaleBy = 0.0
  )
{
  const std::string fileName = getParallelFileName(fileNameBase);
  std::ifstream in_file(fileName.c_str());
  TEST_FOR_EXCEPTION(
    in_file.eof(), std::logic_error
    ,"Error, the file \""<<fileName<<"\" could not be opened for input!"
    );
  Teuchos::RefCountPtr<Thyra::VectorBase<double> >
    vec = Thyra::createMember(vs);
  Thyra::MultiVectorSerialization<double> mvSerializer;
  mvSerializer.unserialize(in_file,&*vec);
  Thyra::Vt_S(&*vec,scaleBy);
  return vec;
}

void writeVectorToFile(
  const Thyra::VectorBase<double>    &vec
  ,const std::string                 &fileNameBase
  )
{
  const std::string fileName = getParallelFileName(fileNameBase);
  std::ofstream out_file(fileName.c_str());
  Thyra::MultiVectorSerialization<double> mvSerializer;
  mvSerializer.serialize(vec,out_file);
}

enum ELOWSFactoryType {
  LOWSF_AMESOS
#ifdef HAVE_AZTECOO_THYRA
  ,LOWSF_AZTECOO
#endif
#ifdef HAVE_BELOS_THYRA
  ,LOWSF_BELOS
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

} // namespace

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

  Teuchos::GlobalMPISession mpiSession(&argc,&argv,NULL);

  const int procRank = mpiSession.getRank();
  const int numProcs = mpiSession.getNProc();
  
  bool dummySuccess = true;

  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

	try {
	
		// Create the solver object
		MoochoSolver  solver;

		//
		// Get options from the command line
		//

    const ELOWSFactoryType LOWSFactoryTypeValues[numLOWSFactoryTypes] = {
      LOWSF_AMESOS
#ifdef HAVE_AZTECOO_THYRA
      ,LOWSF_AZTECOO
#endif
#ifdef HAVE_BELOS_THYRA
      ,LOWSF_BELOS
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

    double              len_x           = 1.0;
    double              len_y           = 1.0;
    int                 local_nx        = 3;
    int                 local_ny        = 4;
    std::string         geomFileBase    = "";
    int                 np              = -1;
    double              beta            = 1.0;
    double              x0              = 0.0;
    double              p0              = 1.0;
    double              reactionRate    = 1.0;
    bool                use_direct      = false;
		bool                do_sim          = false;
    ELOWSFactoryType    lowsFactoryType = LOWSF_AMESOS;
#if defined(HAVE_TEUCHOS_EXTENDED) && defined(HAVE_TEUCHOS_EXPAT)
    std::string         lowsfParamsFile = "";
    std::string         lowsfExtraParams = "";
    std::string         lowsfParamsUsedFile = "";
#endif
    bool                usePrec         = true;
    bool                printOnAllProcs = true;
    bool                dump_all        = false;
    std::string         matchingVecFile = "";
    std::string         stateGuessFileBase  = "";
    double              scaleStateGuess     = 1.0;
    std::string         paramGuessFileBase  = "";
    double              scaleParamGuess     = 1.0;
    std::string         stateSoluFileBase   = "";
    std::string         paramSoluFileBase   = "";

		CommandLineProcessor  clp(false); // Don't throw exceptions

		clp.setOption( "len-x", &len_x, "Mesh dimension in the x direction (Overridden by --geom-file-base)." );
		clp.setOption( "len-y", &len_y, "Mesh dimension in the y direction (Overridden by --geom-file-base)." );
		clp.setOption( "local-nx", &local_nx, "Number of local discretization segments in the x direction (Overridden by --geom-file-base)." );
		clp.setOption( "local-ny", &local_ny, "Number of local discretization segments in the y direction (Overridden by --geom-file-base)." );
		clp.setOption( "geom-file-base", &geomFileBase, "Base name of geometry file to read the mesh from." );
		clp.setOption( "np", &np, "The number of optimization parameters (If < 0 then all of boundary is used)" );
		clp.setOption( "beta", &beta, "Regularization." );
		clp.setOption( "x0", &x0, "Initial guess for the state." );
		clp.setOption( "p0", &p0, "Initial guess or nonminal value for optimization parameters." );
		clp.setOption( "reaction-rate", &reactionRate, "The rate of the reaction" );
		clp.setOption( "use-direct", "use-first-order",  &use_direct, "Flag for if we use the NLPDirect or NLPFirstOrderInfo implementation." );
		clp.setOption( "do-sim", "do-opt",  &do_sim, "Flag for if only the square constraints are solved" );
		clp.setOption( "lowsf", &lowsFactoryType
                   ,numLOWSFactoryTypes,LOWSFactoryTypeValues,LOWSFactoryTypeNames
                   ,"The implementation for the LinearOpWithSolveFactory object used to solve the state linear systems" );
#if defined(HAVE_TEUCHOS_EXTENDED) && defined(HAVE_TEUCHOS_EXPAT)
		clp.setOption( "lowsf-params-file", &lowsfParamsFile, "LOWSF parameters XML file (must be compatible with --lowsf=???" );
		clp.setOption( "lowsf-extra-params", &lowsfExtraParams, "Extra LOWSF parameters specified as a string in XML format (must be compatible with --lowsf=???" );
		clp.setOption( "lowsf-params-used-file", &lowsfParamsUsedFile, "File to write the LOWSF parameters that where actually used to." );
#endif
		clp.setOption( "use-prec", "no-use-prec",  &usePrec, "Flag for if preconditioning is used or not" );
    clp.setOption( "print-on-all-procs", "print-on-root-proc", &printOnAllProcs, "Print on all processors or just the root processor?" );
		clp.setOption( "dump-all", "no-dump-all",  &dump_all, "Flag for if we dump everything to STDOUT" );
		clp.setOption( "q-vec-file", &matchingVecFile, "Base file name to read the objective state matching vector q (i.e. ||x-q||_M in objective)." );
		clp.setOption( "x-guess-file", &stateGuessFileBase, "Base file name to read the guess of the state x from." );
		clp.setOption( "scale-x-guess", &scaleStateGuess, "Amount to scale the guess for x read in by --x-guess-file." );
		clp.setOption( "p-guess-file", &paramGuessFileBase, "Base file name to read the guess of the parameters p from." );
		clp.setOption( "scale-p-guess", &scaleParamGuess, "Amount to scale the guess for p read in by --p-guess-file." );
		clp.setOption( "x-solu-file", &stateSoluFileBase, "Base file name to write the state solution x to." );
		clp.setOption( "p-solu-file", &paramSoluFileBase, "Base file name to write the parameter solution p to." );
    solver.setup_commandline_processor(&clp);

		CommandLineProcessor::EParseCommandLineReturn
			parse_return = clp.parse(argc,argv,&std::cerr);

		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
			return parse_return;

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
    journalOut->copyAllOutputOptions(*out);

    *out
      << "\n***"
      << "\n*** NLPThyraEpetraAdvDiffReactOptMain, numProcs="<<numProcs
      << "\n***\n";

    
		//
		// Create the NLP
		//

    *out << "\nCreate the GLpApp::GLpYUEpetraDataPool object ...\n";

#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    GLpApp::GLpYUEpetraDataPool dat(Teuchos::rcp(&comm,false),beta,len_x,len_y,local_nx,local_ny,geomFileBase.c_str(),false);

    *out << "\nCreate the GLpApp::AdvDiffReactOptModel wrapper object ...\n";

    GLpApp::AdvDiffReactOptModel epetraModel(Teuchos::rcp(&dat,false),np,x0,p0,reactionRate);
    epetraModel.setOStream(journalOut);
    if(dump_all) epetraModel.setVerbLevel(Teuchos::VERB_EXTREME);
    
    Teuchos::RefCountPtr<Thyra::LinearOpWithSolveFactoryBase<Scalar> > lowsFactory;
    switch(lowsFactoryType) {
      case LOWSF_AMESOS:
        *out << "\nCreating a Thyra::AmesosLinearOpWithSolveFactory object ...\n";
        lowsFactory = Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory());
        break;
#ifdef HAVE_AZTECOO_THYRA
      case LOWSF_AZTECOO:
        *out << "\nCreating a Thyra::AztecOOLinearOpWithSolveFactory object ...\n";
        lowsFactory = Teuchos::rcp(new Thyra::AztecOOLinearOpWithSolveFactory());
        break;
#endif // HAVE_AZTECOO_THYRA
#ifdef HAVE_BELOS_THYRA
      case LOWSF_BELOS:
        *out << "\nCreating a Thyra::BelosLinearOpWithSolveFactory object ...\n";
        lowsFactory = Teuchos::rcp(new Thyra::BelosLinearOpWithSolveFactory<Scalar>());
        break;
#endif // HAVE_BELOS_THYRA
      default:
        TEST_FOR_EXCEPT(true); // should never get here!
    }
#ifdef HAVE_IFPACK_THYRA
    if( usePrec && lowsFactory->acceptsPreconditionerFactory() ) {
      *out << "\nCreating a Thyra::IfpackPreconditionerFactory object ...\n";
      lowsFactory->setPreconditionerFactory(Teuchos::rcp(new Thyra::IfpackPreconditionerFactory()),"");
    }
#endif // HAVE_IFPACK_THYRA
    Teuchos::RefCountPtr<Teuchos::ParameterList>
      lowsfPL = Teuchos::rcp(new Teuchos::ParameterList("LOWSF"));
    if(1) {
#if defined(HAVE_TEUCHOS_EXTENDED) && defined(HAVE_TEUCHOS_EXPAT)
      Teuchos::XMLParameterListReader xmlPLReader;
      if(lowsfParamsFile.length()) {
        Teuchos::FileInputSource xmlFile(lowsfParamsFile);
        Teuchos::XMLObject xmlParams = xmlFile.getObject();
        lowsfPL->setParameters(xmlPLReader.toParameterList(xmlParams));
        *journalOut << "\nLOWSF parameters read from the file \""<<lowsfParamsFile<<"\":\n";
        lowsfPL->print(*OSTab(journalOut).getOStream(),0,true);
      }
      if(lowsfExtraParams.length()) {
        Teuchos::StringInputSource xmlStr(lowsfExtraParams);
        Teuchos::XMLObject xmlParams = xmlStr.getObject();
        lowsfPL->setParameters(xmlPLReader.toParameterList(xmlParams));
        *journalOut << "\nExtra LOWSF parameters taken from the command-line:\n";
        lowsfPL->print(*OSTab(journalOut).getOStream(),0,true);
      }
#endif // defined(HAVE_TEUCHOS_EXTENDED) && defined(HAVE_TEUCHOS_EXPAT)
      lowsFactory->setParameterList(lowsfPL);
      *journalOut << "\nList of all valid LOWSF parameters:\n";
      OSTab tab(journalOut);
      *journalOut << lowsFactory->getValidParameters()->name() << " ->\n";
      tab.incrTab();
      lowsFactory->getValidParameters()->print(*journalOut,0,true);
    }
    
    *out << "\nCreate the Thyra::EpetraModelEvaluator wrapper object ...\n";
    
    Thyra::EpetraModelEvaluator thyraModel; // Sets default options!
    thyraModel.setOStream(journalOut);
    thyraModel.initialize(Teuchos::rcp(&epetraModel,false),lowsFactory);

    if(matchingVecFile != "") {
      *out << "\nReading the matching vector \'q\' from the file(s) with base name \""<<matchingVecFile<<"\" ...\n";
      epetraModel.set_q(
        Thyra::get_Epetra_Vector(
          *epetraModel.get_x_map(),readVectorFromFile(matchingVecFile,thyraModel.get_x_space()
            )
          )
        );
    }
    
    if(1) {
      Thyra::ModelEvaluatorBase::InArgs<double> thyraModel_initialGuess = thyraModel.createInArgs();
      if(stateGuessFileBase != "") {
        *out << "\nReading the guess of the state \'x\' from the file(s) with base name \""<<stateGuessFileBase<<"\" ...\n";
        thyraModel_initialGuess.set_x(readVectorFromFile(stateGuessFileBase,thyraModel.get_x_space(),scaleStateGuess));
      }
      if(paramGuessFileBase != "") {
        *out << "\nReading the guess of the parameters \'p\' from the file(s) with base name \""<<paramGuessFileBase<<"\" ...\n";
        thyraModel_initialGuess.set_p(0,readVectorFromFile(paramGuessFileBase,thyraModel.get_p_space(0),scaleParamGuess));
      }
      thyraModel.setInitialGuess(thyraModel_initialGuess);
    }
    
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
    
		// Set the NLP
		solver.set_nlp(nlp);

		// Solve the NLP
		const MoochoSolver::ESolutionStatus	solution_status = solver.solve_nlp();

    if(stateSoluFileBase != "") {
      *out << "\nWriting the state solution \'x\' to the file(s) with base name \""<<stateSoluFileBase<<"\" ...\n";
      writeVectorToFile(*thyraModel.getFinalPoint().get_x(),stateSoluFileBase);
    }
    if( paramSoluFileBase != "" ) {
      *out << "\nWriting the parameter solution \'p\' to the file(s) with base name \""<<paramSoluFileBase<<"\" ...\n";
      writeVectorToFile(*thyraModel.getFinalPoint().get_p(0),paramSoluFileBase);
    }

    // Write the LOWSF parameters that were used:
#if defined(HAVE_TEUCHOS_EXTENDED) && defined(HAVE_TEUCHOS_EXPAT)
    if(lowsfParamsUsedFile != "" && procRank == 0) {
      Teuchos::XMLParameterListWriter plWriter;
      Teuchos::XMLObject xml = plWriter.toXML(*lowsfPL);
      std::ofstream of(lowsfParamsUsedFile.c_str());
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
