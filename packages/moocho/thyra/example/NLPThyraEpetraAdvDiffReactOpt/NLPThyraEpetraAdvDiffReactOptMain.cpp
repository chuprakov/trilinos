#include "GLpApp_AdvDiffReactOptModel.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "NLPInterfacePack_NLPDirectThyraModelEvaluator.hpp"
#include "NLPInterfacePack_NLPFirstOrderThyraModelEvaluator.hpp"
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultFiniteDifferenceModelEvaluator.hpp"
#include "Thyra_DefaultStateEliminationModelEvaluator.hpp"
#include "Thyra_DefaultEvaluationLoggerModelEvaluator.hpp"
#include "Thyra_ParallelMultiVectorFileIO.hpp"
#include "Thyra_DampenedNewtonNonlinearSolver.hpp"
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
//#  include "Teuchos_FileInputSource.hpp"
//#  include "Teuchos_StringInputSource.hpp"
//#  include "Teuchos_XMLParameterListReader.hpp"
//#  include "Teuchos_XMLParameterListWriter.hpp"
#  include "Teuchos_XMLParameterListHelpers.hpp"
#endif

namespace {

typedef AbstractLinAlgPack::value_type  Scalar;

Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >
readVectorFromFile(
  const std::string                                                   &fileNameBase
  ,const Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> >  &vs
  ,const double                                                       scaleBy = 1.0
  )
{
  Thyra::ParallelMultiVectorFileIO<Scalar> fileIO;
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >
    vec = fileIO.readVectorFromFile(fileNameBase,vs);
  Thyra::Vt_S(&*vec,scaleBy);
  return vec;
}

void writeVectorToFile(
  const Thyra::VectorBase<Scalar>    &vec
  ,const std::string                 &fileNameBase
  )
{
  Thyra::ParallelMultiVectorFileIO<Scalar> fileIO;
  fileIO.writeToFile(vec,fileNameBase);
}

//
// Linear Solver option
//

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

//
// Finite difference method
//

const int numFDMethodTypes = 7;

Thyra::DirectionalFiniteDiffCalculatorTypes::EFDMethodType fdMethodTypeValues[numFDMethodTypes] =
{
  Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_ONE
  ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_TWO
  ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_TWO_CENTRAL
  ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_TWO_AUTO
  ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_FOUR
  ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_FOUR_CENTRAL
  ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_FOUR_AUTO
};

const char* fdMethodTypeNames[numFDMethodTypes] =
{
  "order-one"
  ,"order-two"
  ,"order-two-central"
  ,"order-two-auto"
  ,"order-four"
  ,"order-four-central"
  ,"order-four-auto"
};

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
  namespace DFDT = Thyra::DirectionalFiniteDiffCalculatorTypes;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  const int procRank = mpiSession.getRank();
  const int numProcs = mpiSession.getNProc();

  Teuchos::Time timer("");
  
  bool dummySuccess = true;

  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {
  
    // Create the solver object
    MoochoSolver  solver;

    //
    // Get options from the command line
    //

    double              len_x           = 1.0;
    double              len_y           = 1.0;
    int                 local_nx        = 3;
    int                 local_ny        = 4;
    std::string         geomFileBase    = "";
    int                 np              = -1;
    bool                normalizeBasis  = false;
    double              beta            = 1.0;
    double              x0              = 0.0;
    double              p0              = 1.0;
    double              reactionRate    = 1.0;
    bool                do_sim          = false;
    bool                use_direct      = false;
    bool                use_black_box   = false;
    bool                use_finite_diff = false;
    DFDT::EFDMethodType fd_method_type  = DFDT::FD_ORDER_ONE;
    double              fd_step_len     = 1e-4;
    double              fwd_newton_tol  = -1.0; // Determine automatically
    int                 fwd_newton_max_iters  = 20;
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

    int                 numProcsPerCluster  = -1;

    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    clp.setOption( "len-x", &len_x, "Mesh dimension in the x direction (Overridden by --geom-file-base)." );
    clp.setOption( "len-y", &len_y, "Mesh dimension in the y direction (Overridden by --geom-file-base)." );
    clp.setOption( "local-nx", &local_nx, "Number of local discretization segments in the x direction (Overridden by --geom-file-base)." );
    clp.setOption( "local-ny", &local_ny, "Number of local discretization segments in the y direction (Overridden by --geom-file-base)." );
    clp.setOption( "geom-file-base", &geomFileBase, "Base name of geometry file to read the mesh from." );
    clp.setOption( "np", &np, "The number of optimization parameters p (If < 0 then all of boundary is used)" );
    clp.setOption( "normalize-basis", "no-normalize-basis", &normalizeBasis, "Normalize the basis for the parameters p or not." );
    clp.setOption( "beta", &beta, "Regularization." );
    clp.setOption( "x0", &x0, "Initial guess for the state." );
    clp.setOption( "p0", &p0, "Initial guess or nonminal value for optimization parameters." );
    clp.setOption( "reaction-rate", &reactionRate, "The rate of the reaction" );
    clp.setOption( "do-sim", "do-opt",  &do_sim, "Flag for if only the square constraints are solved" );
    clp.setOption( "use-direct", "use-first-order",  &use_direct, "Flag for if we use the NLPDirect or NLPFirstOrderInfo implementation." );
    clp.setOption( "black-box", "all-at-once",  &use_black_box, "Flag for if we do black-box NAND or all-at-once SAND." );
    clp.setOption( "use-finite-diff", "no-use-finite-diff",  &use_finite_diff, "Flag for if we are using finite differences or not." );
    clp.setOption( "finite-diff-type", &fd_method_type
                   ,numFDMethodTypes,fdMethodTypeValues,fdMethodTypeNames
                   ,"The type of finite differences used [---use-finite-diff]" );
    clp.setOption( "fd-step-len", &fd_step_len, "The finite-difference step length ( < 0 determine automatically )." );
    clp.setOption( "fwd-newton-tol", &fwd_newton_tol, "Nonlinear residual tolerance for the forward Newton solve." );
    clp.setOption( "fwd-newton-max-iters", &fwd_newton_max_iters, "Max number of iterations for the forward Newton solve." );
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
    clp.setOption( "num-procs-per-cluster", &numProcsPerCluster, "Number of processors in a cluster (<=0 means only one cluster)." );
    solver.setup_commandline_processor(&clp);

    CommandLineProcessor::EParseCommandLineReturn
      parse_return = clp.parse(argc,argv,&std::cerr);

    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
      return parse_return;

    //
    // Setup the output streams
    //
    
    std::string journalOutName;
    if( numProcs > 1 && out->getOutputToRootOnly() < 0 ) {
      std::ostringstream oss;
      oss << "MoochoJournal."<<std::setfill('0')<<std::setw(4)<<std::right<<procRank<<".out";
      journalOutName = oss.str();
    }
    else {
      journalOutName = "MoochoJournal.out";
    }
    
    Teuchos::RefCountPtr<Teuchos::FancyOStream>
      journalOut = Teuchos::rcp(
        new Teuchos::FancyOStream(
          Teuchos::rcp(new std::ofstream(journalOutName.c_str()))
          ,"  "
          )
        );
    journalOut->copyAllOutputOptions(*out);

    *out
      << "\n***"
      << "\n*** NLPThyraEpetraAdvDiffReactOptMain, Global numProcs = "<<numProcs
      << "\n***\n";

    int clusterRank = -1;
    int numClusters = -1;
#ifdef HAVE_MPI
    MPI_Comm mpiComm = MPI_COMM_WORLD;
    if( numProcsPerCluster > 0 ) {
      *out << "\nCreating communicator for local cluster of "<<numProcsPerCluster<<" processors ...\n";
      numClusters = numProcs/numProcsPerCluster;
      const int remainingProcs = numProcs%numProcsPerCluster;
      TEST_FOR_EXCEPTION(
        remainingProcs!=0,std::logic_error
        ,"Error, The number of processors per cluster numProcsPerCluster="<<numProcsPerCluster
        << " does not divide into the global number of processors numProcs="<<numProcs
        << " and instead has remainder="<<remainingProcs<<"!"
        );
      // Determine which cluster this processor is part of and what the global
      // processor ranges are.
      clusterRank = procRank / numProcsPerCluster; // Integer division!
      *out << "\nclusterRank = " << clusterRank << "\n";
      const int firstClusterProcRank = clusterRank * numProcsPerCluster;
      const int lastClusterProcRank = firstClusterProcRank + numProcsPerCluster - 1;
      *out << "\nclusterProcRange = ["<<firstClusterProcRank<<","<<lastClusterProcRank<<"]\n";
      // Create the group and the communicator for this cluster of processors
      MPI_Group globalGroup = MPI_GROUP_NULL;
      MPI_Comm_group(MPI_COMM_WORLD,&globalGroup);
      int procRanges[1][3];
      procRanges[0][0] = firstClusterProcRank;
      procRanges[0][1] = lastClusterProcRank;
      procRanges[0][2] = 1;
      MPI_Group clusterGroup = MPI_GROUP_NULL;
      MPI_Group_range_incl(globalGroup,1,&procRanges[0],&clusterGroup);
      MPI_Comm clusterComm = MPI_COMM_NULL;
      MPI_Comm_create(MPI_COMM_WORLD,clusterGroup,&clusterComm);
      mpiComm = clusterComm;
    }
#endif
    
    //
    // Create the Thyra::ModelEvaluator object
    //
    
    *out << "\nCreate the GLpApp::AdvDiffReactOptModel wrapper object ...\n";
    
#ifdef HAVE_MPI
    Epetra_MpiComm comm(mpiComm);
#else
    Epetra_SerialComm comm;
#endif
    
    GLpApp::AdvDiffReactOptModel
      epetraModel(
        Teuchos::rcp(&comm,false),beta,len_x,len_y,local_nx,local_ny,geomFileBase.c_str()
        ,np,x0,p0,reactionRate,normalizeBasis
        );
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
      if(lowsfParamsFile.length()) {
        Teuchos::updateParametersFromXmlFile(lowsfParamsFile,&*lowsfPL);
        *journalOut << "\nLOWSF parameters read from the file \""<<lowsfParamsFile<<"\":\n";
        lowsfPL->print(*OSTab(journalOut).getOStream(),0,true);
      }
      if(lowsfExtraParams.length()) {
        Teuchos::updateParametersFromXmlString(lowsfExtraParams,&*lowsfPL);
        *journalOut << "\nUpdated with extra LOWSF parameters taken from the command-line:\n";
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
    
    Teuchos::RefCountPtr<Thyra::EpetraModelEvaluator>
      epetraThyraModel = rcp(new Thyra::EpetraModelEvaluator()); // Sets default options!
    epetraThyraModel->setOStream(journalOut);
    epetraThyraModel->initialize(Teuchos::rcp(&epetraModel,false),lowsFactory);

    *out
      << "\nnx = " << epetraThyraModel->get_x_space()->dim()
      << "\nnp = " << epetraThyraModel->get_p_space(0)->dim() << "\n";

    if(matchingVecFile != "") {
      *out << "\nReading the matching vector \'q\' from the file(s) with base name \""<<matchingVecFile<<"\" ...\n";
      epetraModel.set_q(
        Thyra::get_Epetra_Vector(
          *epetraModel.get_x_map(),readVectorFromFile(matchingVecFile,epetraThyraModel->get_x_space()
            )
          )
        );
    }
    
    if(1) {
      Thyra::ModelEvaluatorBase::InArgs<Scalar> epetraThyraModel_initialGuess = epetraThyraModel->createInArgs();
      if(stateGuessFileBase != "") {
        *out << "\nReading the guess of the state \'x\' from the file(s) with base name \""<<stateGuessFileBase<<"\" ...\n";
        epetraThyraModel_initialGuess.set_x(readVectorFromFile(stateGuessFileBase,epetraThyraModel->get_x_space(),scaleStateGuess));
      }
      if(paramGuessFileBase != "") {
        *out << "\nReading the guess of the parameters \'p\' from the file(s) with base name \""<<paramGuessFileBase<<"\" ...\n";
        epetraThyraModel_initialGuess.set_p(0,readVectorFromFile(paramGuessFileBase,epetraThyraModel->get_p_space(0),scaleParamGuess));
      }
      epetraThyraModel->setInitialGuess(epetraThyraModel_initialGuess);
    }

    Teuchos::RefCountPtr<std::ostream>
      modelEvalLogOut;
    if(procRank==0)
      modelEvalLogOut = rcp(new std::ofstream("ModelEvaluationLog.out"));
    else
      modelEvalLogOut = rcp(new Teuchos::oblackholestream());
    Teuchos::RefCountPtr<Thyra::DefaultEvaluationLoggerModelEvaluator<Scalar> >
      loggerThyraModel
      = rcp(
        new Thyra::DefaultEvaluationLoggerModelEvaluator<Scalar>(
          epetraThyraModel,modelEvalLogOut
          )
        );

    if( numClusters > 0 ) {

      Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> >
        x0 = epetraThyraModel->getNominalValues().get_x();
      double nrm_x0;
      
      *out << "\nTiming a global reduction across just this cluster: ||x0||_1 = ";
      timer.start(true);
      nrm_x0 = Thyra::norm_1(*x0);
      *out << nrm_x0 << "\n";
      timer.stop();
      *out << "\n    time = " << timer.totalElapsedTime() << " seconds\n";
      
      *out << "\nTiming a global reduction across the entire set of processors: ||x0||_1 = ";
      timer.start(true);
      RTOpPack::ROpNorm1<Scalar> norm_1_op;
      Teuchos::RefCountPtr<RTOpPack::ReductTarget> norm_1_targ = norm_1_op.reduct_obj_create();
      const Thyra::VectorBase<Scalar>* vecs[] = { &*x0 };
      Teuchos::dyn_cast<const Thyra::MPIVectorBase<Scalar> >(*x0).applyOp(
        MPI_COMM_WORLD,norm_1_op,1,vecs,0,static_cast<Thyra::VectorBase<Scalar>**>(NULL),&*norm_1_targ
        ,0,0,0
        );
      nrm_x0 = norm_1_op(*norm_1_targ);
      *out << nrm_x0 << "\n";
      timer.stop();
      *out << "\n    time = " << timer.totalElapsedTime() << " seconds\n";
      
    }
    
    //
    // Create the NLP
    //
    
    Teuchos::RefCountPtr<NLP> nlp;

    if(do_sim) {
      nlp = rcp(new NLPFirstOrderThyraModelEvaluator(loggerThyraModel,-1,-1));
    }
    else {
      // Setup finite difference object
      RefCountPtr<Thyra::DirectionalFiniteDiffCalculator<Scalar> > direcFiniteDiffCalculator;
      if(use_finite_diff) {
        direcFiniteDiffCalculator = rcp(new Thyra::DirectionalFiniteDiffCalculator<Scalar>());
        direcFiniteDiffCalculator->fd_method_type(fd_method_type);
        direcFiniteDiffCalculator->fd_step_size(fd_step_len);
      }
      if( use_black_box ) {
        // Create a Thyra::NonlinearSolverBase object to solve and eliminate the
        // state variables and the state equations
        Teuchos::RefCountPtr<Thyra::DampenedNewtonNonlinearSolver<Scalar> >
          stateSolver = rcp(new Thyra::DampenedNewtonNonlinearSolver<Scalar>()); // ToDo: Replace with MOOCHO!
        stateSolver->defaultTol(fwd_newton_tol);
        stateSolver->defaultMaxNewtonIterations(fwd_newton_max_iters);
        // Create the reduced  Thyra::ModelEvaluator object for p -> g_hat(p)
        Teuchos::RefCountPtr<Thyra::DefaultStateEliminationModelEvaluator<Scalar> >
          reducedThyraModel = rcp(new Thyra::DefaultStateEliminationModelEvaluator<Scalar>(loggerThyraModel,stateSolver));
        Teuchos::RefCountPtr<Thyra::ModelEvaluator<Scalar> >
          finalReducedThyraModel;
        if(use_finite_diff) {
          // Create the finite-difference wrapped Thyra::ModelEvaluator object
          Teuchos::RefCountPtr<Thyra::DefaultFiniteDifferenceModelEvaluator<Scalar> >
            fdReducedThyraModel = rcp(
              new Thyra::DefaultFiniteDifferenceModelEvaluator<Scalar>(
                reducedThyraModel,direcFiniteDiffCalculator)
              );
          finalReducedThyraModel = fdReducedThyraModel;
        }
        else {
          finalReducedThyraModel = reducedThyraModel;
        }
        // Wrap the reduced NAND Thyra::ModelEvaluator object in an NLP object
        nlp = rcp(new NLPFirstOrderThyraModelEvaluator(finalReducedThyraModel,0,0));
      }
      else {
        if(use_direct) {
          Teuchos::RefCountPtr<NLPDirectThyraModelEvaluator>
            nlpDirect = rcp(new NLPDirectThyraModelEvaluator(loggerThyraModel,0,0));
          if(use_finite_diff) {
            nlpDirect->set_direcFiniteDiffCalculator(direcFiniteDiffCalculator);
          }
          nlp = nlpDirect;
        }
        else {
          nlp = rcp(new NLPFirstOrderThyraModelEvaluator(loggerThyraModel,0,0));
        }
      }
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

    if( stateSoluFileBase != "" && epetraThyraModel->getFinalPoint().get_x().get() ) {
      *out << "\nWriting the state solution \'x\' to the file(s) with base name \""<<stateSoluFileBase<<"\" ...\n";
      writeVectorToFile(*epetraThyraModel->getFinalPoint().get_x(),stateSoluFileBase);
    }
    if( paramSoluFileBase != "" ) {
      *out << "\nWriting the parameter solution \'p\' to the file(s) with base name \""<<paramSoluFileBase<<"\" ...\n";
      writeVectorToFile(*epetraThyraModel->getFinalPoint().get_p(0),paramSoluFileBase);
    }

    // Write the LOWSF parameters that were used:
#if defined(HAVE_TEUCHOS_EXTENDED) && defined(HAVE_TEUCHOS_EXPAT)
    if(lowsfParamsUsedFile != "" && procRank == 0) {
      Teuchos::writeParameterListToXmlFile(*lowsfPL,lowsfParamsUsedFile);
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
