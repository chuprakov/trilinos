#include "GLpApp_AdvDiffReactOptModelCreator.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_ParallelMultiVectorFileIO.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Thyra_RealLinearOpWithSolveFactoryCreator.hpp"
#include "MoochoPack_ThyraModelEvaluatorSolver.hpp"
#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif

namespace {

typedef AbstractLinAlgPack::value_type  Scalar;

} // namespace

int main( int argc, char* argv[] )
{
  using Teuchos::rcp;
  using Teuchos::RefCountPtr;
  using Teuchos::OSTab;
  using MoochoPack::MoochoSolver;
  using MoochoPack::ThyraModelEvaluatorSolver;
  using Teuchos::CommandLineProcessor;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  const int procRank = mpiSession.getRank();
  const int numProcs = mpiSession.getNProc();

  Teuchos::Time timer("");
  
  bool dummySuccess = true;

  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {
  
    // Create the solver object
    Thyra::RealLinearOpWithSolveFactoryCreator lowsfCreator;
    GLpApp::AdvDiffReactOptModelCreator epetraModelCreator;
    ThyraModelEvaluatorSolver solver;

    //
    // Get options from the command line
    //

    std::string         matchingVecFile     = "";
    bool                dump_all            = false;
    int                 numProcsPerCluster  = -1;

    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    epetraModelCreator.setupCLP(&clp);
    lowsfCreator.setupCLP(&clp);
    solver.setupCLP(&clp);
    clp.setOption( "dump-all", "no-dump-all",  &dump_all, "Flag for if we dump everything to STDOUT" );
    clp.setOption( "q-vec-file", &matchingVecFile, "Base file name to read the objective state matching vector q (i.e. ||x-q||_M in objective)." );
    clp.setOption( "num-procs-per-cluster", &numProcsPerCluster, "Number of processors in a cluster (<=0 means only one cluster)." );

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

    Teuchos::RefCountPtr<Epetra_Comm> comm = Teuchos::null;
#ifdef HAVE_MPI
    comm = Teuchos::rcp(new Epetra_MpiComm(mpiComm));
#else
    comm = Teuchos::rcp(new Epetra_SerialComm());
#endif
    
    //
    // Create the Thyra::ModelEvaluator object
    //
    
    *out << "\nCreate the GLpApp::AdvDiffReactOptModel wrapper object ...\n";
    
    Teuchos::RefCountPtr<GLpApp::AdvDiffReactOptModel>
      epetraModel = epetraModelCreator.createModel(comm);
    epetraModel->setOStream(journalOut);
    if(dump_all) epetraModel->setVerbLevel(Teuchos::VERB_EXTREME);

    *out << "\nCreate the Thyra::LinearOpWithSolveFactory object ...\n";

    Teuchos::RefCountPtr<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
      lowsFactory = lowsfCreator.createLOWSF(OSTab(journalOut).getOStream().get());
    
    *out << "\nCreate the Thyra::EpetraModelEvaluator wrapper object ...\n";
    
    Teuchos::RefCountPtr<Thyra::EpetraModelEvaluator>
      epetraThyraModel = rcp(new Thyra::EpetraModelEvaluator()); // Sets default options!
    epetraThyraModel->setOStream(journalOut);
    epetraThyraModel->initialize(epetraModel,lowsFactory);

    *out
      << "\nnx = " << epetraThyraModel->get_x_space()->dim()
      << "\nnp = " << epetraThyraModel->get_p_space(0)->dim() << "\n";

    if(matchingVecFile != "") {
      *out << "\nReading the matching vector \'q\' from the file(s) with base name \""<<matchingVecFile<<"\" ...\n";
      Thyra::ParallelMultiVectorFileIO<Scalar> fileIO;
      epetraModel->set_q(
        Thyra::get_Epetra_Vector(
          *epetraModel->get_x_map(),fileIO.readVectorFromFile(matchingVecFile,epetraThyraModel->get_x_space())
          )
        );
    }

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
    // Solve the NLP
    //

    // Set the journal file
    solver.getSolver().set_journal_out(journalOut);
    
    // Set the model
    solver.setModel(epetraThyraModel);

    solver.readInitialGuess(out.get());

    // Solve the NLP
    const MoochoSolver::ESolutionStatus	solution_status = solver.solve();

    solver.writeFinalSolution(out.get());
    
    //
    // Return the solution status (0 if sucessfull)
    //

    return solution_status;

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,dummySuccess)

  return MoochoSolver::SOLVE_RETURN_EXCEPTION;
}
