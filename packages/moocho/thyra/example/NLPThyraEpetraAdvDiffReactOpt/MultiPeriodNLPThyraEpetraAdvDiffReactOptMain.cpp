#include "GLpApp_AdvDiffReactOptModelCreator.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_ParallelMultiVectorFileIO.hpp"
#include "Thyra_DefaultClusteredMPIProductVectorSpace.hpp"
#include "Thyra_DefaultMultiPeriodModelEvaluator.hpp"
#include "Thyra_VectorSpaceTester.hpp"
#include "RTOpPack_MPI_apply_op_decl.hpp"
#include "Teuchos_OpaqueWrapper.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_arrayArg.hpp"
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
  using Teuchos::OpaqueWrapper;
  using Teuchos::OSTab;
  using MoochoPack::MoochoSolver;
  using MoochoPack::ThyraModelEvaluatorSolver;
  using Teuchos::CommandLineProcessor;
  typedef Thyra::ModelEvaluatorBase MEB;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  const int procRank = mpiSession.getRank();
  const int numProcs = mpiSession.getNProc();

  Teuchos::Time timer("");
  
  bool dummySuccess = true;

  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {
  
    Thyra::RealLinearOpWithSolveFactoryCreator lowsfCreator;
    GLpApp::AdvDiffReactOptModelCreator epetraModelCreator;

    // Create the solver object
    ThyraModelEvaluatorSolver solver;

    //
    // Get options from the command line
    //

    int                 numProcsPerCluster     = -1;
    double              perturbedParamScaling  = 1.0;
    bool                dumpAll                = false;

    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    epetraModelCreator.setupCLP(&clp);
    lowsfCreator.setupCLP(&clp);
    solver.setupCLP(&clp);
    clp.setOption( "num-procs-per-cluster", &numProcsPerCluster, "Number of processes in a cluster (<=0 means only one cluster)." );
    clp.setOption( "p-perturb-scaling", &perturbedParamScaling, "Scaling for perturbed paramters from the initial forward solve." );
    clp.setOption( "dump-all", "no-dump-all", &dumpAll, "Set to true, then a bunch of debugging output will be created." );

    CommandLineProcessor::EParseCommandLineReturn
      parse_return = clp.parse(argc,argv,&std::cerr);

    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
      return parse_return;

    *out
      << "\n***"
      << "\n*** NLPThyraEpetraAdvDiffReactOptMain, Global numProcs = "<<numProcs
      << "\n***\n";

    int clusterRank = -1;
    int numClusters = -1;
#ifdef HAVE_MPI
    RefCountPtr<OpaqueWrapper<MPI_Comm> >
      intraClusterComm = Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD),
      interClusterComm = Teuchos::null;
    if( numProcsPerCluster > 0 ) {
      *out << "\nCreating communicator for local cluster of "<<numProcsPerCluster<<" processes ...\n";
      numClusters = numProcs/numProcsPerCluster;
      const int remainingProcs = numProcs%numProcsPerCluster;
      TEST_FOR_EXCEPTION(
        remainingProcs!=0,std::logic_error
        ,"Error, The number of processes per cluster numProcsPerCluster="<<numProcsPerCluster
        << " does not divide into the global number of processes numProcs="<<numProcs
        << " and instead has remainder="<<remainingProcs<<"!"
        );
      // Determine which cluster this process is part of and what the global
      // process ranges are.
      clusterRank = procRank / numProcsPerCluster; // Integer division!
      *out << "\nclusterRank = " << clusterRank << "\n";
      const int firstClusterProcRank = clusterRank * numProcsPerCluster;
      const int lastClusterProcRank = firstClusterProcRank + numProcsPerCluster - 1;
      *out << "\nclusterProcRange = ["<<firstClusterProcRank<<","<<lastClusterProcRank<<"]\n";
      // Create the communicator for this cluster of processes
      *out << "\nCreating intraClusterComm ...";
      MPI_Comm rawIntraClusterComm = MPI_COMM_NULL;
      MPI_Comm_split(
        MPI_COMM_WORLD        // comm
        ,clusterRank          // color (will all be put in the same output comm)
        ,0                    // key (not important here)
        ,&rawIntraClusterComm // newcomm
        );
      intraClusterComm = Teuchos::opaqueWrapper(rawIntraClusterComm,MPI_Comm_free);
      if(1) {
        *out << "\nintraClusterComm:";
        Teuchos::OSTab tab(out);
        int rank, size;
        MPI_Comm_size(*intraClusterComm,&size);
        MPI_Comm_rank(*intraClusterComm,&rank);
        *out << "\nsize="<<size;
        *out << "\nrank="<<rank;
        *out << "\n";
      }
      // Create the communicator for just the root process in each cluster
      *out << "\nCreating interClusterComm ...";
      MPI_Comm rawInterClusterComm = MPI_COMM_NULL;
      MPI_Comm_split(
        MPI_COMM_WORLD                                  // comm
        ,procRank==firstClusterProcRank?0:MPI_UNDEFINED // color
        ,0                                              // key
        ,&rawInterClusterComm                           // newcomm
        );
      if(rawInterClusterComm!=MPI_COMM_NULL)
        interClusterComm = Teuchos::opaqueWrapper(rawInterClusterComm,MPI_Comm_free);
      else
        interClusterComm = Teuchos::opaqueWrapper(rawInterClusterComm);
      if(1) {
        *out << "\ninterClusterComm:";
        Teuchos::OSTab tab(out);
        if(*interClusterComm==MPI_COMM_NULL) {
          *out << " NULL\n";
        }
        else {
          int rank, size;
          MPI_Comm_size(*interClusterComm,&size);
          MPI_Comm_rank(*interClusterComm,&rank);
          *out << "\nsize="<<size;
          *out << "\nrank="<<rank;
          *out << "\n";
        }
      }
    }
#endif

    Teuchos::RefCountPtr<Epetra_Comm> comm = Teuchos::null;
#ifdef HAVE_MPI
    comm = Teuchos::rcp(new Epetra_MpiComm(*intraClusterComm));
    Teuchos::set_extra_data(intraClusterComm,"mpiComm",&comm);
#else
    comm = Teuchos::rcp(new Epetra_SerialComm());
#endif
    
    //
    // Create the Thyra::ModelEvaluator object
    //
    
    *out << "\nCreate the GLpApp::AdvDiffReactOptModel wrapper object ...\n";
    
    Teuchos::RefCountPtr<GLpApp::AdvDiffReactOptModel>
      epetraModel = epetraModelCreator.createModel(comm);

    *out << "\nCreate the Thyra::LinearOpWithSolveFactory object ...\n";

    Teuchos::RefCountPtr<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
      lowsFactory = lowsfCreator.createLOWSF(OSTab(out).getOStream().get());
    
    *out << "\nCreate the Thyra::EpetraModelEvaluator wrapper object ...\n";
    
    Teuchos::RefCountPtr<Thyra::EpetraModelEvaluator>
      epetraThyraModel = rcp(new Thyra::EpetraModelEvaluator()); // Sets default options!
    epetraThyraModel->setOStream(out);
    epetraThyraModel->initialize(epetraModel,lowsFactory);

    *out
      << "\nnx = " << epetraThyraModel->get_x_space()->dim()
      << "\nnp = " << epetraThyraModel->get_p_space(0)->dim() << "\n";

    Teuchos::RefCountPtr<Thyra::ModelEvaluator<Scalar> >
      thyraModel = epetraThyraModel;
    
#ifdef HAVE_MPI

    if( numClusters > 0 ) {

      *out << "\nCreate block parallel vector spaces for multi-period model.x and model.f ...\n";
      Teuchos::RefCountPtr<Thyra::DefaultClusteredMPIProductVectorSpace<Scalar> >
        x_bar_space = Teuchos::rcp(
          new Thyra::DefaultClusteredMPIProductVectorSpace<Scalar>(
            intraClusterComm
            ,0 // clusterRootRank
            ,interClusterComm
            ,1 // numBlocks
            ,Teuchos::arrayArg<Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> > >(
              epetraThyraModel->get_x_space()
              )()
            )
          ),
        f_bar_space = Teuchos::rcp(
          new Thyra::DefaultClusteredMPIProductVectorSpace<Scalar>(
            intraClusterComm
            ,0 // clusterRootRank
            ,interClusterComm
            ,1 // numBlocks
            ,Teuchos::arrayArg<Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> > >(
              epetraThyraModel->get_f_space()
              )()
            )
          );

      Thyra::VectorSpaceTester<Scalar> vectorSpaceTester;
      vectorSpaceTester.show_all_tests(true);
      vectorSpaceTester.dump_all(dumpAll);

      RTOpPack::show_mpi_apply_op_dump = dumpAll;
      Thyra::MPIVectorBase<Scalar>::show_dump = dumpAll;

      *out << "\nTesting the vector space x_bar_space ...\n";
      vectorSpaceTester.check(*x_bar_space,&*OSTab(out).getOStream());

      *out << "\nTesting the vector space f_bar_space ...\n";
      vectorSpaceTester.check(*f_bar_space,&*OSTab(out).getOStream());
      
      Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> >
        x0 = epetraThyraModel->getNominalValues().get_x();
      double nrm_x0;
      
      *out << "\nTiming a global reduction across just this cluster: ||x0||_1 = ";
      timer.start(true);
      nrm_x0 = Thyra::norm_1(*x0);
      *out << nrm_x0 << "\n";
      timer.stop();
      *out << "\n    time = " << timer.totalElapsedTime() << " seconds\n";
      
      *out << "\nTiming a global reduction across the entire set of processes: ||x0||_1 = ";
      timer.start(true);
      RTOpPack::ROpNorm1<Scalar> norm_1_op;
      Teuchos::RefCountPtr<RTOpPack::ReductTarget> norm_1_targ = norm_1_op.reduct_obj_create();
      const Thyra::VectorBase<Scalar>* vecs[] = { &*x0 };
      Teuchos::dyn_cast<const Thyra::MPIVectorBase<Scalar> >(*x0).applyOp(
        MPI_COMM_WORLD,norm_1_op,1,vecs,0,static_cast<Thyra::VectorBase<Scalar>**>(NULL),&*norm_1_targ
        ,0,-1,0
        );
      nrm_x0 = norm_1_op(*norm_1_targ);
      *out << nrm_x0 << "\n";
      timer.stop();
      *out << "\n    time = " << timer.totalElapsedTime() << " seconds\n";

      RTOpPack::show_mpi_apply_op_dump = false;
      Thyra::MPIVectorBase<Scalar>::show_dump = false;

      const int N = 1;
      const int z_index = 1;
      Teuchos::Array<Teuchos::RefCountPtr<Thyra::ModelEvaluator<Scalar> > >
        models(N);
      Teuchos::Array<Scalar>
        weights(N);
      Teuchos::Array<Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > >
        z(N);
      for( int i = 0; i < N; ++i ) {
        models[i] = epetraThyraModel;
        weights[i] = 1.0;
        z[i] = epetraThyraModel->getNominalValues().get_p(z_index)->clone_v();
      }
      thyraModel =
        rcp(
          new Thyra::DefaultMultiPeriodModelEvaluator<Scalar>(
            1,&models[0],&weights[0],&z[0],z_index
            ,x_bar_space,f_bar_space
            )
          );

    }

#endif // HAVE_MPI

    MoochoSolver::ESolutionStatus solution_status;

    //
    *out << "\n***\n*** Solving the initial forward problem\n***\n";
    //

    // Set the deliminator for the output files!
    solver.getSolver().set_output_context("fwd-init");

    // Set the solve mode to solve the forward problem
    solver.setDoSim(true);
    
    // Set the model
    solver.setModel(thyraModel);
    
    // Set the initial guess from files (if specified on commandline)
    solver.readInitialGuess(out.get());
    
    // Solve the initial forward problem
    solution_status = solver.solve();

    // Save the solution for model.x and model.p to be used later
    RefCountPtr<Thyra::VectorBase<Scalar> >
      x_opt = solver.getFinalPoint().get_x()->clone_v(),
      x_init = solver.getFinalPoint().get_x()->clone_v(),
      p_init = solver.getFinalPoint().get_p(0)->clone_v();
    
    //
    *out << "\n***\n*** Solving the perturbed forward problem\n***\n";
    //

    // Set the deliminator for the output files!
    solver.getSolver().set_output_context("fwd");

    // Set the solve mode to solve the forward problem
    solver.setDoSim(true);
   
    // Set the model
    solver.setModel(thyraModel);

    // Set the initial guess and the perturbed parameters
    if(1) {
      MEB::InArgs<Scalar> initialGuess = thyraModel->createInArgs();
      initialGuess.setArgs(thyraModel->getNominalValues());
      initialGuess.set_x(x_init);
      Thyra::Vt_S(&*p_init,perturbedParamScaling);
      initialGuess.set_p(0,p_init);
      //*out << "\nInitial Guess:\n" << Teuchos::describe(initialGuess,Teuchos::VERB_EXTREME);
      solver.setInitialGuess(initialGuess);
    }

    // Solve the perturbed forward problem
    solution_status = solver.solve();

    // Save the solution for model.x and model.p to be used later
    x_init = solver.getFinalPoint().get_x()->clone_v();
    p_init = solver.getFinalPoint().get_p(0)->clone_v();
    
    //
    *out << "\n***\n*** Solving the perturbed inverse problem\n***\n";
    //

    // Set the deliminator for the output files!
    solver.getSolver().set_output_context("inv");

    // Set the matching vector
    epetraModel->set_q(Thyra::get_Epetra_Vector(*epetraModel->get_x_map(),x_opt));

    // Set the solve mode to solve the inverse problem
    solver.setDoSim(false);
   
    // Set the model
    solver.setModel(thyraModel);

    // Set the initial guess for model.x and model.p
    if(1) {
      MEB::InArgs<Scalar> initialGuess = thyraModel->createInArgs();
      initialGuess.setArgs(thyraModel->getNominalValues());
      initialGuess.set_x(x_init);
      initialGuess.set_p(0,p_init);
      //*out << "\nInitial Guess:\n" << Teuchos::describe(initialGuess,Teuchos::VERB_EXTREME);
      solver.setInitialGuess(initialGuess);
    }
    
    // Solve the inverse problem
    solution_status = solver.solve();
    
    //
    // Write the final solution
    //
    
    solver.writeFinalSolution(out.get());
    
    //
    // Return the solution status (0 if successful)
    //
    
    return solution_status;
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,dummySuccess)
    
  return MoochoSolver::SOLVE_RETURN_EXCEPTION;
  
}
