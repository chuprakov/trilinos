// //////////////////////////////////////////////////////////////////
// NLPThyraEpetraModelEval4DOptMain.cpp

#include "NLPInterfacePack_NLPFirstOrderThyraModelEvaluator.hpp"
#include "EpetraModelEval4DOpt.hpp"
#include "Thyra_DefaultModelEvaluatorWithSolveFactory.hpp"
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
  using NLPInterfacePack::NLPFirstOrderThyraModelEvaluator;
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
    bool         externalFactory = false;

    CommandLineProcessor  clp(false); // Don't throw exceptions

    clp.setOption( "xt0", &xt0 );
    clp.setOption( "xt1", &xt1 );
    clp.setOption( "pt0", &pt0 );
    clp.setOption( "pt1", &pt1 );
    clp.setOption( "d", &d );
    clp.setOption( "x00", &x00 );
    clp.setOption( "x01", &x01 );
    clp.setOption( "p00", &p00 );
    clp.setOption( "p01", &p01 );
    clp.setOption( "do-sim", "do-opt",  &do_sim, "Flag for if only the square constraints are solved" );
    clp.setOption( "external-lowsf", "internal-lowsf", &externalFactory
                   ,"Determines of the Thyra::LinearOpWithSolveFactory is used externally or internally to the Thyra::EpetraModelEvaluator object"  );
  
    CommandLineProcessor::EParseCommandLineReturn
      parse_return = clp.parse(argc,argv,&std::cerr);

    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
      return parse_return;

    //
    // Create the NLP
    //
    
    // Create the EpetraExt::ModelEvaluator object

    Teuchos::RefCountPtr<EpetraExt::ModelEvaluator>
      epetraModel = rcp(new EpetraModelEval4DOpt(xt0,xt1,pt0,pt1,d,x00,x01,p00,p01));

    // Create the Thyra::EpetraModelEvaluator object

    Teuchos::RefCountPtr<Thyra::LinearOpWithSolveFactoryBase<double> >
      lowsFactory = rcp(new Thyra::AmesosLinearOpWithSolveFactory());

    Teuchos::RefCountPtr<Thyra::EpetraModelEvaluator>
      epetraThyraModel = rcp(new Thyra::EpetraModelEvaluator());
    
    Teuchos::RefCountPtr<Thyra::ModelEvaluator<double> > thyraModel;
    if(externalFactory) {
      epetraThyraModel->initialize(epetraModel,Teuchos::null);
      thyraModel = Teuchos::rcp(
        new Thyra::DefaultModelEvaluatorWithSolveFactory<double>(
          epetraThyraModel
          ,lowsFactory
          )
        );
    }
    else {
      epetraThyraModel->initialize(epetraModel,lowsFactory);
      thyraModel = epetraThyraModel;
    }

    // ToDo: Specify an initial guess not from the Thyra::EpetraModelEvaluator object
    
    NLPFirstOrderThyraModelEvaluator nlp(
      thyraModel
      ,do_sim ? -1 : 0
      ,do_sim ? -1 : 0
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
