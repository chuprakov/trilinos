// ///////////////////////////////
// cxx_main.cpp

#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFCoreEpetraLinearOp.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_LocalMap.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Time.hpp"
#ifdef RTOp_USE_MPI
#  include "TSFCoreSimpleMPIVectorSpace.hpp"
#  include "Epetra_MpiComm.h"
#endif

// Define this if you want to see only Epetra-based computations
//#define EPETRA_ADAPTERS_EPETRA_ONLY

namespace TSFCore {

///
/** Testing program for TSFCore/Epetra adpaters.
 *
 * This testing program shows how you can easily mix and match
 * different implementations of vectors and multi-vectors for serial
 * and SPMD MPI implementations.  This code is worth study to show how
 * this is done.
 */
int main_body( int argc, char* argv[] ) {

	typedef double Scalar;

	using Teuchos::CommandLineProcessor;
	using Teuchos::RefCountPtr;
	using Teuchos::rcp;
	using Teuchos::rcp_static_cast;
	using Teuchos::rcp_const_cast;
	using Teuchos::set_extra_data;
	
	bool verbose = true;
#ifdef RTOp_USE_MPI
	bool useMPI  = true;
#else
	bool useMPI  = false;
#endif
	bool success = true;
	bool result;

	const Scalar err_tol = 1e-10; // Todo: Make this adjustable!

	MPI_Init(&argc,&argv);
	
	try {

		Teuchos::Time timer("");

		//
		// Read options from the commandline
		//

		int   local_dim      = 4;
		int   num_mv_cols    = 4;

		CommandLineProcessor  clp(false); // Don't throw exceptions
		clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
		clp.setOption( "local-dim", &local_dim, "Number of vector elements per process." );
		clp.setOption( "num-mv-cols", &num_mv_cols, "Number columns in each multi-vector (>=4)." );
		CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFULL ) return parse_return;

		TEST_FOR_EXCEPTION(
			num_mv_cols < 4, std::logic_error
			,"Error, num-mv-cols must be >= 4!"
			);

		//
		// Create two different vector spaces (one Epetra and one non-Epetra)
		// that should be compatible
		//
		MPI_Comm mpiComm = MPI_COMM_NULL;
		int numProc = 1;
#ifdef RTOp_USE_MPI
		mpiComm = MPI_COMM_WORLD;
		MPI_Comm_size( mpiComm, &numProc );
#endif
		RefCountPtr<const Epetra_Comm> epetra_comm;
		RefCountPtr<const Epetra_Map> epetra_map;
		RefCountPtr<const VectorSpace<Scalar> > epetra_vs;
#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
		RefCountPtr<const VectorSpace<Scalar> > non_epetra_vs;
#endif
#ifdef RTOp_USE_MPI
		if(useMPI) {
			//
			// Create parallel vector spaces using compatible maps
			//
			epetra_comm = rcp(new Epetra_MpiComm(mpiComm));
			epetra_map = rcp(new Epetra_Map(-1,local_dim,0,*epetra_comm));
			epetra_vs = rcp(new EpetraVectorSpace(epetra_map));
#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
			// Non-Epetra vector space
			non_epetra_vs = rcp(new SimpleMPIVectorSpace<Scalar>(mpiComm,local_dim));
#endif
		}
		else {
#endif
			//
			// Create serial vector spaces (i.e. VectorSpace::isInCore()==true)
			//
			// Epetra vector space
			epetra_comm = rcp(new Epetra_SerialComm);
			epetra_map = rcp(new Epetra_LocalMap(local_dim,0,*epetra_comm));
			epetra_vs = rcp(new EpetraVectorSpace(epetra_map));
#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
			// Non-Epetra vector space
			non_epetra_vs = rcp(new SerialVectorSpace<Scalar>(local_dim));
#endif
#ifdef RTOp_USE_MPI
		}
#endif

		//
		// Create vectors and multi-vectors from each vector space
		//

		RefCountPtr<Vector<Scalar> >
			ev1 = epetra_vs->createMember(),
			ev2 = epetra_vs->createMember();
#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
		RefCountPtr<Vector<Scalar> >
			nev1 = non_epetra_vs->createMember(),
			nev2 = non_epetra_vs->createMember();
#endif

		RefCountPtr<MultiVector<Scalar> >
			eV1 = epetra_vs->createMembers(num_mv_cols),
			eV2 = epetra_vs->createMembers(num_mv_cols);
#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
		RefCountPtr<MultiVector<Scalar> >
			neV1 = non_epetra_vs->createMembers(num_mv_cols),
			neV2 = non_epetra_vs->createMembers(num_mv_cols);
#endif

		//
		// Check for compatibility of the vector and Multi-vectors
		// w.r.t. RTOps
		//

		if(verbose)
			std::cout
				<< "\n***"
				<< "\n*** Testing individual vector/multi-vector RTOps"
				<< "\n***\n";

		assign( &*ev1, 0.0 );
		assign( &*ev2, 1.0 );
#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
		assign( &*nev1, 0.0 );
		assign( &*nev2, 1.0 );
#endif
		assign( &*eV1, 0.0 );
		assign( &*eV2, 1.0 );
#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
		assign( &*neV1, 0.0 );
		assign( &*neV2, 1.0 );
#endif

		Scalar
			ev1_nrm = norm_1(*ev1),
			eV1_nrm = norm_1(*eV1);
#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
		Scalar
			nev1_nrm = norm_1(*nev1),
			neV1_nrm = norm_1(*neV1);
#endif
		
 		result = (ev1_nrm == 0.0);
		if(verbose) std::cout << "\nnorm_1(ev1) == 0 : " << (result ? "passed!" : "failed!") << std::endl;
		if(!result) success = false;

#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
		result = (nev1_nrm == 0.0);
		if(verbose) std::cout << "\nnorm_1(nev1) == 0 : " << (result ? "passed!" : "failed!") << std::endl;
		if(!result) success = false;
#endif
		
 		result = (eV1_nrm == 0.0);
		if(verbose) std::cout << "\nnorm_1(eV1) == 0 : " << (result ? "passed!" : "failed!") << std::endl;
		if(!result) success = false;

#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
		result = (neV1_nrm == 0.0);
		if(verbose) std::cout << "\nnorm_1(neV1) == 0 : " << (result ? "passed!" : "failed!") << std::endl;
		if(!result) success = false;
#endif

		if(verbose)
			std::cout
				<< "\n***"
				<< "\n*** Testing Test RTOps with two or more arguments"
				<< "\n***\n";

 		assign( &*ev1, *ev2 );
		assign( &*eV1, *eV2 );
#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
		assign( &*nev1, *nev2 );
		assign( &*neV1, *neV2 );
#endif

		if(verbose) std::cout << "\nnorm_1(env1) = " << norm_1(*ev1) << std::endl;
		if(verbose) std::cout << "\nnorm_1(eV1) = " << norm_1(*eV1) << std::endl;
#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
		if(verbose) std::cout << "\nnorm_1(nev1) = " << norm_1(*nev1) << std::endl;
		if(verbose) std::cout << "\nnorm_1(neV1) = " << norm_1(*neV1) << std::endl;
#endif

		if(verbose)
			std::cout
				<< "\n***"
				<< "\n*** Test RTOps with two or more arguments with mixed types"
				<< "\n***\n";

#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
 		assign( &*ev1, *nev2 );
		assign( &*eV1, *neV2 );
		assign( &*nev1, *ev2 );
		assign( &*neV1, *eV2 );
#endif

		if(verbose) std::cout << "\nnorm_1(env1) = " << norm_1(*ev1) << std::endl;
		if(verbose) std::cout << "\nnorm_1(eV1) = " << norm_1(*eV1) << std::endl;
#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
		if(verbose) std::cout << "\nnorm_1(neV1) = " << norm_1(*neV1) << std::endl;
		if(verbose) std::cout << "\nnorm_1(nev1) = " << norm_1(*nev1) << std::endl;
#endif

		if(verbose)
			std::cout
				<< "\n***"
				<< "\n*** Test MultiVector::apply(...)"
				<< "\n***\n";

		RefCountPtr<MultiVector<Scalar> >
			T = eV1->domain()->createMembers(num_mv_cols);

		if(verbose) std::cout << "\nPerforming eV1'*eV2 ...\n";
		timer.start();
		eV1->apply( TRANS, *eV2, &*T );
		timer.stop();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(verbose) std::cout << "  norm_1(eV1'*eV2) = " << norm_1(*T) << std::endl;

#ifndef EPETRA_ADAPTERS_EPETRA_ONLY

		if(verbose) std::cout << "\nPerforming neV1'*eV2 ...\n";
		timer.start();
		neV1->apply( TRANS, *eV2, &*T );
		timer.stop();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(verbose) std::cout << "  norm_1(neV1'*eV2) = " << norm_1(*T) << std::endl;

		if(verbose) std::cout << "\nPerforming eV1'*neV2 ...\n";
		timer.start();
		eV1->apply( TRANS, *neV2, &*T );
		timer.stop();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(verbose) std::cout << "  norm_1(eV1'*neV2) = " << norm_1(*T) << std::endl;

		if(verbose) std::cout << "\nPerforming neV1'*neV2 ...\n";
		timer.start();
		neV1->apply( TRANS, *neV2, &*T );
		timer.stop();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(verbose) std::cout << "  norm_1(neV1'*neV2) = " << norm_1(*T) << std::endl;

#endif

		if(verbose)
			std::cout
				<< "\n***"
				<< "\n*** Creating a diagonal Epetra_Operator"
				<< "\n***\n";

		RefCountPtr<Epetra_Operator>  epetra_op;

		if(1) {
			// Create a diagonal matrix with 0.5 on the diagonal
			RefCountPtr<Epetra_CrsMatrix>
				epetra_mat = rcp(new Epetra_CrsMatrix(::Copy,*epetra_map,1));
			double values[1] = { 0.5 };
			int indices[1];
			const int IB = epetra_map->IndexBase();
			for( int k = 0; k < local_dim; ++k ) {
				indices[0] = k + IB;
				epetra_mat->InsertGlobalValues(
					k+IB                // GlobalRow
					,1                  // NumEntries
					,values             // Values
					,indices            // Indices
					);
			}
			epetra_mat->FillComplete();
			epetra_op = epetra_mat;
		}

		RefCountPtr<const LinearOp<Scalar> >
			Op = rcp(new EpetraLinearOp(epetra_op));
		
		if(verbose)
			std::cout
				<< "\n***"
				<< "\n*** Mix and match vector and Multi-vectors with Epetra opeator"
				<< "\n***\n";
		
		RefCountPtr<Vector<Scalar> >
			ey  = epetra_vs->createMember();
		RefCountPtr<MultiVector<Scalar> >
			eY  = epetra_vs->createMembers(num_mv_cols);
#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
		RefCountPtr<Vector<Scalar> >
			ney = non_epetra_vs->createMember();
		RefCountPtr<MultiVector<Scalar> >
			neY = non_epetra_vs->createMembers(num_mv_cols);
#endif

		if(verbose) std::cout << "\nPerforming ey = 2*Op*ev1 ...\n";
		timer.start();
		Op->apply( NOTRANS, *ev1, &*ey, 2.0 );
		timer.stop();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(verbose) std::cout << "  norm_1(ey) = " << norm_1(*ey) << std::endl;

		if(verbose) std::cout << "\nPerforming eY = 2*Op*eV1 ...\n";
		timer.start();
		Op->apply( NOTRANS, *eV1, &*eY, 2.0 );
		timer.stop();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(verbose) std::cout << "  norm_1(eY) = " << norm_1(*eY) << std::endl;

#ifndef EPETRA_ADAPTERS_EPETRA_ONLY

		if(verbose) std::cout << "\nPerforming ney = 2*Op*ev1 ...\n";
		timer.start();
		Op->apply( NOTRANS, *ev1, &*ney, 2.0 );
		timer.stop();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(verbose) std::cout << "  norm_1(ney) = " << norm_1(*ney) << std::endl;

		if(verbose) std::cout << "\nPerforming neY = 2*Op*eV1 ...\n";
		timer.start();
		Op->apply( NOTRANS, *eV1, &*neY, 2.0 );
		timer.stop();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(verbose) std::cout << "  norm_1(neY) = " << norm_1(*neY) << std::endl;

		if(verbose) std::cout << "\nPerforming ey = 2*Op*nev1 ...\n";
		timer.start();
		Op->apply( NOTRANS, *nev1, &*ey, 2.0 );
		timer.stop();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(verbose) std::cout << "  norm_1(ey) = " << norm_1(*ey) << std::endl;

		if(verbose) std::cout << "\nPerforming eY = 2*Op*neV1 ...\n";
		timer.start();
		Op->apply( NOTRANS, *neV1, &*eY, 2.0 );
		timer.stop();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(verbose) std::cout << "  norm_1(eY) = " << norm_1(*eY) << std::endl;

		if(verbose) std::cout << "\nPerforming ney = 2*Op*nev1 ...\n";
		timer.start();
		Op->apply( NOTRANS, *nev1, &*ney, 2.0 );
		timer.stop();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(verbose) std::cout << "  norm_1(ney) = " << norm_1(*ney) << std::endl;

		if(verbose) std::cout << "\nPerforming neY = 2*Op*neV1 ...\n";
		timer.start();
		Op->apply( NOTRANS, *neV1, &*neY, 2.0 );
		timer.stop();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(verbose) std::cout << "  norm_1(neY) = " << norm_1(*neY) << std::endl;

#endif

		if(verbose)
			std::cout
				<< "\n***"
				<< "\n*** Testing Multi-vector views with Epetra operator"
				<< "\n***\n";

		const Range1D col_rng(1,2);
		const int numCols = 2;
		const int cols[] = { 3, 4 };

		RefCountPtr<const MultiVector<Scalar> >
			eV1_v1  = rcp_static_cast<const MultiVector<Scalar> >(eV1)->subView(col_rng),
			eV1_v2  = rcp_static_cast<const MultiVector<Scalar> >(eV1)->subView(numCols,cols);
#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
		RefCountPtr<const MultiVector<Scalar> >
			neV1_v1  = rcp_static_cast<const MultiVector<Scalar> >(neV1)->subView(col_rng),
			neV1_v2  = rcp_static_cast<const MultiVector<Scalar> >(neV1)->subView(numCols,cols);
#endif

		if(verbose) std::cout << "\nPerforming eY_v1 = 2*Op*eV1_v1 ...\n";
		timer.start();
		Op->apply( NOTRANS, *eV1_v1, &*eY->subView(col_rng), 2.0 );
		timer.stop();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(verbose) std::cout << "  norm_1(eY) = " << norm_1(*eY->subView(col_rng)) << std::endl;

		if(verbose) std::cout << "\nPerforming eY_v2 = 2*Op*eV1_v2 ...\n";
		timer.start();
		Op->apply( NOTRANS, *eV1_v2, &*eY->subView(numCols,cols), 2.0 );
		timer.stop();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(verbose) std::cout << "  norm_1(eY) = " << norm_1(*eY->subView(numCols,cols)) << std::endl;

#ifndef EPETRA_ADAPTERS_EPETRA_ONLY

		if(verbose) std::cout << "\nPerforming neY_v1 = 2*Op*eV1_v1 ...\n";
		timer.start();
		Op->apply( NOTRANS, *eV1_v1, &*neY->subView(col_rng), 2.0 );
		timer.stop();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(verbose) std::cout << "  norm_1(eY) = " << norm_1(*neY->subView(col_rng)) << std::endl;

		if(verbose) std::cout << "\nPerforming eY_v1 = 2*Op*neV1_v1 ...\n";
		timer.start();
		Op->apply( NOTRANS, *neV1_v1, &*eY->subView(col_rng), 2.0 );
		timer.stop();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(verbose) std::cout << "  norm_1(eY) = " << norm_1(*eY->subView(col_rng)) << std::endl;

		if(verbose) std::cout << "\nPerforming neY_v2 = 2*Op*eV1_v2 ...\n";
		timer.start();
		Op->apply( NOTRANS, *eV1_v2, &*neY->subView(numCols,cols), 2.0 );
		timer.stop();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(verbose) std::cout << "  norm_1(eY) = " << norm_1(*neY->subView(numCols,cols)) << std::endl;

		if(verbose) std::cout << "\nPerforming eY_v2 = 2*Op*neV1_v2 ...\n";
		timer.start();
		Op->apply( NOTRANS, *neV1_v2, &*eY->subView(numCols,cols), 2.0 );
		timer.stop();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(verbose) std::cout << "  norm_1(eY) = " << norm_1(*eY->subView(numCols,cols)) << std::endl;

#endif

		if(verbose)
			std::cout << "\nAll of the tests seem to have run sucessfully!\n";

	} // end try
	catch( const std::exception &excpt ) {
		if(verbose)
			std::cerr << "*** Caught a standard exception : " << excpt.what() << std::endl;
		success = -1;
	}
	catch( ... ) {
		if(verbose)
			std::cerr << "*** Caught an unknown exception!\n";
		success = -1;
	}

 	MPI_Finalize();

	return (success ? 0 : -1);

}

} // namespace TSFCore

int main( int argc, char* argv[] ) {
	return TSFCore::main_body(argc,argv);
}
