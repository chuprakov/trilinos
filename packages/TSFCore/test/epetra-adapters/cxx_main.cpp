// ///////////////////////////////
// cxx_main.cpp

#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFCoreEpetraLinearOp.hpp"
#include "TSFCoreEpetraMultiVector.hpp"
#include "TSFCoreTestingTools.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_LocalMap.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "dynamic_cast_verbose.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_oblackholestream.hpp"
#ifdef RTOp_USE_MPI
#  include "TSFCoreSimpleMPIVectorSpace.hpp"
#  include "Epetra_MpiComm.h"
#endif

// Define this if you want to see only Epetra-based computations
//#define EPETRA_ADAPTERS_EPETRA_ONLY

// Define this if you want to exclude opeations with Epetra_Operator
//#define EPETRA_ADAPTERS_EXCLUDE_EPETRA_OPERATOR

//
// Some helper functions
//

namespace {

double relErr( const double& v1, const double& v2 )
{
	return std::fabs(v1 - v2) / ( 1.0 + std::fabs(v1) + std::fabs(v2) );
}

bool testRelErr(
	const std::string    &v1_name
	,const double        &v1
	,const std::string   &v2_name
	,const double        &v2
	,const std::string   &maxRelErr_name
	,const double        &maxRelErr
	,bool                verbose
	,std::ostream        &out
	)
{
	const double rel_err = relErr( v1, v2 );
	const bool success = ( rel_err <= maxRelErr );
	if(verbose) {
		out << "\nCheck: rel_err(" << v1_name << "," << v2_name << ")\n"
			<< "       = rel_err(" << v1 << "," << v2 << ") "
			<< "= " << rel_err
			<< " <= " << maxRelErr_name << " = " << maxRelErr << " : "
			<<  (success ? "passed!" : "failed!") << std::endl;
	}
	return success;
}

void print_performance_stats(
	const int        num_time_samples
	,const double    raw_epetra_time
	,const double    tsfcore_wrapped_time
	,bool            verbose
	,std::ostream    &out
	)
{
	if(verbose)
		out << "\nAverage times (out of " << num_time_samples << " samples):\n"
			<< "  Raw Epetra              = " << (raw_epetra_time/num_time_samples) << std::endl
			<< "  TSFCore Wrapped Epetra  = " << (tsfcore_wrapped_time/num_time_samples) << std::endl
			<< "\nRelative performance of TSFCore wrapped verses raw Epetra:\n"
			<< "  ( raw epetra time / tsfcore wrapped time ) = ( " << raw_epetra_time << " / " << tsfcore_wrapped_time << " ) = "
			<< (raw_epetra_time/tsfcore_wrapped_time) << std::endl;
}

} // namespace

namespace TSFCore {

///
/** Testing program for TSFCore/Epetra adpaters.
 *
 * This testing program shows how you can easily mix and match
 * different implementations of vectors and multi-vectors for serial
 * and SPMD MPI implementations.  This code is worth study to show how
 * this is done.
 *
 * Note that the tests performed do not prove that the Epetra adapters
 * (or Epetra itself) perform correctly as only a few post conditions
 * are checked.  Because of the simple nature of these computations it
 * would be possible to put in more exactly component-wise tests if
 * that is needed in the future.
 */
int main_body( int argc, char* argv[] ) {

	typedef double Scalar;

	using DynamicCastHelperPack::dyn_cast;
	using Teuchos::CommandLineProcessor;
	using Teuchos::RefCountPtr;
	using Teuchos::rcp;
	using Teuchos::rcp_static_cast;
	using Teuchos::rcp_const_cast;
	using Teuchos::set_extra_data;
	
	bool verbose = true;
	bool dumpAll = false;
#ifdef RTOp_USE_MPI
	bool useMPI  = true;
#else
	bool useMPI  = false;
#endif
	bool success = true;
	bool result;

	int procRank = 0;

	const Scalar err_tol = 1e-10; // Todo: Make this adjustable!

	Scalar rel_err;

	MPI_Init(&argc,&argv);
	
	try {

		Teuchos::Time timer("");

		//
		// Read options from the commandline
		//

		int     local_dim         = 1000;
		int     num_mv_cols       = 4;
		double  max_rel_err       = 1e-13;
		double  scalar            = 1.5;
		double  max_flop_rate     = 2.0e8;
		CommandLineProcessor  clp(false); // Don't throw exceptions
		clp.setOption( "verbose", "quiet", &verbose, "Determines if any output is printed or not." );
		clp.setOption( "dump-all", "no-dump", &dumpAll, "Determines if quantities are dumped or not." );
		clp.setOption( "local-dim", &local_dim, "Number of vector elements per process." );
		clp.setOption( "num-mv-cols", &num_mv_cols, "Number columns in each multi-vector (>=4)." );
		clp.setOption( "max-rel-err", &max_rel_err, "Maximum relative error for tests." );
		clp.setOption( "scalar", &scalar, "A scalar used in all computations." );
		clp.setOption( "max-flop-rate", &max_flop_rate, "Approx flop rate used for loop timing." );
		CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFULL ) return parse_return;

		TEST_FOR_EXCEPTION(
			num_mv_cols < 4, std::logic_error
			,"Error, num-mv-cols must be >= 4!"
			);

		//
		// Get basic MPI info
		//

		MPI_Comm mpiComm = MPI_COMM_NULL;
		int numProc = 1;
#ifdef RTOp_USE_MPI
		mpiComm = MPI_COMM_WORLD;
		MPI_Comm_size( mpiComm, &numProc );
		MPI_Comm_rank( mpiComm, &procRank );
#endif

		//
		// Setup the output stream
		//
		
		Teuchos::oblackholestream black_hole_out;
		std::ostream &out = ( procRank == 0 ? std::cout : black_hole_out );

		if(verbose)
			out
				<< "\n***"
				<< "\n*** (A) Creating two vector spaces (an Epetra-based and a non-Epetra-based)"
				<< "\n***\n";

		//
		// Create two different vector spaces (one Epetra and one non-Epetra)
		// that should be compatible
		//
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
			if(verbose)
				out << "\nCreating TSFCore::EpetraVectorSpace using Epetra_MpiComm ...\n";
			epetra_comm = rcp(new Epetra_MpiComm(mpiComm));
			epetra_map = rcp(new Epetra_Map(-1,local_dim,0,*epetra_comm));
			epetra_vs = rcp(new EpetraVectorSpace(epetra_map));
#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
			// Non-Epetra vector space
			if(verbose)
				out << "\nCreating TSFCore::SimpleMPIVectorSpace ...\n";
			non_epetra_vs = rcp(new SimpleMPIVectorSpace<Scalar>(mpiComm,local_dim));
#endif
		}
		else {
#endif
			//
			// Create serial vector spaces (i.e. VectorSpace::isInCore()==true)
			//
			// Epetra vector space
			if(verbose)
				out << "\nCreating TSFCore::EpetraVectorSpace using Epetra_SerialComm ...\n";
			epetra_comm = rcp(new Epetra_SerialComm);
			epetra_map = rcp(new Epetra_LocalMap(local_dim,0,*epetra_comm));
			epetra_vs = rcp(new EpetraVectorSpace(epetra_map));
#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
			if(verbose)
				out << "\nCreating TSFCore::SerialVectorSpace ...\n";
			non_epetra_vs = rcp(new SerialVectorSpace<Scalar>(local_dim));
#endif
#ifdef RTOp_USE_MPI
		}
#endif

		const int global_dim = local_dim * numProc;

		if(verbose)
			out
				<< "\nscalar              = " << scalar
				<< "\nlocal_dim           = " << local_dim
				<< "\nglobal_dim          = " << global_dim
				<< "\nnum_mv_cols         = " << num_mv_cols
				<< "\nepetra_vs.dim()     = " << epetra_vs->dim()
				<< "\nnon_epetra_vs.dim() = " << non_epetra_vs->dim()
				<< std::endl;

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

		if(verbose)
			out
				<< "\n***"
				<< "\n*** (B) Testing Epetra and non-Epetra TSFCore wrapped objects"
				<< "\n***\n";

		//
		// Check for compatibility of the vector and Multi-vectors
		// w.r.t. RTOps
		//

		if(verbose)
			out
				<< "\n*** (B.1) Testing individual vector/multi-vector RTOps\n";

		assign( &*ev1, 0.0 );
		assign( &*ev2, scalar );
#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
		assign( &*nev1, 0.0 );
		assign( &*nev2, scalar );
#endif
		assign( &*eV1, 0.0 );
		assign( &*eV2, scalar );
#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
		assign( &*neV1, 0.0 );
		assign( &*neV2, scalar );
#endif

		Scalar
			ev1_nrm = norm_1(*ev1),
			ev2_nrm = norm_1(*ev2),
			eV1_nrm = norm_1(*eV1),
			eV2_nrm = norm_1(*eV2);
#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
		Scalar
			nev1_nrm = norm_1(*nev1),
			nev2_nrm = norm_1(*nev2),
			neV1_nrm = norm_1(*neV1),
			neV2_nrm = norm_1(*neV2);
#endif

		const std::string s1_n = "fabs(scalar)*global_dim";
		const double s1 = fabs(scalar)*global_dim;
		
		testRelErr("norm_1(ev1)",ev1_nrm,"0",0,"max_rel_err",max_rel_err,verbose,out) || (success=false);
		testRelErr("norm_1(ev2)",ev2_nrm,s1_n,s1,"max_rel_err",max_rel_err,verbose,out) || (success=false);
#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
		testRelErr("norm_1(nev1)",nev1_nrm,"0",0,"max_rel_err",max_rel_err,verbose,out) || (success=false);
		testRelErr("norm_1(nev2)",nev2_nrm,s1_n,s1,"max_rel_err",max_rel_err,verbose,out) || (success=false);
#endif
		testRelErr("norm_1(eV1)",eV1_nrm,"0",0,"max_rel_err",max_rel_err,verbose,out) || (success=false);
		testRelErr("norm_1(eV2)",eV2_nrm,s1_n,s1,"max_rel_err",max_rel_err,verbose,out) || (success=false);
#ifndef EPETRA_ADAPTERS_EPETRA_ONLY
		testRelErr("norm_1(neV1)",neV1_nrm,"0",0,"max_rel_err",max_rel_err,verbose,out) || (success=false);
		testRelErr("norm_1(neV2)",neV2_nrm,s1_n,s1,"max_rel_err",max_rel_err,verbose,out) || (success=false);
#endif

		if(verbose)
			out
				<< "\n*** (B.2) Test RTOps with two or more arguments\n";

		if(verbose) out << "\nPerforming ev1 = ev2 ...\n";
		timer.start(true);
 		assign( &*ev1, *ev2 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(ev1)",norm_1(*ev1),"norm_1(ev2)",ev2_nrm,"max_rel_err",max_rel_err,verbose,out) || (success=false);

		if(verbose) out << "\nPerforming eV1 = eV2 ...\n";
		timer.start(true);
 		assign( &*eV1, *eV2 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(eV1)",norm_1(*eV1),"norm_1(eV2)",eV2_nrm,"max_rel_err",max_rel_err,verbose,out) || (success=false);

#ifndef EPETRA_ADAPTERS_EPETRA_ONLY

		if(verbose) out << "\nPerforming ev1 = nev2 ...\n";
		timer.start(true);
 		assign( &*ev1, *nev2 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(ev1)",norm_1(*ev1),"norm_1(nev2)",nev2_nrm,"max_rel_err",max_rel_err,verbose,out) || (success=false);

		if(verbose) out << "\nPerforming nev1 = ev2 ...\n";
		timer.start(true);
 		assign( &*nev1, *ev2 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(nev1)",norm_1(*nev1),"norm_1(ev2)",ev2_nrm,"max_rel_err",max_rel_err,verbose,out) || (success=false);

		if(verbose) out << "\nPerforming nev1 = nev2 ...\n";
		timer.start(true);
 		assign( &*nev1, *nev2 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(nev1)",norm_1(*nev1),"norm_1(nev2)",nev2_nrm,"max_rel_err",max_rel_err,verbose,out) || (success=false);

		if(verbose) out << "\nPerforming eV1 = neV2 ...\n";
		timer.start(true);
 		assign( &*eV1, *neV2 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(eV1)",norm_1(*eV1),"norm_1(neV2)",neV2_nrm,"max_rel_err",max_rel_err,verbose,out) || (success=false);

		if(verbose) out << "\nPerforming neV1 = eV2 ...\n";
		timer.start(true);
 		assign( &*neV1, *eV2 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(neV1)",norm_1(*neV1),"norm_1(eV2)",eV2_nrm,"max_rel_err",max_rel_err,verbose,out) || (success=false);

		if(verbose) out << "\nPerforming neV1 = neV2 ...\n";
		timer.start(true);
 		assign( &*neV1, *neV2 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(neV1)",norm_1(*neV1),"norm_1(neV2)",neV2_nrm,"max_rel_err",max_rel_err,verbose,out) || (success=false);

#endif

		const std::string s2_n = "scalar^2*global_dim*num_mv_cols";
		const double s2 = scalar*scalar*global_dim*num_mv_cols;

		RefCountPtr<MultiVector<Scalar> >
			T = eV1->domain()->createMembers(num_mv_cols);

		if(verbose)
			out
				<< "\n*** (B.3) Test MultiVector::apply(...)\n";

		if(verbose) out << "\nPerforming eV1'*eV2 ...\n";
		timer.start(true);
		eV1->apply( TRANS, *eV2, &*T );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(eV1'*eV2)",norm_1(*T),s2_n,s2,"max_rel_err",max_rel_err,verbose,out) || (success=false);
		if(verbose && dumpAll) out << "\neV1'*eV2 =\n" << *T;

#ifndef EPETRA_ADAPTERS_EPETRA_ONLY

		if(verbose) out << "\nPerforming neV1'*eV2 ...\n";
		timer.start(true);
		neV1->apply( TRANS, *eV2, &*T );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(neV1'*eV2)",norm_1(*T),s2_n,s2,"max_rel_err",max_rel_err,verbose,out) || (success=false);
		if(verbose && dumpAll) out << "\nneV1'*eV2 =\n" << *T;

		if(verbose) out << "\nPerforming eV1'*neV2 ...\n";
		timer.start(true);
		eV1->apply( TRANS, *neV2, &*T );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(eV1'*neV2)",norm_1(*T),s2_n,s2,"max_rel_err",max_rel_err,verbose,out) || (success=false);
		if(verbose && dumpAll) out << "\neV1'*neV2 =\n" << *T;

		if(verbose) out << "\nPerforming neV1'*neV2 ...\n";
		timer.start(true);
		neV1->apply( TRANS, *neV2, &*T );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(neV1'*neV2)",norm_1(*T),s2_n,s2,"max_rel_err",max_rel_err,verbose,out) || (success=false);
		if(verbose && dumpAll) out << "\nneV1'*neV2 =\n" << *T;

#endif

#ifndef EPETRA_ADAPTERS_EXCLUDE_EPETRA_OPERATOR

		if(verbose)
			out
				<< "\n*** (B.4) Creating a diagonal Epetra_Operator\n";

		RefCountPtr<Epetra_Operator>  epetra_op;

		if(1) {
			// Create a diagonal matrix with scalar on the diagonal
			RefCountPtr<Epetra_CrsMatrix>
				epetra_mat = rcp(new Epetra_CrsMatrix(::Copy,*epetra_map,1));
			double values[1] = { scalar };
			int indices[1];
			const int IB = epetra_map->IndexBase(), offset = procRank*local_dim;
			for( int k = 0; k < local_dim; ++k ) {
				indices[0] = offset + k + IB;  // global column
				epetra_mat->InsertGlobalValues(
					offset + k + IB     // GlobalRow
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

		if(verbose)
			out
				<< "\n*** (B.5) Mix and match vector and Multi-vectors with Epetra opeator\n";

		const std::string s3_n = "2*scalar^2*global_dim";
		const double s3 = 2*scalar*scalar*global_dim;
		
		if(verbose) out << "\nPerforming ey = 2*Op*ev1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *ev1, &*ey, 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(ey)",norm_1(*ey),s3_n,s3,"max_rel_err",max_rel_err,verbose,out) || (success=false);

		if(verbose) out << "\nPerforming eY = 2*Op*eV1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *eV1, &*eY, 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(eY)",norm_1(*eY),s3_n,s3,"max_rel_err",max_rel_err,verbose,out) || (success=false);

#ifndef EPETRA_ADAPTERS_EPETRA_ONLY

		if(verbose) out << "\nPerforming ney = 2*Op*ev1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *ev1, &*ney, 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(ney)",norm_1(*ney),s3_n,s3,"max_rel_err",max_rel_err,verbose,out) || (success=false);

		if(verbose) out << "\nPerforming neY = 2*Op*eV1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *eV1, &*neY, 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(neY)",norm_1(*neY),s3_n,s3,"max_rel_err",max_rel_err,verbose,out) || (success=false);

		if(verbose) out << "\nPerforming ey = 2*Op*nev1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *nev1, &*ey, 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(ey)",norm_1(*ey),s3_n,s3,"max_rel_err",max_rel_err,verbose,out) || (success=false);

		if(verbose) out << "\nPerforming eY = 2*Op*neV1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *neV1, &*eY, 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(eY)",norm_1(*eY),s3_n,s3,"max_rel_err",max_rel_err,verbose,out) || (success=false);

		if(verbose) out << "\nPerforming ney = 2*Op*nev1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *nev1, &*ney, 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(ney)",norm_1(*ney),s3_n,s3,"max_rel_err",max_rel_err,verbose,out) || (success=false);

		if(verbose) out << "\nPerforming neY = 2*Op*neV1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *neV1, &*neY, 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(neY)",norm_1(*neY),s3_n,s3,"max_rel_err",max_rel_err,verbose,out) || (success=false);

#endif

		if(verbose)
			out
				<< "\n*** (B.6) Testing Multi-vector views with Epetra operator\n";

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

		if(verbose) out << "\nPerforming eY_v1 = 2*Op*eV1_v1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *eV1_v1, &*eY->subView(col_rng), 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(eY_v1)",norm_1(*eY->subView(col_rng)),s3_n,s3,"max_rel_err",max_rel_err,verbose,out) || (success=false);

		if(verbose) out << "\nPerforming eY_v2 = 2*Op*eV1_v2 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *eV1_v2, &*eY->subView(numCols,cols), 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(eY_v2)",norm_1(*eY->subView(numCols,cols)),s3_n,s3,"max_rel_err",max_rel_err,verbose,out) || (success=false);

#ifndef EPETRA_ADAPTERS_EPETRA_ONLY

		if(verbose) out << "\nPerforming neY_v1 = 2*Op*eV1_v1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *eV1_v1, &*neY->subView(col_rng), 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(neY_v1)",norm_1(*neY->subView(col_rng)),s3_n,s3,"max_rel_err",max_rel_err,verbose,out) || (success=false);

		if(verbose) out << "\nPerforming eY_v1 = 2*Op*neV1_v1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *neV1_v1, &*eY->subView(col_rng), 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(eY_v1)",norm_1(*eY->subView(col_rng)),s3_n,s3,"max_rel_err",max_rel_err,verbose,out) || (success=false);

		if(verbose) out << "\nPerforming neY_v2 = 2*Op*eV1_v2 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *eV1_v2, &*neY->subView(numCols,cols), 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(neY_v2)",norm_1(*neY->subView(numCols,cols)),s3_n,s3,"max_rel_err",max_rel_err,verbose,out) || (success=false);

		if(verbose) out << "\nPerforming eY_v2 = 2*Op*neV1_v2 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *neV1_v2, &*eY->subView(numCols,cols), 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		testRelErr("norm_1(eY_v2)",norm_1(*eY->subView(numCols,cols)),s3_n,s3,"max_rel_err",max_rel_err,verbose,out) || (success=false);

#endif

#endif // EPETRA_ADAPTERS_EXCLUDE_EPETRA_OPERATOR

		if(verbose)
			out
				<< "\n***"
				<< "\n*** (C) Comparing the speed of TSFCore adapted Epetra objects verses raw Epetra objects"
				<< "\n***\n";

		//
		// Setup the number of timing loops to get good timings
		//
		// Here we try to shoot for timing about 1 second's worth of
		// computations and adjust the number of evaluation loops
		// accordingly.  Let X be the approximate number of flops per
		// loop (or per evaluation).  We then compute the number of
		// loops as:
		//
		// 1.0 sec |  num CPU flops |   1 loop  |
		// --------|----------------|-----------|
		//         |       sec      |   X flops |
		//
		// This just comes out to be:
		//
		//   num_time_loops_X =  max_flop_rate / (X flops per loop)
		//
		// In this computation we ignore extra overhead that will be
		// an issue when local_dim is small.
		//

		double raw_epetra_time, tsfcore_wrapped_time;
		
		if(verbose)
			out
				<< "\n*** (C.1) Comparing the speed of RTOp verses raw Epetra_Vector operations\n";

		const double flop_adjust_factor_1 = 3.0;
		const int num_time_loops_1 = int( max_flop_rate / ( flop_adjust_factor_1 * local_dim * num_mv_cols ) ) + 1;

		if(1) {
				
			// Get constant references to Epetra_MultiVector objects in eV1 and eV2
			const Epetra_MultiVector
				&eeV2
				= *dyn_cast<const EpetraMultiVector>(const_cast<const MultiVector<Scalar>&>(*eV2)).epetra_multi_vec();
			
			// Get the Epetra_MultiVector object inside of Y to be modified.  Note that the following is the recommended
			// way to do this since it gives the greatest flexibility in the implementation of TSFCore::EpetraMultiVector.
			RefCountPtr<Epetra_MultiVector>  eeV1;
			RefCountPtr<const EpetraVectorSpace> eV1_range, eV1_domain;
			dyn_cast<EpetraMultiVector>(*eV1).setUninitialized(&eeV1,&eV1_range,&eV1_domain);
			
			if(verbose)
				out << "\nPerforming eeV1 = eeV2 (using raw Epetra_MultiVector::operator=(...)) " << num_time_loops_1 << " times ...\n";
			timer.start(true);
			for(int k = 0; k < num_time_loops_1; ++k ) {
				*eeV1 = eeV2;
			}
			timer.stop();
			raw_epetra_time = timer.totalElapsedTime();
			if(verbose) out << "  total time = " << raw_epetra_time << " sec\n";
			
			dyn_cast<EpetraMultiVector>(*eV1).initialize(eeV1,eV1_range,eV1_domain);
				
		}
		
		if(verbose)
			out << "\nPerforming eV1 = eV2 (using TSFCore::MPIMultiVectorBase::applyOp(...)) " << num_time_loops_1 << " times ...\n";
		timer.start(true);
		for(int k = 0; k < num_time_loops_1; ++k ) {
			assign( &*eV1, *eV2 );
		}
		timer.stop();
		tsfcore_wrapped_time = timer.totalElapsedTime();
		if(verbose) out << "  total time = " << tsfcore_wrapped_time << " sec\n";
		
		print_performance_stats( num_time_loops_1, raw_epetra_time, tsfcore_wrapped_time, verbose, out );

		// RAB: 2004/01/05: Note, the above relative performance is
		// likely to be the worst of all of the others since RTOp
		// operators are applied seperately column by column but the
		// relative performance should go to about 1.0 when local_dim
		// is sufficiently large!  However, because
		// Epetra_MultiVector::Assign(...) is implemented with a bad
		// algorithm (as of this 2004/01/05) with lots of cache misses
		// for the wrong problem sizes, the column-by-column RTOp
		// implementation used with the TSFCore adapters is actually
		// much faster in some cases.  However, the extra overhead of
		// RTOp is much worse for very very small (order 10) sizes.

		if(verbose)
			out
				<< "\n*** (C.2) Comparing TSFCore::MPIMultiVectorBase::apply() verses raw Epetra_MultiVector::Multiply()\n";

		const double flop_adjust_factor_2 = 2.0;
		const int num_time_loops_2 = int( max_flop_rate / ( flop_adjust_factor_2* local_dim * num_mv_cols * num_mv_cols ) ) + 1;

		if(1) {
			
			// Get constant references to Epetra_MultiVector objects in eV1 and eV2
			const Epetra_MultiVector
				&eeV1
				= *dyn_cast<const EpetraMultiVector>(const_cast<const MultiVector<Scalar>&>(*eV1)).epetra_multi_vec();
			const Epetra_MultiVector
				&eeV2
				= *dyn_cast<const EpetraMultiVector>(const_cast<const MultiVector<Scalar>&>(*eV2)).epetra_multi_vec();
			
			// Get the Epetra_MultiVector object inside of T to be modified.  Note that the following is the recommended
			// way to do this since it gives the greatest flexibility in the implementation of TSFCore::EpetraMultiVector.
			RefCountPtr<Epetra_MultiVector>  eT;
			RefCountPtr<const EpetraVectorSpace> eT_range, eT_domain;
			dyn_cast<EpetraMultiVector>(*T).setUninitialized(&eT,&eT_range,&eT_domain);
			
			if(verbose)
				out << "\nPerforming eeV1'*eeV2 (using raw Epetra_MultiVector::Multiply(...)) "	<< num_time_loops_2 << " times ...\n";
			timer.start(true);
			for(int k = 0; k < num_time_loops_2; ++k ) {
				eT->Multiply( 'T', 'N', 1.0, eeV1, eeV2, 0.0 );
			}
			timer.stop();
			raw_epetra_time = timer.totalElapsedTime();
			if(verbose) out << "  total time = " << raw_epetra_time << " sec\n";
			
			dyn_cast<EpetraMultiVector>(*T).initialize(eT,eT_range,eT_domain);
			
		}
		
		if(verbose)
			out << "\nPerforming eV1'*eV2 (using TSFCore::MPIMultiVectorBase::apply(...)) "	<< num_time_loops_2 << " times ...\n";
		timer.start(true);
		for(int k = 0; k < num_time_loops_2; ++k ) {
			eV1->apply( TRANS, *eV2, &*T );
		}
		timer.stop();
		tsfcore_wrapped_time = timer.totalElapsedTime();
		if(verbose) out << "  total time = " << tsfcore_wrapped_time << " sec\n";
	
		print_performance_stats( num_time_loops_2, raw_epetra_time, tsfcore_wrapped_time, verbose, out );
		
		// RAB: 2004/01/05: Note, even though the TSFCore adapter does
		// not actually call Epetra_MultiVector::Multiply(...), the
		// implementation in TSFCore::MPIMultiVectorBase::apply(...)
		// performs almost exactly the same flops and calls dgemm(...)
		// as well.  Herefore, except for some small overhead, the raw
		// Epetra and the TSFCore wrapped computations should give
		// almost identical times in almost all cases.

		if(verbose)
			out
				<< "\n*** (C.3) Comparing TSFCore::EpetraLinearOp::apply() verses raw Epetra_Operator::apply()\n";

		const double flop_adjust_factor_3 = 10.0; // lots of indirect addressing
		const int num_time_loops_3 = int( max_flop_rate / ( flop_adjust_factor_3 * local_dim * num_mv_cols ) ) + 1;

		if(1) {
			
			// Get constant references to Epetra_MultiVector objects in eV1 and eV2
			const Epetra_MultiVector
				&eeV1
				= *dyn_cast<const EpetraMultiVector>(const_cast<const MultiVector<Scalar>&>(*eV1)).epetra_multi_vec();
			
			// Get the Epetra_MultiVector object inside of Y to be modified.  Note that the following is the recommended
			// way to do this since it gives the greatest flexibility in the implementation of TSFCore::EpetraMultiVector.
			RefCountPtr<Epetra_MultiVector>  eeY;
			RefCountPtr<const EpetraVectorSpace> eY_range, eY_domain;
			dyn_cast<EpetraMultiVector>(*eY).setUninitialized(&eeY,&eY_range,&eY_domain);
			
			if(verbose)
				out << "\nPerforming eeY = 2*eOp*eeV1 (using raw Epetra_Operator::apply(...)) " << num_time_loops_3 << " times ...\n";
			epetra_op->SetUseTranspose(false);
			timer.start(true);
			for(int k = 0; k < num_time_loops_3; ++k ) {
				epetra_op->Apply( eeV1, *eeY );
				eeY->Scale(2.0);
			}
			timer.stop();
			raw_epetra_time = timer.totalElapsedTime();
			if(verbose) out << "  total time = " << raw_epetra_time << " sec\n";
			
			dyn_cast<EpetraMultiVector>(*eY).initialize(eeY,eY_range,eY_domain);
			
		}
		
		if(verbose)
			out << "\nPerforming eY = 2*Op*eV1 (using TSFCore::EpetraLinearOp::apply(...)) " << num_time_loops_3 << " times ...\n";
		timer.start(true);
		for(int k = 0; k < num_time_loops_3; ++k ) {
			Op->apply( NOTRANS, *eV1, &*eY, 2.0 );
		}
		timer.stop();
		tsfcore_wrapped_time = timer.totalElapsedTime();
		if(verbose) out << "  total time = " << tsfcore_wrapped_time << " sec\n";
		
		print_performance_stats( num_time_loops_3, raw_epetra_time, tsfcore_wrapped_time, verbose, out );

		// RAB: 2004/01/05: Note, the above Epetra adapter is a true
		// adapter and simply calls Epetra_Operator::apply(...) so
		// except for some small overhead, the raw Epetra and the
		// TSFCore wrapped computations should give about exactly the
		// same runtime for almost all cases.

		if(verbose) {
			if(success)
				out << "\nCongratulations! All of the tests seem to have run sucessfully!\n";
			else
				out << "\nOh no! at least one of the tests did not check out!\n";
		}

	} // end try
	catch( const std::exception &excpt ) {
		if(verbose)
			std::cerr << "*** Caught a standard exception (procRank="<<procRank<<"): " << excpt.what() << std::endl;
		success = -1;
	}
	catch( ... ) {
		if(verbose)
			std::cerr << "*** Caught an unknown exception (procRank="<<procRank<<")!\n";
		success = -1;
	}

 	MPI_Finalize();

	return (success ? 0 : -1);

}

} // namespace TSFCore

int main( int argc, char* argv[] ) {
	return TSFCore::main_body(argc,argv);
}
