// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// ///////////////////////////////
// cxx_main.cpp

#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFCoreEpetraLinearOp.hpp"
#include "TSFCoreEpetraMultiVector.hpp"
#include "TSFCoreTestingTools.hpp"
#include "TSFCoreLinearOpTester.hpp"
#include "TSFCoreMPIVectorSpaceStd.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_LocalMap.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_oblackholestream.hpp"
#ifdef RTOp_USE_MPI
#  include "TSFCoreSimpleMPIVectorSpace.hpp"
#  include "Epetra_MpiComm.h"
#endif

//
// Some helper functions
//

namespace {

void print_performance_stats(
	const int        num_time_samples
	,const double    raw_epetra_time
	,const double    tsfcore_wrapped_time
	,bool            verbose
	,std::ostream    &out
	)
{
	if(verbose)
		out
			<< "\nAverage times (out of " << num_time_samples << " samples):\n"
			<< "  Raw Epetra              = " << (raw_epetra_time/num_time_samples) << std::endl
			<< "  TSFCore Wrapped Epetra  = " << (tsfcore_wrapped_time/num_time_samples) << std::endl
			<< "\nRelative performance of TSFCore wrapped verses raw Epetra:\n"
			<< "  ( raw epetra time / tsfcore wrapped time ) = ( " << raw_epetra_time << " / " << tsfcore_wrapped_time << " ) = "
			<< (raw_epetra_time/tsfcore_wrapped_time) << std::endl;
}

} // namespace

/* Testing program for TSFCore/Epetra adpaters.
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
int main( int argc, char* argv[] ) {

	typedef double Scalar;
	typedef Teuchos::ScalarTraits<Scalar> ST;

	using Teuchos::dyn_cast;
	using Teuchos::CommandLineProcessor;
	using Teuchos::RefCountPtr;
	using Teuchos::rcp;
	using Teuchos::rcp_static_cast;
	using Teuchos::rcp_const_cast;

	using TSFCore::testRelErr;
	using TSFCore::NOTRANS;
	using TSFCore::TRANS;
	
	bool verbose = true;
	bool dumpAll = false;
	bool success = true;
	bool result;

	int procRank = 0;

	MPI_Init(&argc,&argv);
	
	try {

		Teuchos::Time timer("");

		//
		// Read options from the commandline
		//

		int     local_dim            = 1000;
		int     num_mv_cols          = 4;
		double  max_rel_err          = 1e-13;
		double  max_rel_warn         = 1e-15;
		double  scalar               = 1.5;
		double  max_flop_rate        = 2.0e8;
		bool    use_mpi_vec_spc_std  = true;
#ifdef RTOp_USE_MPI
		bool    useMPI               = true;
#endif
		CommandLineProcessor  clp(false); // Don't throw exceptions
		clp.setOption( "verbose", "quiet", &verbose, "Determines if any output is printed or not." );
		clp.setOption( "dump-all", "no-dump", &dumpAll, "Determines if quantities are dumped or not." );
		clp.setOption( "local-dim", &local_dim, "Number of vector elements per process." );
		clp.setOption( "num-mv-cols", &num_mv_cols, "Number columns in each multi-vector (>=4)." );
		clp.setOption( "max-rel-err-tol", &max_rel_err, "Maximum relative error tolerance for tests." );
		clp.setOption( "max-rel-warn-tol", &max_rel_warn, "Maximum relative warning tolerance for tests." );
		clp.setOption( "scalar", &scalar, "A scalar used in all computations." );
		clp.setOption( "max-flop-rate", &max_flop_rate, "Approx flop rate used for loop timing." );
		clp.setOption( "use-mpi-vec-spc-std", "no-use-mpi-vec-spc-std", &use_mpi_vec_spc_std
									 ,"Use MPIVectorSpaceStd or SerialVectorSpaceStd (serial), SimpleMPIVectorSpace (parallel)." );
#ifdef RTOp_USE_MPI
		clp.setOption( "use-mpi", "no-use-mpi", &useMPI, "Actually use MPI or just run independent serial programs." );
#endif
		CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

		TEST_FOR_EXCEPTION(
			num_mv_cols < 4, std::logic_error
			,"Error, num-mv-cols must be >= 4!"
			);

		//
		// Get basic MPI info
		//

#ifdef RTOp_USE_MPI
		MPI_Comm mpiComm;
		int numProc;
		if(useMPI) {
			mpiComm = MPI_COMM_WORLD;
			MPI_Comm_size( mpiComm, &numProc );
			MPI_Comm_rank( mpiComm, &procRank );
		}
		else {
			mpiComm = MPI_COMM_NULL;
			numProc = 1;
			procRank = 0;
		}
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
		// Create two different vector spaces (one Epetra and one
		// non-Epetra) that should be compatible.
		//
		RefCountPtr<const Epetra_Comm> epetra_comm;
		RefCountPtr<const Epetra_Map> epetra_map;
		RefCountPtr<const TSFCore::VectorSpace<Scalar> > epetra_vs;
		RefCountPtr<const TSFCore::VectorSpace<Scalar> > non_epetra_vs;
#ifdef RTOp_USE_MPI
		if(useMPI) {
			//
			// Create parallel vector spaces with compatible maps
			//
			// Epetra vector space
			if(verbose) out << "\nCreating TSFCore::EpetraVectorSpace using Epetra_MpiComm ...\n";
			epetra_comm = rcp(new Epetra_MpiComm(mpiComm));
			epetra_map = rcp(new Epetra_Map(-1,local_dim,0,*epetra_comm));
			epetra_vs = rcp(new TSFCore::EpetraVectorSpace(epetra_map));
			// Non-Epetra vector space
			if(use_mpi_vec_spc_std) {
				if(verbose) out << "\nCreating TSFCore::MPIVectorSpaceStd ...\n";
				non_epetra_vs = rcp(new TSFCore::MPIVectorSpaceStd<Scalar>(mpiComm,local_dim,-1));
			}
			else {
				if(verbose) out << "\nCreating TSFCore::SimpleMPIVectorSpace ...\n";
				non_epetra_vs = rcp(new TSFCore::SimpleMPIVectorSpace<Scalar>(mpiComm,local_dim));
			}
		}
		else {
#endif
			//
			// Create serial vector spaces (i.e. VectorSpace::isInCore()==true)
			//
			// Epetra vector space
			if(verbose) out << "\nCreating TSFCore::EpetraVectorSpace using Epetra_SerialComm ...\n";
			epetra_comm = rcp(new Epetra_SerialComm);
			epetra_map = rcp(new Epetra_LocalMap(local_dim,0,*epetra_comm));
			epetra_vs = rcp(new TSFCore::EpetraVectorSpace(epetra_map));
			// Non-Epetra vector space
			if(use_mpi_vec_spc_std) {
				if(verbose) out << "\nCreating TSFCore::MPIVectorSpaceStd ...\n";
				non_epetra_vs = rcp(new TSFCore::MPIVectorSpaceStd<Scalar>(MPI_COMM_NULL,local_dim,-1));
			}
			else {
				if(verbose) out << "\nCreating TSFCore::SerialVectorSpaceStd ...\n";
				non_epetra_vs = rcp(new TSFCore::SerialVectorSpaceStd<Scalar>(local_dim));
			}
#ifdef RTOp_USE_MPI
		}
#endif // end create vector spacdes [Doxygen looks for this!]

#ifdef RTOp_USE_MPI
		const int global_dim = local_dim * numProc;
#else
		const int global_dim = local_dim;
#endif

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

		RefCountPtr<TSFCore::Vector<Scalar> >
			ev1 = epetra_vs->createMember(),
			ev2 = epetra_vs->createMember();
		RefCountPtr<TSFCore::Vector<Scalar> >
			nev1 = non_epetra_vs->createMember(),
			nev2 = non_epetra_vs->createMember();

		RefCountPtr<TSFCore::MultiVector<Scalar> >
			eV1 = epetra_vs->createMembers(num_mv_cols),
			eV2 = epetra_vs->createMembers(num_mv_cols);
		RefCountPtr<TSFCore::MultiVector<Scalar> >
			neV1 = non_epetra_vs->createMembers(num_mv_cols),
			neV2 = non_epetra_vs->createMembers(num_mv_cols);

		if(verbose)
			out
				<< "\n***"
				<< "\n*** (B) Testing Epetra and non-Epetra TSFCore wrapped objects"
				<< "\n***\n";

		//
		// Check for compatibility of the vector and Multi-vectors
		// w.r.t. RTOps
		//

		if(verbose) out << "\n*** (B.1) Testing individual vector/multi-vector RTOps\n";

		TSFCore::assign( &*ev1, 0.0 );
		TSFCore::assign( &*ev2, scalar );
		TSFCore::assign( &*nev1, 0.0 );
		TSFCore::assign( &*nev2, scalar );
		TSFCore::assign( &*eV1, 0.0 );
		TSFCore::assign( &*eV2, scalar );
		TSFCore::assign( &*neV1, 0.0 );
		TSFCore::assign( &*neV2, scalar );

		Scalar
			ev1_nrm = TSFCore::norm_1(*ev1),
			ev2_nrm = TSFCore::norm_1(*ev2),
			eV1_nrm = TSFCore::norm_1(*eV1),
			eV2_nrm = TSFCore::norm_1(*eV2),
			nev1_nrm = TSFCore::norm_1(*nev1),
			nev2_nrm = TSFCore::norm_1(*nev2),
			neV1_nrm = TSFCore::norm_1(*neV1),
			neV2_nrm = TSFCore::norm_1(*neV2);

		const std::string s1_n = "fabs(scalar)*global_dim";
		const Scalar s1 = fabs(scalar)*global_dim;
		
		if(!testRelErr("TSFCore::norm_1(ev1)",ev1_nrm,"0",Scalar(0),"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
		if(verbose && dumpAll) out << "\nev1 =\n" << *ev1;
		if(!testRelErr("TSFCore::norm_1(ev2)",ev2_nrm,s1_n,s1,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
		if(verbose && dumpAll) out << "\nev2 =\n" << *ev2;
		if(!testRelErr("TSFCore::norm_1(nev1)",nev1_nrm,"0",Scalar(0),"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
		if(verbose && dumpAll) out << "\nnev2 =\n" << *ev1;
		if(!testRelErr("TSFCore::norm_1(nev2)",nev2_nrm,s1_n,s1,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
		if(verbose && dumpAll) out << "\nnev2 =\n" << *nev2;
		if(!testRelErr("TSFCore::norm_1(eV1)",eV1_nrm,"0",Scalar(0),"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
		if(verbose && dumpAll) out << "\neV1 =\n" << *eV1;
		if(!testRelErr("TSFCore::norm_1(eV2)",eV2_nrm,s1_n,s1,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
		if(verbose && dumpAll) out << "\neV2 =\n" << *eV2;
		if(!testRelErr("TSFCore::norm_1(neV1)",neV1_nrm,"0",Scalar(0),"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
		if(verbose && dumpAll) out << "\nneV1 =\n" << *neV1;
		if(!testRelErr("TSFCore::norm_1(neV2)",neV2_nrm,s1_n,s1,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
		if(verbose && dumpAll) out << "\nneV2 =\n" << *neV2;

		if(verbose) out << "\n*** (B.2) Test RTOps with two or more arguments\n";

		if(verbose) out << "\nPerforming ev1 = ev2 ...\n";
		timer.start(true);
 		TSFCore::assign( &*ev1, *ev2 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(ev1)",TSFCore::norm_1(*ev1),"TSFCore::norm_1(ev2)",ev2_nrm,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
		if(verbose && dumpAll) out << "\nev1 =\n" << *ev1;

		if(verbose) out << "\nPerforming eV1 = eV2 ...\n";
		timer.start(true);
 		TSFCore::assign( &*eV1, *eV2 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(eV1)",TSFCore::norm_1(*eV1),"TSFCore::norm_1(eV2)",eV2_nrm,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
		if(verbose && dumpAll) out << "\neV1 =\n" << *eV1;

		if(verbose) out << "\nPerforming ev1 = nev2 ...\n";
		timer.start(true);
 		TSFCore::assign( &*ev1, *nev2 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(ev1)",TSFCore::norm_1(*ev1),"TSFCore::norm_1(nev2)",nev2_nrm,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
		if(verbose && dumpAll) out << "\nev1 =\n" << *ev1;

		if(verbose) out << "\nPerforming nev1 = ev2 ...\n";
		timer.start(true);
 		TSFCore::assign( &*nev1, *ev2 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(nev1)",TSFCore::norm_1(*nev1),"TSFCore::norm_1(ev2)",ev2_nrm,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
		if(verbose && dumpAll) out << "\nnev1 =\n" << *nev1;

		if(verbose) out << "\nPerforming nev1 = nev2 ...\n";
		timer.start(true);
 		TSFCore::assign( &*nev1, *nev2 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(nev1)",TSFCore::norm_1(*nev1),"TSFCore::norm_1(nev2)",nev2_nrm,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
		if(verbose && dumpAll) out << "\nnev1 =\n" << *nev1;

		if(verbose) out << "\nPerforming eV1 = neV2 ...\n";
		timer.start(true);
 		TSFCore::assign( &*eV1, *neV2 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(eV1)",TSFCore::norm_1(*eV1),"TSFCore::norm_1(neV2)",neV2_nrm,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
		if(verbose && dumpAll) out << "\neV1 =\n" << *eV1;

		if(verbose) out << "\nPerforming neV1 = eV2 ...\n";
		timer.start(true);
 		TSFCore::assign( &*neV1, *eV2 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(neV1)",TSFCore::norm_1(*neV1),"TSFCore::norm_1(eV2)",eV2_nrm,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
		if(verbose && dumpAll) out << "\nneV1 =\n" << *neV1;

		if(verbose) out << "\nPerforming neV1 = neV2 ...\n";
		timer.start(true);
 		TSFCore::assign( &*neV1, *neV2 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(neV1)",TSFCore::norm_1(*neV1),"TSFCore::norm_1(neV2)",neV2_nrm,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
		if(verbose && dumpAll) out << "\nneV1 =\n" << *neV1;

		TSFCore::LinearOpTester<Scalar> linearOpTester(max_rel_warn,max_rel_err);

		if(verbose) out << "\n*** (B.3) Test Vector linear operator interface\n";

		if(verbose) out << "\nChecking out linear operator interface of ev1 ...\n";
		result = linearOpTester.check(*ev1,verbose?&out:NULL);
		if(!result) success = false;

		if(verbose) out << "\nChecking out linear operator interface of nev1 ...\n";
		result = linearOpTester.check(*nev1,verbose?&out:NULL);
		if(!result) success = false;

		if(verbose) out << "\n*** (B.4) Test MultiVector linear operator interface\n";

		if(verbose) out << "\nChecking out linear operator interface of eV1 ...\n";
		result = linearOpTester.check(*eV1,verbose?&out:NULL);
		if(!result) success = false;

		if(verbose) out << "\nChecking out linear operator interface of neV1 ...\n";
		result = linearOpTester.check(*neV1,verbose?&out:NULL);
		if(!result) success = false;

		const std::string s2_n = "scalar^2*global_dim*num_mv_cols";
		const Scalar s2 = scalar*scalar*global_dim*num_mv_cols;

		RefCountPtr<TSFCore::MultiVector<Scalar> >
			T = eV1->domain()->createMembers(num_mv_cols);

		if(verbose) out << "\n*** (B.5) Test MultiVector::apply(...)\n";

		if(verbose) out << "\nPerforming eV1'*eV2 ...\n";
		timer.start(true);
		eV1->apply( TRANS, *eV2, &*T );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(eV1'*eV2)",TSFCore::norm_1(*T),s2_n,s2,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
		if(verbose && dumpAll) out << "\neV1'*eV2 =\n" << *T;

		if(verbose) out << "\nPerforming neV1'*eV2 ...\n";
		timer.start(true);
		neV1->apply( TRANS, *eV2, &*T );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(neV1'*eV2)",TSFCore::norm_1(*T),s2_n,s2,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
		if(verbose && dumpAll) out << "\nneV1'*eV2 =\n" << *T;

		if(verbose) out << "\nPerforming eV1'*neV2 ...\n";
		timer.start(true);
		eV1->apply( TRANS, *neV2, &*T );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(eV1'*neV2)",TSFCore::norm_1(*T),s2_n,s2,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
		if(verbose && dumpAll) out << "\neV1'*neV2 =\n" << *T;

		if(verbose) out << "\nPerforming neV1'*neV2 ...\n";
		timer.start(true);
		neV1->apply( TRANS, *neV2, &*T );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(neV1'*neV2)",TSFCore::norm_1(*T),s2_n,s2,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
		if(verbose && dumpAll) out << "\nneV1'*neV2 =\n" << *T;

		if(verbose) out << "\n*** (B.6) Creating a diagonal Epetra_Operator Op\n";

		RefCountPtr<Epetra_Operator>  epetra_op;

		if(1) {
			// Create a diagonal matrix with scalar on the diagonal
			RefCountPtr<Epetra_CrsMatrix>
				epetra_mat = rcp(new Epetra_CrsMatrix(::Copy,*epetra_map,1));
			Scalar values[1] = { scalar };
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
		} // end epetra_op

		RefCountPtr<const TSFCore::LinearOp<Scalar> >
			Op = rcp(new TSFCore::EpetraLinearOp(epetra_op));

		if(verbose) out << "\n*** (B.7) Test EpetraLinearOp linear operator interface\n";

		if(verbose) out << "\nChecking out linear operator interface of Op ...\n";
		result = linearOpTester.check(*Op,verbose?&out:NULL);
		if(!result) success = false;

		RefCountPtr<TSFCore::Vector<Scalar> >
			ey  = epetra_vs->createMember();
		RefCountPtr<TSFCore::MultiVector<Scalar> >
			eY  = epetra_vs->createMembers(num_mv_cols);
		RefCountPtr<TSFCore::Vector<Scalar> >
			ney = non_epetra_vs->createMember();
		RefCountPtr<TSFCore::MultiVector<Scalar> >
			neY = non_epetra_vs->createMembers(num_mv_cols);

		if(verbose) out << "\n*** (B.8) Mix and match vector and Multi-vectors with Epetra opeator\n";

		const std::string s3_n = "2*scalar^2*global_dim";
		const Scalar s3 = 2*scalar*scalar*global_dim;
		
		if(verbose) out << "\nPerforming ey = 2*Op*ev1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *ev1, &*ey, 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(ey)",TSFCore::norm_1(*ey),s3_n,s3,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;

		if(verbose) out << "\nPerforming eY = 2*Op*eV1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *eV1, &*eY, 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(eY)",TSFCore::norm_1(*eY),s3_n,s3,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;

		if(verbose) out << "\nPerforming ney = 2*Op*ev1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *ev1, &*ney, 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(ney)",TSFCore::norm_1(*ney),s3_n,s3,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;

		if(verbose) out << "\nPerforming neY = 2*Op*eV1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *eV1, &*neY, 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(neY)",TSFCore::norm_1(*neY),s3_n,s3,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;

		if(verbose) out << "\nPerforming ey = 2*Op*nev1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *nev1, &*ey, 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(ey)",TSFCore::norm_1(*ey),s3_n,s3,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;

		if(verbose) out << "\nPerforming eY = 2*Op*neV1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *neV1, &*eY, 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(eY)",TSFCore::norm_1(*eY),s3_n,s3,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;

		if(verbose) out << "\nPerforming ney = 2*Op*nev1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *nev1, &*ney, 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(ney)",TSFCore::norm_1(*ney),s3_n,s3,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;

		if(verbose) out << "\nPerforming ney = 2*Op*nev1 through MultiVector interface ...\n";
		timer.start(true);
		Op->apply( NOTRANS, static_cast<const TSFCore::MultiVector<Scalar>&>(*nev1), static_cast<TSFCore::MultiVector<Scalar>*>(&*ney), 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(ney)",TSFCore::norm_1(*ney),s3_n,s3,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;

		if(verbose) out << "\nPerforming neY = 2*Op*neV1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *neV1, &*neY, 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(neY)",TSFCore::norm_1(*neY),s3_n,s3,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;

		if(verbose) out << "\n*** (B.9) Testing Multi-vector views with Epetra operator\n";

		const TSFCore::Range1D col_rng(1,2);
		const int numCols = 2;
		const int cols[] = { 3, 4 };

		RefCountPtr<const TSFCore::MultiVector<Scalar> >
			eV1_v1  = rcp_static_cast<const TSFCore::MultiVector<Scalar> >(eV1)->subView(col_rng),
			eV1_v2  = rcp_static_cast<const TSFCore::MultiVector<Scalar> >(eV1)->subView(numCols,cols);
		RefCountPtr<const TSFCore::MultiVector<Scalar> >
			neV1_v1  = rcp_static_cast<const TSFCore::MultiVector<Scalar> >(neV1)->subView(col_rng),
			neV1_v2  = rcp_static_cast<const TSFCore::MultiVector<Scalar> >(neV1)->subView(numCols,cols);

		if(verbose) out << "\nPerforming eY_v1 = 2*Op*eV1_v1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *eV1_v1, &*eY->subView(col_rng), 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(eY_v1)",TSFCore::norm_1(*eY->subView(col_rng)),s3_n,s3,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;

		if(verbose) out << "\nPerforming eY_v2 = 2*Op*eV1_v2 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *eV1_v2, &*eY->subView(numCols,cols), 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(eY_v2)",TSFCore::norm_1(*eY->subView(numCols,cols)),s3_n,s3,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;

		if(verbose) out << "\nPerforming neY_v1 = 2*Op*eV1_v1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *eV1_v1, &*neY->subView(col_rng), 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(neY_v1)",TSFCore::norm_1(*neY->subView(col_rng)),s3_n,s3,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;

		if(verbose) out << "\nPerforming eY_v1 = 2*Op*neV1_v1 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *neV1_v1, &*eY->subView(col_rng), 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(eY_v1)",TSFCore::norm_1(*eY->subView(col_rng)),s3_n,s3,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;

		if(verbose) out << "\nPerforming neY_v2 = 2*Op*eV1_v2 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *eV1_v2, &*neY->subView(numCols,cols), 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(neY_v2)",TSFCore::norm_1(*neY->subView(numCols,cols)),s3_n,s3,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;

		if(verbose) out << "\nPerforming eY_v2 = 2*Op*neV1_v2 ...\n";
		timer.start(true);
		Op->apply( NOTRANS, *neV1_v2, &*eY->subView(numCols,cols), 2.0 );
		timer.stop();
		if(verbose) out << "  time = " << timer.totalElapsedTime() << " sec\n";
		if(!testRelErr("TSFCore::norm_1(eY_v2)",TSFCore::norm_1(*eY->subView(numCols,cols)),s3_n,s3,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;

		if(verbose) out << "\n*** (B.10) Testing Vector and MultiVector view creation functions\n";

    if(1) {

      const std::string s_n = "fabs(scalar)*num_mv_cols";
      const Scalar s = fabs(scalar)*num_mv_cols;

      std::vector<Scalar>  t_raw_values( num_mv_cols );
      RTOpPack::MutableSubVectorT<Scalar> t_raw( 0, num_mv_cols, &t_raw_values[0], 1 );

      std::fill_n( t_raw_values.begin(), t_raw_values.size(), ST::zero() );
			TSFCore::assign( &*T->range()->createMemberView(t_raw), scalar );
      Teuchos::RefCountPtr<const TSFCore::Vector<Scalar> > t_view = T->range()->createMemberView(static_cast<RTOpPack::SubVectorT<Scalar>&>(t_raw));
      Scalar t_nrm = TSFCore::norm_1(*t_view);
      if(!testRelErr("TSFCore::norm_1(t_view)",t_nrm,s_n,s,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
      if(verbose && dumpAll) out << "\nt_view =\n" << *t_view;

#ifndef __sun // The sun compiler Forte Developer 5.4 does not destory temporaries properly and this does not work
      std::fill_n( t_raw_values.begin(), t_raw_values.size(), ST::zero() );
      TSFCore::assign( &*T->range()->TSFCore::VectorSpace<Scalar>::createMemberView(t_raw), scalar );
      t_view = T->range()->TSFCore::VectorSpace<Scalar>::createMemberView(static_cast<RTOpPack::SubVectorT<Scalar>&>(t_raw));
      t_nrm = TSFCore::norm_1(*t_view);
      if(!testRelErr("TSFCore::norm_1(t_view)",t_nrm,s_n,s,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
      if(verbose && dumpAll) out << "\nt_view =\n" << *t_view;
#endif

      std::vector<Scalar>  T_raw_values( num_mv_cols * num_mv_cols );
      RTOpPack::MutableSubMultiVectorT<Scalar> T_raw( 0, num_mv_cols, 0, num_mv_cols, &T_raw_values[0], num_mv_cols );

      std::fill_n( T_raw_values.begin(), T_raw_values.size(), ST::zero() );
      TSFCore::assign( &*T->range()->createMembersView(T_raw), scalar );
      Teuchos::RefCountPtr<const TSFCore::MultiVector<Scalar> >
				T_view = T->range()->createMembersView(static_cast<RTOpPack::SubMultiVectorT<Scalar>&>(T_raw));
      Scalar T_nrm = TSFCore::norm_1(*T_view);
      if(!testRelErr("TSFCore::norm_1(T_view)",T_nrm,s_n,s,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
      if(verbose && dumpAll) out << "\nT_view =\n" << *T_view;

#ifndef __sun // The sun compiler Forte Developer 5.4 does not destory temporaries properly and this does not work
      std::fill_n( T_raw_values.begin(), T_raw_values.size(), ST::zero() );
      TSFCore::assign( &*T->range()->TSFCore::VectorSpace<Scalar>::createMembersView(T_raw), scalar );
      T_view = T->range()->TSFCore::VectorSpace<Scalar>::createMembersView(static_cast<RTOpPack::SubMultiVectorT<Scalar>&>(T_raw));
      T_nrm = TSFCore::norm_1(*T_view);
      if(!testRelErr("TSFCore::norm_1(T_view)",T_nrm,s_n,s,"max_rel_err",max_rel_err,verbose?&out:NULL)) success=false;
      if(verbose && dumpAll) out << "\nT_view =\n" << *T_view;
#endif

    }

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
		
		if(verbose) out << "\n*** (C.1) Comparing the speed of RTOp verses raw Epetra_Vector operations\n";

		const double flop_adjust_factor_1 = 3.0;
		const int num_time_loops_1 = int( max_flop_rate / ( flop_adjust_factor_1 * local_dim * num_mv_cols ) ) + 1;

		if(1) {
				
			// Get constant references to Epetra_MultiVector objects in eV1 and eV2
			const Epetra_MultiVector
				&eeV2
				= *dyn_cast<const TSFCore::EpetraMultiVector>(const_cast<const TSFCore::MultiVector<Scalar>&>(*eV2)).epetra_multi_vec();
			
			// Get the Epetra_MultiVector object inside of Y to be modified.  Note that the following is the recommended
			// way to do this since it gives the greatest flexibility in the implementation of TSFCore::EpetraMultiVector.
			RefCountPtr<Epetra_MultiVector>  eeV1;
			RefCountPtr<const TSFCore::EpetraVectorSpace> eV1_range;
#ifdef TSFCORE_EPETRA_USE_EPETRA_DOMAIN_VECTOR_SPACE
      RefCountPtr<const TSFCore::EpetraVectorSpace> eV1_domain;
#else
      RefCountPtr<const TSFCore::ScalarProdVectorSpaceBase<Scalar> > eV1_domain;
#endif
			dyn_cast<TSFCore::EpetraMultiVector>(*eV1).setUninitialized(&eeV1,&eV1_range,&eV1_domain);
			
			if(verbose)
				out << "\nPerforming eeV1 = eeV2 (using raw Epetra_MultiVector::operator=(...)) " << num_time_loops_1 << " times ...\n";
			timer.start(true);
			for(int k = 0; k < num_time_loops_1; ++k ) {
				*eeV1 = eeV2;
			}
			timer.stop();
			raw_epetra_time = timer.totalElapsedTime();
			if(verbose) out << "  total time = " << raw_epetra_time << " sec\n";
			
			dyn_cast<TSFCore::EpetraMultiVector>(*eV1).initialize(eeV1,eV1_range,eV1_domain);
				
		}
		
		if(verbose)
			out << "\nPerforming eV1 = eV2 (using TSFCore::MPIMultiVectorBase::applyOp(...)) " << num_time_loops_1 << " times ...\n";
		timer.start(true);
		for(int k = 0; k < num_time_loops_1; ++k ) {
			TSFCore::assign( &*eV1, *eV2 );
		}
		timer.stop();
		tsfcore_wrapped_time = timer.totalElapsedTime();
		if(verbose) out << "  total time = " << tsfcore_wrapped_time << " sec\n";
		
		print_performance_stats( num_time_loops_1, raw_epetra_time, tsfcore_wrapped_time, verbose, out );

		// RAB: 2004/01/05: Note, the above relative performance is likely
		// to be the worst of all of the others since RTOp operators are
		// applied seperately column by column but the relative
		// performance should go to about 1.0 when local_dim is
		// sufficiently large!  However, because
		// Epetra_MultiVector::TSFCore::Assign(...) is implemented using double
		// pointer indexing, the RTOp implementation used with the TSFCore
		// adapters is actually faster in some cases.  However, the extra
		// overhead of RTOp is much worse for very very small (order 10)
		// sizes.

		if(verbose)
			out
				<< "\n*** (C.2) Comparing TSFCore::MPIMultiVectorBase::apply() verses raw Epetra_MultiVector::Multiply()\n";

		const double flop_adjust_factor_2 = 2.0;
		const int num_time_loops_2 = int( max_flop_rate / ( flop_adjust_factor_2* local_dim * num_mv_cols * num_mv_cols ) ) + 1;

		if(1) {
			
			// Get constant references to Epetra_MultiVector objects in eV1 and eV2
			const Epetra_MultiVector
				&eeV1
				= *dyn_cast<const TSFCore::EpetraMultiVector>(const_cast<const TSFCore::MultiVector<Scalar>&>(*eV1)).epetra_multi_vec();
			const Epetra_MultiVector
				&eeV2
				= *dyn_cast<const TSFCore::EpetraMultiVector>(const_cast<const TSFCore::MultiVector<Scalar>&>(*eV2)).epetra_multi_vec();
			
      Epetra_LocalMap eT_map(T->range()->dim(),0,*epetra_comm);
			Epetra_MultiVector eT(eT_map,T->domain()->dim());
			
			if(verbose)
				out << "\nPerforming eeV1'*eeV2 (using raw Epetra_MultiVector::Multiply(...)) "	<< num_time_loops_2 << " times ...\n";
			timer.start(true);
			for(int k = 0; k < num_time_loops_2; ++k ) {
				eT.Multiply( 'T', 'N', 1.0, eeV1, eeV2, 0.0 );
			}
			timer.stop();
			raw_epetra_time = timer.totalElapsedTime();
			if(verbose) out << "  total time = " << raw_epetra_time << " sec\n";
			
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

		if(verbose) out << "\n*** (C.3) Comparing TSFCore::EpetraLinearOp::apply() verses raw Epetra_Operator::apply()\n";

		const double flop_adjust_factor_3 = 10.0; // lots of indirect addressing
		const int num_time_loops_3 = int( max_flop_rate / ( flop_adjust_factor_3 * local_dim * num_mv_cols ) ) + 1;

		if(1) {
			
			// Get constant references to Epetra_MultiVector objects in eV1 and eV2
			const Epetra_MultiVector
				&eeV1
				= *dyn_cast<const TSFCore::EpetraMultiVector>(const_cast<const TSFCore::MultiVector<Scalar>&>(*eV1)).epetra_multi_vec();
			
			// Get the Epetra_MultiVector object inside of Y to be modified.  Note that the following is the recommended
			// way to do this since it gives the greatest flexibility in the implementation of TSFCore::EpetraMultiVector.
			RefCountPtr<Epetra_MultiVector>  eeY;
			RefCountPtr<const TSFCore::EpetraVectorSpace> eY_range;
#ifdef TSFCORE_EPETRA_USE_EPETRA_DOMAIN_VECTOR_SPACE
      RefCountPtr<const TSFCore::EpetraVectorSpace> eY_domain;
#else
      RefCountPtr<const TSFCore::ScalarProdVectorSpaceBase<Scalar> > eY_domain;
#endif
			dyn_cast<TSFCore::EpetraMultiVector>(*eY).setUninitialized(&eeY,&eY_range,&eY_domain);
			
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
			
			dyn_cast<TSFCore::EpetraMultiVector>(*eY).initialize(eeY,eY_range,eY_domain);
			
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
		// adapter and simply calls Epetra_Operator::apply(...) so, except
		// for some small overhead, the raw Epetra and the TSFCore wrapped
		// computations should give about exactly the same runtime for
		// almost all cases.

		if(verbose) {
			if(success) out << "\nCongratulations! All of the tests seem to have run sucessfully!\n";
			else        out << "\nOh no! at least one of the tests did not check out!\n";
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

} // end main()
