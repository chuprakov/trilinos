// @HEADER
// ***********************************************************************
// 
//      TSFCoreUtils: Trilinos Solver Framework Utilities Package 
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

#include "RTOpPack_RTOpT.hpp"
#include "RTOpPack_MPI_apply_op.hpp"
#include "RTOpPack_ROpCountNanInf.hpp"
#include "RTOpPack_ROpSum.hpp"
#include "RTOpPack_TOpAssignScalar.hpp"
//#include "RTOpPack_RTOpC.hpp"
//#include "RTOp_ROp_sum.h"
//#include "RTOp_TOp_assign_scalar.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_arrayArg.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_dyn_cast.hpp"

//
// Simple templated testing function
//

namespace {

template<class Scalar>
void test_do_stuff( MPI_Comm mpiComm, const int n, std::ostream &out )
{
  using Teuchos::arrayArg;
  using RTOpPack::SubVectorT;
  using RTOpPack::MutableSubVectorT;
  using RTOpPack::ReductTarget;
  out << "\n*** Entering test_do_stuff<"<<Teuchos::ScalarTraits<Scalar>::name()<<">) ...\n";
  std::vector<Scalar> z0_store(n);
  MutableSubVectorT<Scalar> z0(0,n,&z0_store[0],1);
  out << "\nPerforming: z0 = 1.0\n";
  RTOpPack::TOpAssignScalar<Scalar> assign_scalar_op(Scalar(1.0));
  RTOpPack::MPI_apply_op(
    mpiComm, assign_scalar_op, -1
    ,0,(const RTOpPack::SubVectorT<Scalar>*)NULL
    ,1,arrayArg(z0)()
    ,NULL
    );
  out << "\nPerforming sum(z0) = ";
  RTOpPack::ROpSum<Scalar> sum_op;
  Teuchos::RefCountPtr<ReductTarget> sum_targ = sum_op.reduct_obj_create();
  RTOpPack::MPI_apply_op(
    mpiComm, sum_op, -1
    ,1,arrayArg(z0)()
    ,0,(const RTOpPack::MutableSubVectorT<Scalar>*)NULL
    ,&*sum_targ
    );
  out << sum_op(*sum_targ) << std::endl;
  out << "\nisnaninf(sum(z0)) = " << Teuchos::ScalarTraits<Scalar>::isnaninf(sum_op(*sum_targ)) << std::endl;;
  out << "\nPerforming z0 = nan\n";
  assign_scalar_op.alpha(Teuchos::ScalarTraits<Scalar>::nan());
  RTOpPack::MPI_apply_op(
    mpiComm, assign_scalar_op, -1
    ,0,(const RTOpPack::SubVectorT<Scalar>*)NULL
    ,1,arrayArg(z0)()
    ,NULL
    );
/*
  out << "\ncountNanInf(z0) = ";
  RTOpPack::ROpCountNanInf<Scalar> countNanInf_op;
  Teuchos::RefCountPtr<ReductTarget> countNanInf_targ = countNanInf_op.reduct_obj_create();
  RTOpPack::MPI_apply_op(
    mpiComm, countNanInf_op, -1
    ,1,arrayArg(z0)()
    ,0,(const RTOpPack::MutableSubVectorT<Scalar>*)NULL
    ,&*sum_targ
    );
  out << countNanInf_op(*countNanInf_targ) << std::endl;
*/
  out << "\nPerforming sum(z0) = ";
  sum_op.reduct_obj_reinit(&*sum_targ);
  RTOpPack::MPI_apply_op(
    mpiComm, sum_op, -1
    ,1,arrayArg(z0)()
    ,0,(const RTOpPack::MutableSubVectorT<Scalar>*)NULL
    ,&*sum_targ
    );
  out << sum_op(*sum_targ) << std::endl;
  out << "\nisnaninf(sum(z0)) = " << Teuchos::ScalarTraits<Scalar>::isnaninf(sum_op(*sum_targ)) << std::endl;;
  out << "\n*** Leaving test_do_stuff<"<<Teuchos::ScalarTraits<Scalar>::name()<<">) ...\n";
}

} // namespace

int main( int argc, char* argv[] ) {

	using Teuchos::CommandLineProcessor;


	MPI_Init(&argc,&argv);

	bool success = true;
	bool verbose = true;
	int procRank = 0;

	try {

    //
    // Get basic MPI information
		MPI_Comm mpiComm = MPI_COMM_NULL;
#ifdef RTOp_USE_MPI
		mpiComm = MPI_COMM_WORLD;
		MPI_Comm_rank( mpiComm, &procRank );
#endif

		//
		// Setup the output stream
		//
		
		Teuchos::oblackholestream black_hole_out;
		std::ostream &out = ( procRank == 0 ? std::cout : black_hole_out );

		//
		// Read options from the commandline
		//

    int   n = 4;

		CommandLineProcessor  clp(false); // Don't throw exceptions
		clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
		clp.setOption( "n", &n, "Number of local vector elements." );
		CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    //
    // Use fully templated operators
    //

    test_do_stuff<float>(mpiComm,n,out);
    test_do_stuff<double>(mpiComm,n,out);
#if defined(HAVE_COMPLEX) && defined(HAVE_TEUCHOS_COMPLEX)
    test_do_stuff<std::complex<float> >(mpiComm,n,out);
    test_do_stuff<std::complex<double> >(mpiComm,n,out);
#endif

/*

    //
    // Use adapters for C RTOp_RTOp operator subclasses
    //
    
    if(1) {
      typedef RTOp_value_type Scalar;
      using Teuchos::arrayArg;
      using RTOpPack::SubVectorT;
      using RTOpPack::MutableSubVectorT;
      using RTOpPack::ReductTarget;
      out << "\n*** Testing with RTOpC wrapped RTOp_RTOp subclasses ...\n";
      std::vector<Scalar>  z0_store(n);
      MutableSubVectorT<Scalar> z0(0,n,&z0_store[0],1);
      out << "\nPerforming: z0 = 1.0\n";
      RTOpPack::RTOpC  assign_scalar_op;
      TEST_FOR_EXCEPTION(
        0!=RTOp_TOp_assign_scalar_construct(1.0,&assign_scalar_op.op())
        ,std::logic_error,"Error!"
        );
      RTOpPack::MPI_apply_op(
        mpiComm, assign_scalar_op, -1
        ,0,(const RTOpPack::SubVectorT<Scalar>*)NULL
        ,1,arrayArg(z0)()
        ,NULL
        );
      out << "\nPerforming sum(z0) = ";
      RTOpPack::RTOpC  sum_op;
      TEST_FOR_EXCEPTION(
        0!=RTOp_ROp_sum_construct(&sum_op.op())
        ,std::logic_error,"Error!"
        );
      Teuchos::RefCountPtr<ReductTarget> sum_targ = sum_op.reduct_obj_create();
      RTOpPack::MPI_apply_op(
        mpiComm, sum_op, -1
        ,1,arrayArg(z0)()
        ,0,(const RTOpPack::MutableSubVectorT<Scalar>*)NULL
        ,&*sum_targ
        );
      out << RTOp_ROp_sum_val(sum_op(*sum_targ)) << std::endl;
    }

*/

    // Do more stuff, Blah blah blah ...

		if(verbose)
			out << "\nAll of the tests seem to have run sucessfully!\n";

	} // end try
	catch( const std::exception &excpt ) {
		if(verbose)
			std::cerr << "*** Caught a standard exception : " << excpt.what() << std::endl;
		return -1;
	}
	catch( ... ) {
		if(verbose)
			std::cerr << "*** Caught an unknown exception!\n";
		success = -1;
	}

 	MPI_Finalize();

	return (success ? 0 : -1);

}
