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

// ///////////////////////////////////////////////////////////////////
// TSFCoreLinearOpTester.hpp

#ifndef TSFCORE_LINEAR_OP_TESTER_HPP
#define TSFCORE_LINEAR_OP_TESTER_HPP

#include "TSFCoreLinearOpTesterDecl.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "TSFCoreTestingTools.hpp"
#include "TSFCoreLinearOp.hpp"

namespace TSFCore {

template<class Scalar>
LinearOpTester<Scalar>::LinearOpTester(
	const ScalarMag      warning_tol
	,const ScalarMag     error_tol
	)
	:warning_tol_(warning_tol)
	,error_tol_(error_tol)
{}

template<class Scalar>
bool LinearOpTester<Scalar>::check(
	const LinearOp<Scalar>      &op
	,std::ostream               *out
	) const
{
	typedef Teuchos::ScalarTraits<Scalar> ST;
	bool success = true, result;
	const Scalar one = ST::one();
	const Scalar half = Scalar(0.5)*one;
	if(out) {
		*out
			<< "\n*** Entering LinearOpTester<Scalar>::check(...)\n"
			<< "\nTesting an operator of concrete type \'" << typeid(op).name() << "\' ...\n";
	}

	if(out)
		*out << "\nChecking the domain and range spaces ...\n";

	Teuchos::RefCountPtr<const VectorSpace<Scalar> >
		domain = op.domain(),
		range  = op.range();

	if(out)	*out << "op.domain().get() != NULL ? ";
	result = domain.get() != NULL;
	if(!result) success = false;
	if(out)	*out << (result ? "passed" : "failed" ) << std::endl;

	if(out)	*out << "op.range().get() != NULL ? ";
	result = range.get() != NULL;
	if(!result) success = false;
	if(out)	*out << (result ? "passed" : "failed" ) << std::endl;
	
	if(out)
		*out
			<< "\nChecking that the transpose agrees with the non-transpose operator as:"
			<< "\n  <0.5*op*v2,v1> == <v2,0.5*op*v1>"
			<< "\n   \\_______/            \\_______/"
			<< "\n      v4                   v3"
			<< "\n"
			<< "\n         <v4,v2> == <v2,v3>"
			<< std::endl;

	if(out) *out << "\nv1 = randomize(-1,+1); ...\n" ;
	Teuchos::RefCountPtr<Vector<Scalar> >	v1 = domain->createMember();
	TSFCore::randomize( -one, +one, &*v1 );

	if(out) *out << "\nv2 = randomize(-1,+1); ...\n" ;
	Teuchos::RefCountPtr<Vector<Scalar> >	v2 = range->createMember();
	TSFCore::randomize( -one, +one, &*v2 );

	if(out) *out << "\nv3 = 0.5*op*v1 ...\n" ;
	Teuchos::RefCountPtr<Vector<Scalar> >	v3 = range->createMember();
	op.apply( NOTRANS, *v1, &*v3, half );

	if(out) *out << "\nv4 = 0.5*op'*v2 ...\n" ;
	Teuchos::RefCountPtr<Vector<Scalar> >	v4 = domain->createMember();
	op.apply( TRANS, *v2, &*v4, half );

	const Scalar
		prod1 = domain->scalarProd(*v1,*v4),
		prod2 = range->scalarProd(*v2,*v3);
	
	if(out)
		*out
			<< "\nprod1 = <v1,v4> = " << prod1
			<< "\nprod2 = <v2,v3> = " << prod2 << std::endl;
		
	const Scalar
		rel_err = relErr(prod1,prod2);
	
	if(out)
		*out << "\nrel_err(prod1,prod2) = " << rel_err;

	result = rel_err <= error_tol();
	if(!result) success = false;

	if( !result ) {
		if(out)
			*out << " > error_tol() = " << error_tol() << ": Error!\n";
	}
	else {
		result = rel_err <= warning_tol();
		if(out) {
			if(!result)
				*out << " > warning_tol() = " << warning_tol() << ": Warning!\n";
			else
				*out << ": passed!\n";
		}
				
	}

	if(out)
		*out << "\n*** Leaving LinearOpTester<Scalar>::check(...)\n";

	return success;
}

} // namespace TSFCore

#endif // TSFCORE_LINEAR_OP_TESTER_HPP
