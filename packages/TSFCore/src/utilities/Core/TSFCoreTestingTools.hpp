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

// /////////////////////////////////////////////////////////////////////////
// TSFCoreTestingTools.hpp

#ifndef TSFCORE_TESTING_TOOLS_HPP
#define TSFCORE_TESTING_TOOLS_HPP

#include "TSFCoreTestingToolsDecl.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "TSFCoreLinearOpHandle.hpp"

template <class Scalar>
Scalar TSFCore::relErr( const Scalar &s1, const Scalar &s2 )
{
	return
		Teuchos::ScalarTraits<Scalar>::magnitude( s1 - s2 )
		/
		(
			Teuchos::ScalarTraits<Scalar>::eps()
			+ std::max( Teuchos::ScalarTraits<Scalar>::magnitude(s1), Teuchos::ScalarTraits<Scalar>::magnitude(s1) )
			);
}

template<class Scalar>
bool TSFCore::testRelErr(
	const std::string    &v1_name
	,const Scalar        &v1
	,const std::string   &v2_name
	,const Scalar        &v2
	,const std::string   &maxRelErr_name
	,const Scalar        &maxRelErr
	,std::ostream        *out
	)
{
	typedef Teuchos::ScalarTraits<Scalar> ST;
	const Scalar rel_err = relErr( v1, v2 );
	const bool success = ( ST::magnitude(rel_err) <= ST::magnitude(maxRelErr) );
	if(out) {
		*out
			<< "\nCheck: rel_err(" << v1_name << "," << v2_name << ")\n"
			<< "       = rel_err(" << v1 << "," << v2 << ") "
			<< "= " << rel_err
			<< " <= " << maxRelErr_name << " = " << maxRelErr << " : "
			<<  passfail(success) << std::endl;
	}
	return success;
}

template<class Scalar>
std::ostream& TSFCore::operator<<( std::ostream& o, const Vector<Scalar>& v )
{
	return v.describe(o,Teuchos::VERB_HIGH);
}

template<class Scalar>
std::ostream& TSFCore::operator<<( std::ostream& o, const LinearOp<Scalar>& M )
{
	return o << TSFCore::LinearOpHandle<Scalar>(Teuchos::rcp(&M,false));
}

template<class Scalar>
std::ostream& TSFCore::operator<<( std::ostream& o, const LinearOpHandle<Scalar>& M )
{
	return M.describe(o,Teuchos::VERB_EXTREME);
}

#endif // TSFCORE_TESTING_TOOLS_HPP
