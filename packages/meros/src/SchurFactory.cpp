// @HEADER
// ***********************************************************************
// 
//              Meros: Segregated Preconditioning Package
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

#include "SchurFactory.h"

#include "SchurFactoryBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFLinearOperatorBase.h"
#include "TSFOperatorSourceBase.h"


using namespace SPP;
using namespace TSF;

using std::string;

SchurFactory::SchurFactory()
	: ptr_(0)
{}


SchurFactory
::SchurFactory(SchurFactoryBase* ptr)
	: ptr_(ptr)
{}


TSFLinearOperator SchurFactory::getSchurInvApprox(const TSFOperatorSource& opSrc) const
{
	if (ptr_.get()==0) return TSFLinearOperator(); 
	return ptr_->getSchurInvApprox(opSrc);
}

// TSFLinearOperator SchurFactory::getSchurInvApprox(const KayLoghinRightOperatorSource& opSrc) const
// {
// 	if (ptr_.get()==0) return TSFLinearOperator(); 
// 	return ptr_->getSchurInvApprox(opSrc);
// }


string SchurFactory::toString() const
{
	return ptr_->toString();
}

namespace SPP
{
	ostream& operator<<(ostream& os, const SchurFactory& x)
	{
		return os << x.toString();
	}

	string toString(const SchurFactory& x)
	{
			return x.toString();
		}
}
