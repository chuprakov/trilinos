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

#ifndef SCHURFACTORY_H
#define SCHURFACTORY_H

#include "SchurFactoryBase.h"
#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include "TSFLinearOperator.h"
#include "TSFOperatorSource.h"
#include "TSFOperatorSourceBase.h"
//#include "RightBlockNSOperatorSource.h"
//#include "KayLoghinRightOperatorSource.h"
#include "TSFArray.h"
#include <iostream>
#include <string>

namespace Meros
{

  using std::string;
  using std::ostream;
  using namespace TSF;

  /** \ingroup Preconditioner
   * SchurFactory builds an implementation-specific problem
   * from an abstract specification.
   *
   */

  class SchurFactory
    {
    public:
      /** empty ctor. Constructs a null vector */
      SchurFactory();

      /** assume control of a pointer */
      SchurFactory(SchurFactoryBase* ptr);

      /** get a concrete linear operator for Xinv */
      virtual TSFLinearOperator getSchurInvApprox(const TSFOperatorSource& opSrc) const;

      /** get a concrete linear operator for Xinv */
      //      TSFLinearOperator getSchurInvApprox(const KayLoghinRightOperatorSource& opSrc) const;

      /** write to a string */
      string toString() const ;

    private:
      TSFSmartPtr<SchurFactoryBase> ptr_;
    };

  /** \relates SchurFactory
   * write to an output stream
   */
  ostream& operator<<(ostream& os, const SchurFactory& x);

  /** \relates SchurFactory
   * write to a string
   */
  string toString(const SchurFactory& x);

}


#endif
