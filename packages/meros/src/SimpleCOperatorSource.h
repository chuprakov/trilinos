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

#ifndef SIMPLECOPERATORSOURCE_H
  #define SIMPLECOPERATORSOURCE_H

  #include "TSFDefs.h"
  #include "TSFSmartPtr.h"
  #include "TSFLinearOperator.h"
  #include "TSFOperatorSource.h"
  #include "RightBlockNSOperatorSource.h"
  #include "TSFArray.h"
  #include <iostream>
  #include <string>

  namespace SPP
  {
  using namespace TSF;
  using std::string;
  using std::ostream;

  // Need a better name for this class.

  class SimpleCOperatorSource : public RightBlockNSOperatorSource
    {
    public:
      /** ctor */
      /** ctor (no Mp) */
      SimpleCOperatorSource(TSFLinearOperator& S);

      /** virtual destructor */
      virtual ~SimpleCOperatorSource(){;}

      /** get a concrete linear operator */
      TSFLinearOperator getOp() const;

      /** get pressure Poisson operator */
      TSFLinearOperator getDinv() const;

      /** write to a string */
      string toString() const ;

    private:
      TSFLinearOperator S_;
      mutable TSFLinearOperator Dinv_;
      mutable bool hasDinv_;
    };
   }


 #endif
