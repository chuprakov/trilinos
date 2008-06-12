// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
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

#ifndef TPETRA_VECTORDATA_DECL_HPP
#define TPETRA_VECTORDATA_DECL_HPP

#include <Teuchos_BLAS.hpp>
#include <Tpetra_VectorSpace.hpp>

namespace Tpetra {

  template<typename OrdinalType, typename ScalarType>
  class VectorData : public Object {

  friend class Vector<OrdinalType, ScalarType>;

  public:
    VectorData(VectorSpace<OrdinalType, ScalarType> const& VectorSpace, 
           OrdinalType length, ScalarType seed); 

    ~VectorData();

  protected:
    Teuchos::BLAS<OrdinalType, ScalarType> BLAS_;
    VectorSpace<OrdinalType, ScalarType> VectorSpace_;
    std::vector<ScalarType> scalarArray_;
    ScalarType seed_;
  
  private:
    //! Copy constructor (declared but not defined, do not use)
    VectorData(VectorData<OrdinalType, ScalarType> const& Source);
    //! Assignment operator (declared but not defined, do not use)
    VectorData<OrdinalType, ScalarType>& operator = (VectorData<OrdinalType, ScalarType> const& Source);

  }; // class VectorData

} // namespace Tpetra

#endif // TPETRA_VECTORDATA_DECL_HPP
