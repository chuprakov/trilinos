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

// //////////////////////////////////////////////////////////////
// TSFProductVectorSpace.hpp

#ifndef TSFPRODUCTVECTORSPACE_HPP
#define TSFPRODUCTVECTORSPACE_HPP

#include "Teuchos_Array.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFDescribableByTypeID.hpp"
#include "TSFVectorSpace.hpp"

namespace TSFExtended
{
  using namespace Teuchos;

/** Implementation of a product vector space.
 *
 *  This class allows the construction and use of a product vector
 *  space.  the space can be constructed by either specifying the
 *  entire set of spaces or by constructing it on the fly.  This
 *  latter ability is useful in constructing Block Operators and
 *  Product Vectors on the fly. It also checks to make sure that
 *  things are compatible as they are built.
 *
 */





  template<class Scalar>
  class ProductVectorSpace : public TSFCore::VectorSpace<Scalar>, 
                             public DescribableByTypeID
  {
  public:

    /** Constructor that specifies the complete collection of spaces
     *  making up the product.  If a space is not set, an error is
     *  thrown.
     *
     * @param vecSpaces: Teuchos Array of vector spaces
     */

    ProductVectorSpace(const Teuchos::Array<const VectorSpace<Scalar> > 
		       &vecSpaces);
//     ProductVectorSpace(const Teuchos::Array<const int>
// 		       &vecSpaces);
  




    /** Empty constructor to allow the space to be build on the fly */ 
    ProductVectorSpace()
    {
      vecspaces_.resize(0);
      isSet_.resize(0);
      isFinal_ = false;
    }
  
  
    /** Method to indicate that a product vector space being
     *  constructed on the fly is complete.  The spaces are checked.
     */
    void finalize();



    /** returns the number of blocks */
    int numBlocks() const
    {
      return numBlocks_;
    }


    /** Returns the dimension */
    int dim() const
    {
      return dim_;
    }


    /** Checks to see if the current space is compatible with other
     *
     * @param other: VectorSpace to be checked with this
     */
    bool isCompatible(const VectorSpace<Scalar> 
		      &other ) const;


    /** Returns the kth block vector space */
    VectorSpace<Scalar> getBlock(const int k) const
    {
      TEST_FOR_EXCEPTION( k < 0 || k > numBlocks_, std::out_of_range,
			  "The value of k = " << k << " is out of range"
			  << endl);
      return vecSpaces_[k];
    }



    /** Test equality between two spaces.
     *
     *@param other: Vector space to be compared with this
     */
    bool operator==(const VectorSpace<Scalar>& other) const;


    /** Test inequality between two spaces.
     *
     *@param other: Vector space to be compared with this
     */
    bool operator!=(const VectorSpace<Scalar>& other) const;



    /** Set Block Sets the ith block to subSp if not already set.  If
     *  already set, checks to see if subSp is equal to vecSpaces_[i]
     *
     * @param i int: index of block to be set
     * @param subSp VectorSpace: the block to be set
     */
    void setBlock(const int i, const VectorSpace<Scalar> &subSp);
  



    Vector<Scalar> createMember() const;



    bool isInCore() const
    {
      return isInCore_;
    }

  private:

    Teuchos::Array<VectorSpace<Scalar> > vecSpaces_;
    int numBlocks_;
    int dim_;
    bool isFinal_;
    bool isInCore_;

  };

} // namespace TSFExtended

#endif // TSF_PRODUCT_VECTOR_SPACE_HPP
