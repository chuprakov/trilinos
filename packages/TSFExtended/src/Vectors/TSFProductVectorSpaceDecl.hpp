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

#ifndef TSFPRODUCTVECTORSPACEDECL_HPP
#define TSFPRODUCTVECTORSPACEDECL_HPP

#include "Teuchos_Array.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "TSFCoreProductVectorSpaceDecl.hpp"
#include "TSFVecSpaceDescribableByTypeID.hpp"
#include "TSFVectorSpaceDecl.hpp"

namespace TSFExtended
{
  using Teuchos::Array;

  /** Implementation of a product vector space. */

 



  template<class Scalar>
  class ProductVectorSpace : public TSFCore::ProductVectorSpace<Scalar>, 
                             public DescribableByTypeID,
                             public Handleable<const TSFCore::VectorSpace<Scalar> >
  {
  public:
    GET_RCP(const TSFCore::VectorSpace<Scalar>);

    /** Constructor that specifies the complete collection of spaces
     *  making up the product.  
     *
     * @param vecSpaces: Teuchos Array of vector spaces
     */

    ProductVectorSpace(const Teuchos::Array<const VectorSpace<Scalar> > 
		       &vecSpaces);
    //     ProductVectorSpace(const Teuchos::Array<const int>
    // 		       &vecSpaces);
  


    /** Convenience method for block of two spaces  */
    ProductVectorSpace(const VectorSpace<Scalar>& s1,
		       const VectorSpace<Scalar>& s2);


  


    /** Returns the kth block vector space */
    VectorSpace<Scalar> getBlock(const int& k) const
    {
      TEST_FOR_EXCEPTION( k < 0 || k > numBlocks_, std::out_of_range,
			  "The value of k = " << k << " is out of range"
			  << endl);
      return vecSpaces_[k];
    }




    /** Finalize the space  */
    void finalize();

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





    RefCountPtr<TSFCore::Vector<Scalar> > createMember() const;



    string describe(int k) const;

  private:
    /** Empty constructor not to be used.  */
    ProductVectorSpace(){;}

    Teuchos::Array<VectorSpace<Scalar> > vecSpaces_;
    Teuchos::Array<int> isSet_;
    int numBlocks_;
    int dim_;
    bool isFinal_;
    bool isInCore_;

    /**  Private method to set up the underlying Core Space  */
    void setUpCore();
    

  };

} // namespace TSFExtended



						 
#endif // TSF_PRODUCT_VECTOR_SPACE_HPP
