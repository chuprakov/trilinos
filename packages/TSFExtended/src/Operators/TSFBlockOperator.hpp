/* @HEADER@ */
/* ***********************************************************************
//
//           TSFExtended: Trilinos Solver Framework Extended
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
// **********************************************************************/
/* @HEADER@ */

#ifndef TSFBLOCKOPERATOR_HPP
#define TSFBLOCKOPERATOR_HPP

#include "TSFConfigDefs.hpp"
#include "Teuchos_Array.hpp"
#include "TSFOpDescribableByTypeID.hpp"
 //#include "TSFVectorSpace.hpp"
#include "TSFLinearOperatorDecl.hpp"
#include "TSFProductVectorSpace.hpp"





namespace TSFExtended
{
  using namespace Teuchos;

  /**
   * Class BlockOperator provides an abstract interface for configuration
   * and building block operators
   *
   * @author Paul T Boggs (ptboggs@sandia.gov)
   */
  template <class Scalar>
  class BlockOperator : public OpDescribableByTypeID<Scalar>
  {
  public:

    /** ctor with domain and range specified.  The blocks must be
     *	specified later and all filled before use.
     */
    BlockOperator(const VectorSpace<Scalar> &domain, 
		  const VectorSpace<Scalar> &range);
    
    /** Empty ctor to allow the BlockOperator to be built on the fly.
     *	Note that the user must call finalize when the specification
     *	of the blocks is finished.  The method, finalize, allows
     *	unspecified blocks to be automatically filled with
     *	ZeroOperators
     */
    BlockOperator();

    /** Virtual dtor */
    virtual ~BlockOperator(){;}




    /** get the number of block rows */
    int numBlockRows() const {return nBlockRows_;}

    /** get the number of block columns */
    int numBlockCols() const {return nBlockCols_;}

    /** get the (i,j)-th submatrix */
    LinearOperator<Scalar> getBlock(const int &i, const int &j) const ;

    /** set the (i,j)-th submatrix 
	If the domain and/or the range are not set, then we
	are building the operator
    */
    void setBlock(const int &i, const int &j, 
		  const LinearOperator<Scalar>& sub);

    /** Finalize the matrix by setting the range and domain and 
     *  filling the non-set blocks with zero if specified.  It also
     *  checks to be sure all is correct.
     *
     *  @param zeroFill bool set to true if the non-set blocks are
     *  to be set to zero
     */
    void finalize(const bool &zeroFill);

    /** 
     * Compute alpha*M*x + beta*y, where M=*this.
     * @param M_trans specifies whether the operator is transposed:
     *                op(M) = M, for M_trans == NOTRANS
     *                op(M) = M', for M_trans == TRANS
     * @param x       vector of length this->domain()->dim()
     * @param y       vector of length this->range()->dim()
     * @param alpha   scalar multiplying M*x (default is 1.0)
     * @param beta    scalar multiplying y (default is 0.0)
     */
    virtual void apply(
                       const TSFCore::ETransp            M_trans
                       ,const TSFCore::Vector<Scalar>    &x
                       ,TSFCore::Vector<Scalar>          *y
                       ,const Scalar            //alpha = 1.0
                       ,const Scalar           // beta  = 0.0
                       ) const;
    

    /** apply operator to a vector in the domain space, returning a vector
     * in the range space 
     */
    void apply(const Vector<Scalar>& in, const Scalar alpha, 
	       Vector<Scalar>& out, const Scalar beta) const ;



    /** apply adjoint operator to a vector in the domain space, returning
     * a vector in the range space. The default implementation throws an
     * exception */
    virtual void applyTranspose(const Vector<Scalar>& in, const Scalar alpha, 
			      Vector<Scalar>& out, const Scalar beta) const ;



    /**  create the transpose */
    LinearOperator<Scalar>* formTranspose();

    /**
     * Write to a stream
     */
    void print(ostream& os) const;


    /** Overwrites describe(int depth) to handle the block structure 
     *
     * @param depth int to specify how many tabs to indent the line
     */
    string describe(const int &depth) const;
    

  private:
    ProductVectorSpace<Scalar> domain_;
    ProductVectorSpace<Scalar> range_;
    int nBlockRows_;
    int nBlockCols_;
    bool isUnspecified_;
    bool isFinal_;

    Array<Array<LinearOperator<Scalar> > > sub_;
    Array<Array<int> > isSet_;
    LinearOperator<Scalar> opTrp_;

    /** Private method to build the VectorSpaces when building this
     *  operator on the fly.
     *
     * @param i int specifying the row
     * @param j int specifying the col
     * @param sub LinearOperator to be inserted into the (i, j) position.
     *
     */
    bool buildSpaces(const int &i, const int &j, 
		     const LinearOperator<Scalar>& sub);
    


    /** Private method to check compatibility of the spaces before a
     *  block is set.  
     *
     * @param i int specifying the row
     * @param j int specifying the col
     * @param sub LinearOperator to be inserted into the (i, j) position.
     */
    void chkSpaces(const int &i, const int &j, 
		   const LinearOperator<Scalar>& sub) const;
    


    /** Private method to check that all blocks are set */
    bool chkFinal() const;
    

    /** Private method to put ZeroOperators in the (i, j) location.
     *
     * @param i int specifying the row
     * @param j int specifying the col
     */
    void zeroFill(const int &i, const int &j);
    



  }; 
}

#endif
