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

#ifndef TSFEPETRAMATRIX_HPP
#define TSFEPETRAMATRIX_HPP

#include "TSFCoreEpetraLinearOp.hpp"
#include "TSFEpetraVectorSpace.hpp"
#include "TSFLoadableMatrix.hpp"
#include "TSFLinearOperatorDecl.hpp"
#include "TSFRowAccessibleOp.hpp"
//#include "TSFExplicitlyTransposeableOp.hpp"
#include "TSFHandleable.hpp"
#include "TSFPrintable.hpp"
#include "TSFDescribableByTypeID.hpp"
#include "TSFILUFactorizableOp.hpp"
#include "Epetra_CrsMatrix.h"

namespace TSFExtended
{
  using namespace Teuchos;

  /** */
  class EpetraMatrix : public TSFCore::EpetraLinearOp,
                       public LoadableMatrix<double>,
		       public RowAccessibleOp<double>,
                       // public ExplicitlyTransposeableOp<double>,
                       public Printable,
                       public DescribableByTypeID,
                       public ILUFactorizableOp<double>,
                       public Handleable<TSFCore::LinearOp<double> >
  {
  public:
    /** Construct an uninitialized EpetraMatrix */
    EpetraMatrix(const RefCountPtr<const EpetraVectorSpace>& domain,
                 const RefCountPtr<const EpetraVectorSpace>& range);

    /** Construct an uninitialized EpetraMatrix */
    EpetraMatrix(const RefCountPtr<const EpetraVectorSpace>& domain,
                 const RefCountPtr<const EpetraVectorSpace>& range,
                 const int* numEntriesPerRow);

    /** Virtual dtor */
    virtual ~EpetraMatrix(){;}

    /** */
    virtual void configure(int lowestRow,
                           const std::vector<std::set<int> >& nonzeros);


    /** */
    virtual void configure(int lowestRow,
                           const std::vector<std::vector<int> >& nonzeros);

    /** */
    virtual void configure(int lowestRow,
                           const std::vector<int>& rowPtrs,
                           const std::vector<int>& nnzPerRow,
                           const std::vector<int>& data);

    /** 
     * Set the locations of all my nonzero elements. 
     * @param nLocalRows number of locally-owned rows
     * @param globalRowIndex array of global indices of the local rows
     * @param numNonzeros array of number of nonzeros for each row
     * @param array of arrays of column indices for each row
     */
    virtual void setGraph(int nLocalRows,
                          const int* globalRowIndex,
                          const int* numNonzeros,
                          const int** columnIndices) ;

    /** Insert a set of elements in a row, overwriting any previously
     * existing values. 
     * @param globalRowIndex the global index of the row to which these
     * elements belong.
     * @param nElemsToInsert the number of elements being inserted in this
     * step
     * @param globalColumnIndices array of column indices. Must 
     * be nElemsToInsert in length. 
     * @param elements array of element values. Must be nElemsToInsert in
     * length
     */
    virtual void setRowValues(int globalRowIndex,
                              int nElemsToInsert,
                              const int* globalColumnIndices,
                              const double* elementValues)  ;

    /** Insert a set of elements in a row, adding to any previously
     * existing values. 
     * @param globalRowIndex the global index of the row to which these
     * elements belong.
     * @param nElemsToInsert the number of elements being inserted in this
     * step
     * @param globalColumnIndices array of column indices. Must 
     * be nElemsToInsert in length. 
     * @param elements array of element values. Must be nElemsToInsert in
     * length
     */
    virtual void addToRow(int globalRowIndex,
                          int nElemsToInsert,
                          const int* globalColumnIndices,
                          const double* elementValues) ;


    /** 
     *
     */
    virtual void addElementBatch(int numRows, 
                                 int rowBlockSize,
                                 const int* globalRowIndices,
                                 int numColumnsPerRow,
                                 const int* globalColumnIndices,
                                 const Scalar* values,
                                 const int* skipRow);

    /** Set all elements to zero, preserving the existing structure */
    virtual void zero() ;

    /** Finalize values of the matrix.  */
    virtual void freezeValues() ;


    /** \name incomplete factorization preconditioning interface */
    //@{
    /** create an incomplete factorization. 
     * @param fillLevels number of levels of fill on the local processor
     * @param overlapFill number of levels of fill on remote processors
     * @param relaxationValue fraction of dropped values to be added to the
     * diagonal
     * @param relativeThreshold relative diagonal perutrbation
     * @param absoluteThreshold absolute diagonal perturbation
     * @param leftOrRight whether this preconditioner is to be applied
     * from the left or right 
     * @param rtn newly created preconditioner, returned 
     * by reference argument.
     */
    virtual void getILUKPreconditioner(int fillLevels,
                                       int overlapFill,
                                       double relaxationValue,
                                       double relativeThreshold,
                                       double absoluteThreshold,
                                       LeftOrRight leftOrRight,
                                       Preconditioner<double>& rtn) const ;

    /** Describable interface */
    virtual string description() const ;

    /** Printable interface */
    virtual void print(ostream& os) const ;

    GET_RCP(TSFCore::LinearOp<double>);


    string describe(int depth) const
    {
      string ret = "";
      for (int i = 0; i < depth; i++)
	{
	  ret.append("   ");
	}
      ret.append(typeName());
      ret.append(" of dimension " + toString(range()->dim()) + " by "
		 + toString(domain()->dim()));
      return ret;
    }



    /** */
    static Epetra_CrsMatrix& getConcrete(const LinearOperator<double>& A);

    /** 
     * Read-only access to the underlying crs matrix. Needed for Ifpack.
     */
    const Epetra_CrsMatrix* crsMatrix() const ;
  protected:
     /** \name Allocators for domain and range spaces */
  //@{
  /** Allocate the domain space of the operator. Purpose: In
   * TSFExtended, both EpetraLinearOp and EpetraVectorSpace are
   * extended from the TSFCore versions by inheritance, and the
   * TSFExtended operator subclasses expect to work with an extended
   * vector space subclass. Thus, it is necessary for the base
   * operator class to never directly allocate vector space objects,
   * and allocation is delegated to a virtual allocator function. 
   * KRL and RAB, 2/18/04. */
  virtual RefCountPtr<const TSFCore::EpetraVectorSpace> 
  allocateDomain(const RefCountPtr<Epetra_Operator>  &op 
                 ,TSFCore::ETransp  op_trans 
                 )  const ; 
  
  /** Allocate the range space of the operator. Purpose: In
   * TSFExtended, both EpetraLinearOp and EpetraVectorSpace are
   * extended from the TSFCore versions by inheritance, and the
   * TSFExtended operator subclasses expect to work with an extended
   * vector space subclass. Thus, it is necessary for the base
   * operator class to never directly allocate vector space objects,
   * and allocation is delegated to a virtual allocator function. 
   * KRL and RAB, 2/18/04. */
  virtual RefCountPtr<const TSFCore::EpetraVectorSpace> allocateRange( 
    const RefCountPtr<Epetra_Operator>  &op 
    ,TSFCore::ETransp  op_trans 
    )  const ; 




    /** Get the specified row as defined by RowAccessible  */
    void getRow(const int& row, 
		Teuchos::Array<int>& indices, 
		Teuchos::Array<Scalar>& values) const;
    


  private:
    Epetra_CrsMatrix* crsMatrix();

    
  };
}

#endif
