/* @HEADER@ */
/* @HEADER@ */

#ifndef TSFEPETRAMATRIX_HPP
#define TSFEPETRAMATRIX_HPP

#include "TSFCoreEpetraLinearOp.hpp"
#include "TSFEpetraVectorSpace.hpp"
#include "TSFLoadableMatrix.hpp"
//#include "TSFRowAccessibleOp.hpp"
//#include "TSFExplicitlyTransposeableOp.hpp"
#include "TSFHandleable.hpp"
#include "TSFPrintable.hpp"
#include "TSFDescribable.hpp"
#include "Epetra_CrsMatrix.h"

namespace TSFExtended
{
  using namespace Teuchos;

  /** */
  class EpetraMatrix : public TSFCore::EpetraLinearOp,
                       public LoadableMatrix<double>,
                       //   public RowAccessibleOp<double>,
                       // public ExplicitlyTransposeableOp<double>,
                       public Printable,
                       public Describable,
                       public Handleable<TSFCore::LinearOp<double> >
  {
  public:
    /** Construct an uninitialized EpetraMatrix */
    EpetraMatrix(const RefCountPtr<const EpetraVectorSpace>& domain,
                 const RefCountPtr<const EpetraVectorSpace>& range);

    /** Virtual dtor */
    virtual ~EpetraMatrix(){;}

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

    /** Set all elements to zero, preserving the existing structure */
    virtual void zero() ;

    /** Finalize values of the matrix.  */
    virtual void freezeValues() ;

    /** Describable interface */
    virtual string describe() const ;

    /** Printable interface */
    virtual void print(ostream& os) const ;

    /** Handleable interface */
    virtual RefCountPtr<TSFCore::LinearOp<double> > getRcp()
    {return rcp(this);}

  private:
    Epetra_CrsMatrix* crsMatrix();

    const Epetra_CrsMatrix* crsMatrix() const ;
    
  };
}

#endif
