#ifndef TSFEPETRAVECTOR_HPP
#define TSFEPETRAVECTOR_HPP

#include "TSFConfigDefs.hpp"
#include "TSFPrintable.hpp"
#include "TSFDescribable.hpp"
#include "TSFCoreEpetraVector.hpp"
#include "TSFAccessibleVector.hpp"
#include "TSFLoadableVector.hpp"
#include "TSFVector.hpp"
#include "Epetra_FEVector.h"
#include "Epetra_Vector.h"


namespace TSFExtended
{
  using TSFCore::Index;
  using namespace Teuchos;
  /**
   * TSF extension of TSFCore::EpetraVector, implementing the LoadableVector
   * interface allowing an application to access elements. This class derives
   * from TSFCore::EpetraVector, so it can be used seamlessly in any 
   * TSFCore-based code.
   */
  class EpetraVector : public TSFCore::EpetraVector,
                       public LoadableVector<double>,
                       public AccessibleVector<double>,
                       public Describable,
                       public Printable
    {
    public:
      /** Construct with a smart pointer to an Epetra FE vector. */
      EpetraVector(const RefCountPtr<Epetra_Vector>& vec,
                   const RefCountPtr<const TSFCore::EpetraVectorSpace>& map);

      /** virtual dtor */
      virtual ~EpetraVector() {;}

      /** \name LoadableVector interface */
      //@{
      /** set a single element */
      void setElement(Index globalIndex, const double& value);

      /** add to a single element */
      void addToElement(Index globalIndex, const double& value);

      /** set a group of elements */
      void setElements(size_t numElems, const Index* globalIndices, 
                       const double* values);


      /** add to a group of elements */
      void addToElements(size_t numElems, const Index* globalIndices, 
                         const double* values);

      /** */
      void finalizeAssembly();
      //@}

      /** \name AccessibleVector interface */
      //@{
      /** */
      const double& getElement(Index globalIndex) const ;
      //@}
      

      /** \name Printable interface */
      //@{
      /** Write to a stream  */
      void print(ostream& os) const 
      {
        epetra_vec()->Print(os);
      }
      //@}

      /** \name Describable interface */
      //@{
      /** Write a brief description */
      string describe() const 
      {
        return "EpetraVector";
      }
      //@}


      /** Get a read-only Epetra_Vector */
      static const Epetra_Vector& getConcrete(const Vector<double>& tsfVec);
      /** Get a read-write Epetra_Vector */
      static Epetra_Vector& getConcrete(Vector<double>& tsfVec);
    };
  
}

#endif
