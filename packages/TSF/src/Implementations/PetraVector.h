#ifndef PETRAVECTOR_H
#define PETRAVECTOR_H

#include "TSFDefs.h"
#include "TSFVectorSpace.h"
#include "TSFMPIVector.h"
#include "PetraVectorSpace.h"
#include "DenseSerialVector.h"
#include <typeinfo>

#if HAVE_PETRA

#include "Epetra_Vector.h"




namespace TSF
{

  using std::string;


  /** \ingroup Petra
   * Vector subtype for Petra.
   *
   * A PetraVector object stores two Epetra_Vector objects
   * providing two representations of the same vector: a complete
   * version allowing access to local values plus ghost points, and
   * also a local version allowing access to local values only. For
   * efficiency, the local version is created with Petra's View mode,
   * that is, it is a view of the same chunk of memory used by the
   * complete version.
   */
  class PetraVector : public TSFMPIVector
    {
    public:
      /** construct with a specification of on-processor rows. Ghost values
       * are considered invalid if this construction is used.
       */
      PetraVector(const TSFVectorSpace& space,
                  const TSFSmartPtr<Epetra_Vector>& localValues);


      /** construct with a specification of on-processor rows plus remotely-owned
       * ghost rows. */
      PetraVector(const TSFVectorSpace& space,
                  const TSFSmartPtr<Epetra_Vector>& localValues,
                  const TSFSmartPtr<Epetra_Vector>& allValues,
                  bool ghostValuesAreValid=true);

      /** the usual virtual dtor */
      virtual ~PetraVector(){;}

#ifdef HAVE_RTOP
      ///
        void apply_reduction(
                             const RTOp_RTOp &op, int num_vecs, const TSFVectorBase* vecs[]
                             ,int num_targ_vecs, TSFVectorBase* targ_vecs[], RTOp_ReductTarget reduct_obj
                             ,const RTOp_index_type first_ele = 1, const RTOp_index_type sub_dim = 0, const RTOp_index_type global_offset = 0
                             ) const;
        ///
          void apply_transformation(
                                    const RTOp_RTOp &op, int num_vecs, const TSFVectorBase* vecs[]
                                    ,int num_targ_vecs, TSFVectorBase* targ_vecs[], RTOp_ReductTarget reduct_obj
                                    ,const RTOp_index_type first_ele = 1, const RTOp_index_type sub_dim = 0, const RTOp_index_type global_offset = 0
                                    );
#endif

          /** \name Element access */
          //@{
          /** read-write access */
          virtual TSFReal& setElement(int globalIndex) ;

          /** read-only access */
          virtual const TSFReal& getElement(int globalIndex) const ;

          /** set a block of elements */
          virtual void setElements(int n, const int* globalIndices,
                                   const TSFReal* values) ;

          /** get a block of elements */
          virtual void getElements(int n, const int* globalIndices, TSFReal* values) const ;

          /** add to a block of elements */
          virtual void addToElements(int n, const int* globalIndices,
                                     const TSFReal* values) ;
          //@}

          /** \name maintenance methods */
          //@{
          /** virtual copy ctor */
          virtual TSFVectorBase* deepCopy() const ;

          /** print, called by stream output */
          virtual ostream& print(ostream& os) const ;


          //@}

          /** gather valid ghost values from other procs */
          virtual void synchronizeGhostValues() const ;

          /** mark ghost values as invalid, meaning that they need to be
           * synchronized */
          virtual void invalidateGhostValues() {ghostValuesAreValid_=false;}

          /** extract the local Epetra_Vector from a TSF vector, for use in mvmults */
          static const Epetra_Vector& getLocalValues(const TSFVector& v);
          /** extract the local Epetra_Vector from a TSF vector, for use in mvmults */
          static Epetra_Vector& getLocalValues(TSFVector& v);


          /** extract the derived type from a TSFVector */
          static const PetraVector& getConcrete(const TSFVector& v);

          /** extract the derived type from a TSFVector */
          static PetraVector& getConcrete(TSFVector& v);

          /** */
          static const TSFSmartPtr<Epetra_Vector>& getEpetraVector(const TSFVector& v);

          friend class AZTECSolver;

    protected:

          /** return a read-only pointer to the vector data */
          virtual const TSFReal* dataPointer() const {return &((*allValues_)[0]);}

          /** return a writeable pointer to the vector data */
          virtual TSFReal* dataPointer() {return &((*allValues_)[0]);}

          /** return the number of local elements, which is equal to the total number of
           * elements for a serial vector */
          virtual int nLocal() const {return getLocalMap().NumMyElements();}

          /** return the number of ghost points */
          virtual int nGhost() const {return getMap().NumMyElements() - nLocal();}

    private:


          /** vector object representing all points, local and ghost. */
          TSFSmartPtr<Epetra_Vector> allValues_;

          /** Vector object representing local points only. This vector is constructed with
           * Petra's view mode as a view of the allValues_ vector. */
          TSFSmartPtr<Epetra_Vector> localValues_;

          /** flag to indicate whether ghost points are currently valid. */
          mutable bool ghostValuesAreValid_;

          /** helper function that tests return value of a petra function */
          void petraCheck(const string& method, int code) const ;

          /** returns the petra vector object containing the local values */
          const Epetra_Vector& localValues() const {return *localValues_;}

          /** returns the petra vector object containing the local values */
          Epetra_Vector& localValues() {return *localValues_;}

          /** returns the petra vector object to be used in a unary vector operation: local if
           * ghost points are not valid, all values otherwise */
          const Epetra_Vector& values() const ;

          /** returns the petra vector object to be used in a unary vector operation: local if
           * ghost points are not valid, all values otherwise */
          Epetra_Vector& values() ;

          /** returns the petra vector object to be used in a binary vector operation: local if
           * ghost points are not valid, all values otherwise */
          const Epetra_Vector& values(bool otherValuesAreValid) const ;

          /** returns the petra vector object to be used in a binary vector operation: local if
           * ghost points are not valid, all values otherwise */
          Epetra_Vector& values(bool otherValuesAreValid) ;

          /** returns an importer that can be used to gather ghost values */
          const Epetra_Import& importer() const ;

          /** returns a map that specifies the locally-owned rows */
          const Epetra_Map& getLocalMap() const ;

          /** returns a map that specifies both local and ghost rows */
          const Epetra_Map& getMap() const ;

          /** maps a global index to a local index */
          int getLocalIndex(int globalIndex) const ;

    };
};

#endif /* HAVE_PETRA */

#endif
