#ifndef TSFSERIALVECTOR_H
#define TSFSERIALVECTOR_H

#include "TSFConfig.h"
#include "TSFInCoreVector.h"
#include "DenseSerialVector.h"
#include <typeinfo>


namespace TSF
{
	
	using std::string;

	/** \ingroup ConcreteVectors
	 * A vector contained on a single processor.
	 */

	class TSFSerialVector : public TSFInCoreVector
		{
		public: 
			/** construct with a given space */
			TSFSerialVector(const TSFVectorSpace& space);

			/** the usual virtual dtor */
			virtual ~TSFSerialVector(){;}


            /* return type name */
            string name(){return typeid(*this).name();}

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
			virtual TSFReal& setElement(int globalIndex) {return x_[globalIndex];}
			
			/** read-only access */
			virtual const TSFReal& getElement(int globalIndex) const {return x_[globalIndex];}

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

			/** write to stream */
			virtual ostream& print(ostream& os) const ;
			//@}


			/** read-only access */
			static const DenseSerialVector& getConcrete(const TSFVector& x) ;
			/** read-write access */
			static DenseSerialVector& getConcrete(TSFVector& x) ;
			
		protected:

			/** return a read-only pointer to the vector data */
			virtual const TSFReal* dataPointer() const {return &(x_[0]);}

			/** return a writeable pointer to the vector data */
			virtual TSFReal* dataPointer() {return &(x_[0]);}

			/** return the number of local elements, which is equal to the total number of 
			 * elements for a serial vector */
			virtual int nLocal() const {return x_.length();}

			/** return the number of ghost points, which is zero for a serial vector */
			virtual int nGhost() const {return 0;}

			/** \name reduction methods */
			/** Reduce a calculation across processors.  */
			TSFReal reduce(const TSFReal& localValue, TSFReductionOp /* op */) const 
				{return localValue;}

			DenseSerialVector x_;

		private:
#ifdef HAVE_RTOP
  void apply_op(
	  const TSFSerialVector* const_this, TSFSerialVector* nonconst_this
	  ,const RTOp_RTOp& op
	  ,const int num_vecs,      const TSFVectorBase**   vecs
	  ,const int num_targ_vecs, TSFVectorBase**         targ_vecs
	  ,RTOp_ReductTarget reduct_obj
	  ,const int first_ele  , const int sub_dim  , const int global_offset
	  ) const;
#endif // HAVE_RTOP

		};
};

#endif
