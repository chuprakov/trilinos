#ifndef TSFBLOCKVECTOR_H
#define TSFBLOCKVECTOR_H

#include "TSFConfig.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"

namespace TSF
{
	/** \ingroup VectorSubtypes
	 * Block vectors
	 */

	class TSFBlockVector : public TSFVectorBase
		{
		public: 
			/** construct with a given space */
			TSFBlockVector(const TSFVectorSpace& space);

			/** the usual virtual dtor */
			virtual ~TSFBlockVector(){;}
	
			/** \name Block access */
			//@{
			/** return the number of subvector blocks */
			virtual int numBlocks() const ;

            /** states that I am a BlockVector */
            virtual bool isBlockVector() const
              {
                return true;
              }

			/** return the i-th subvector */
			virtual void getBlock(int i, const TSFVector& self, 
														TSFVector& sub) const ;

			/** set the i-th subvector */
			virtual void setBlock(int i, const TSFVector& sub);
			//@}

			/** \name mathematical methods */
			//@{
			/** x = a*x + y */
			virtual void axpy(const TSFReal& a, const TSFVector& y) ;

			/** multiply by a scalar */
			virtual void scalarMult(const TSFReal& a) ;

			/** element-by-element multiplication (this = this .* y) */
			virtual void dotStar(const TSFVector& x) ;

			/** element-by-element division (this = this ./ y) */
			virtual void dotSlash(const TSFVector& x) ; 

			/** dot product */
			virtual TSFReal dot(const TSFVector& other) const ;

			/** infinity-norm */
			virtual TSFReal normInf() const ;

			/** one-norm */
			virtual TSFReal norm1() const ;

			/** sum all elements */
			virtual TSFReal sumElements() const ;

			/** set all elements to a scalar value */
			virtual void setScalar(const TSFReal& a) ;

			/** find an extreme value, possibly subject to a constraint */
			virtual TSFReal findExtremeValue(MinOrMax type, TSFGeneralizedIndex& location, 
																			 const TSFReal& tol) const ;

			/** elementwise absolute value */
			virtual void abs() ;
			
			//@}	


		
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
					
			/** \name access to elements */
			//@{
			/** read-only access to a single element */
			virtual const TSFReal& getElement(int globalIndex) const ;

			/** read-write access to a single element */
			virtual TSFReal& setElement(int globalIndex) ;

			/** set a block of elements */
			virtual void setElements(int n, const int* globalIndices, 
															 const TSFReal* values) ;

			/** get a block of elements */
			virtual void getElements(int n, const int* globalIndices, 
															 TSFReal* values) const ;
			
			/** add to a block of elements */
			virtual void addToElements(int n, const int* globalIndices, 
																 const TSFReal* values) ;
			//@}	

		

			/** \name generating random vectors */
			//@{
			/** Fill a vector with random elements */
			virtual void randomize(const TSFRandomNumberGenerator& r) ;
			//@}	

			/** \name Hooks for parallel support */
			//@{
			/** gather valid ghost values from other procs. Default is a no-op */
			virtual void synchronizeGhostValues() const ;

			/** mark ghost values as invalid, meaning that they need to be
			 * synchronized. Default is a no-op.  */
			virtual void invalidateGhostValues() ;
			//@}

			/** \name maintenance methods */
			//@{
			/** virtual copy ctor */
			virtual TSFVectorBase* deepCopy() const ;

			/** copy the data from another vector into this vector */
			virtual void acceptCopyOf(const TSFVector& x) ;

			/** write to stream */
			virtual ostream& print(ostream& os) const ;
			//@}
		protected:
			/* array of subvectors */
			TSFArray<TSFVector> subvectors_;
		private:
#ifdef HAVE_RTOP
			// Note: *this is first element of vec[] or targ_vec[]
			void apply_op(
				const RTOp_RTOp &op, int num_vecs, const TSFBlockVector* vecs[]
				,int num_targ_vecs, TSFBlockVector* targ_vecs[], RTOp_ReductTarget reduct_obj
				,const RTOp_index_type first_ele, const RTOp_index_type sub_dim, const RTOp_index_type global_offset
				) const;
#endif
			
		};
};

#endif
