#ifndef TSFINCOREVECTOR_H
#define TSFINCOREVECTOR_H

#include "TSFDefs.h"
#include "TSFVector.h"
#include "TSFAccessibleVector.h"
#include <typeinfo>

namespace TSF
{

  using std::string;

  /** \ingroup VectorSubtypes
   * TSFInCoreVector is a base class for vector implementations that have vector
   * values living as a contiguous array in memory. We can access the array with
   * a single call to getNextChunk().
   */

  class TSFInCoreVector : public TSFAccessibleVector
    {
    public:
      /** construct with a given space */
      TSFInCoreVector(const TSFVectorSpace& space);

      /** the usual virtual dtor */
      virtual ~TSFInCoreVector(){;}

      /** \name Accessing chunks of elements */
      //@{
      /** Set the chunk iterator back to start */
      virtual void rewindChunkIterator() const {hasMoreChunks_ = true;}

      /** Say whether there are more chunks to be processed */
      virtual bool hasMoreChunks() const {return hasMoreChunks_;}

      /** Get the next read-only chunk, returning the size via reference argument */
      virtual const TSFReal* getNextChunk(int& chunkSize) const ;

      /** Get the next read-write chunk, returning the size via reference argument */
      virtual TSFReal* getNextChunk(int& chunkSize) ;

      /** Get the next read-only chunk of on-processor elements,
       * returning the size via reference argument */
      virtual const TSFReal* getNextLocalChunk(int& chunkSize) const ;

      /** Restore a chunk after use. This is needed for out-of-core vectors and
       * some in-core types such as PETSC */
      virtual void restoreChunk(const TSFReal* /* chunk */,
                                int /* chunkSize */) const {;}
      //@}

    protected:

      /** Read-only access to the vector data */
      virtual const TSFReal* dataPointer() const = 0 ;

      /** Read-write access to the vector data */
      virtual TSFReal* dataPointer() = 0 ;

      /** Return number of local elements */
      virtual int nLocal() const = 0 ;

      /** Return number of off-processor elements copied locally */
      virtual int nGhost() const = 0 ;

      /** */
      mutable bool hasMoreChunks_;

      /** */
      static TSFReal dummyElement_;
    };
};

#endif
