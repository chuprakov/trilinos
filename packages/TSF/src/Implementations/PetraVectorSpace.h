#ifndef PETRAVECTORSPACE_H
#define PETRAVECTORSPACE_H

#include "TSFConfig.h"

#if HAVE_PETRA

#define PETRA_BOOL_SUPPORTED
#include "TSFMPIVectorSpace.h"
#include "Epetra_Import.h"
#include "Epetra_Map.h"
#include "TSFSmartPtr.h"

namespace TSF
{
	using std::string;

	/** \ingroup Petra
	 * Base class for vector spaces 
	 */

	class PetraVectorSpace : public TSFMPIVectorSpace
		{
		public: 
			/**
			 * construct with a Epetra_Map specifying the dimension and processor
			 * distribution of the elements.
			 */ 
			PetraVectorSpace(const TSFSmartPtr<Epetra_Map>& localMap);
			/**
			 * construct with a Epetra_Map specifying the dimension and processor
			 * distribution of the elements.
			 */ 
			PetraVectorSpace(const TSFSmartPtr<Epetra_Map>& localMap,
											 const TSFSmartPtr<Epetra_Map>& ghostMap,
											 const TSFSmartPtr<Epetra_Import>& importer);

			/** the usual virtual dtor */
			virtual ~PetraVectorSpace(){;}
	
			/** virtual copy ctor */
			virtual TSFVectorSpaceBase* deepCopy() const ;
	
			/** return dimension of space */
			virtual int dim() const ;

			/** create a vector that is a member of this space */
			virtual TSFVectorBase* createMember(const TSFVectorSpace& handle) const;

			/** write to stream */
			virtual ostream& print(ostream& os) const ;

			/** get the petra map from a vector space. Throw an exception
			 * if the space is not a petra vector space */
			static TSFSmartPtr<Epetra_Map> getLocalMap(const TSFVectorSpace& space);

			static TSFSmartPtr<Epetra_Map> getMap(const TSFVectorSpace& space);

			static TSFSmartPtr<Epetra_Import> getImporter(const TSFVectorSpace& space);
		protected:
			static const PetraVectorSpace* getPtr(const TSFVectorSpace& space);

			TSFSmartPtr<Epetra_Map> localMap_;
			TSFSmartPtr<Epetra_Map> ghostMap_;
			TSFSmartPtr<Epetra_Import> importer_;
		};
	
}

#endif // PETRA
#endif
