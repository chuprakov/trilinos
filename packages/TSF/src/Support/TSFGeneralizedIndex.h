#ifndef TSFGENERALIZEDINDEX_H
#define TSFGENERALIZEDINDEX_H

#include "TSFConfig.h"
#include "TSFError.h"
#include "TSFStack.h"

namespace TSF
{
	using std::exception;
	using std::string;

	/** 
	 *
	 */
	class TSFGeneralizedIndex
		{
		public:
			/** */
			TSFGeneralizedIndex(const TSFGeneralizedIndex& index, int i);

			/** */
			TSFGeneralizedIndex(int i);

			/** */
			TSFGeneralizedIndex();

			/** */
			int index() const {return indices_.peek();}

			/** */
			TSFGeneralizedIndex remainder() const ;

		private:
			/** */
			TSFGeneralizedIndex(const TSFStack<int>& stack);

			
			/** */
			mutable TSFStack<int> indices_;

			
			
		};
}
#endif
