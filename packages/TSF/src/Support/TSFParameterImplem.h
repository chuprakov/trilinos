#ifndef TSFPARAMETERIMPLEM_H
#define TSFPARAMETERIMPLEM_H

#include "TSFConfig.h"
#include <string>
#include "TSFArray.h"

namespace TSF 
{

	using namespace std;
	/** 
	 *	TSFParameterImplem:  low-level implementation of TSF parameter class.
	 */
	class TSFParameterImplem 
		{
		public:
			/* Default copy ctor and dtor are used. This is safe because we use
			 * well-behaved containers (string and vector) internally. */

			/** construct with a label */
			TSFParameterImplem(const string& label);
			
			/** TSFParameterImplem character constructor. */
			TSFParameterImplem(const string& label, char parameter);
			
			/** TSFParameterImplem character string constructor. */
			TSFParameterImplem(const string& label, const string& parameter);
			
			/** TSFParameterImplem int constructor. */
			TSFParameterImplem(const string& label, int parameter);
			
			/** TSFParameterImplem int array constructor. */
			TSFParameterImplem(const string& label, const TSFArray<int>& parameter);
			
			/** TSFParameterImplem double constructor. */
			TSFParameterImplem(const string& label, double parameter);
			
			/** TSFParameterImplem double array constructor. */
			TSFParameterImplem(const string& label, const TSFArray<double>& parameter);
			
			/** set character value */
			void set(char parameter);

			/** set string value */
			void set(const string& parameter);

			/** set integer value */
			void set(int parameterImplem);

			/** set int array value */
			void set(const TSFArray<int>& parameter);

			/** set double value */
			void set(double parameter);

			/** set double array value */
			void set(const TSFArray<double>& parameter);

			/** True if parameter is character */
			bool isChar() const { return(isChar_);};
			/** True if parameter is character string */
			bool isCharString() const { return(isCharString_);};
			/** True if parameter is integer */
			bool isInt() const { return(isInt_);};
			/** True if parameter is integer array */
			bool isIntArray() const { return(isIntArray_);};
			/** True if parameter is double */
			bool isDouble() const { return(isDouble_);};
			/** True if parameter is double array */
			bool isDoubleArray() const { return(isDoubleArray_);};

			/** Return character parameter. */
			char getChar() const {mustBeChar(); return(charValue_);};

			/** Return character string parameter. */
			const string& getCharString() const {mustBeCharString(); return charStringValue_;}

			/** Return integer parameter. */
			int  getInt() const {mustBeInt(); return intValue_;}

			/** Return integer array parameter. */
			const TSFArray<int>&  getIntArray() const {mustBeIntArray(); return intArrayValue_;}

			/** Return double parameter. */
			double getDouble() const {mustBeDouble(); return doubleValue_;}

			/** Return double array parameter. */
			const TSFArray<double>& getDoubleArray() const 
				{mustBeDoubleArray(); return doubleArrayValue_;}

			//@}

			/** Returns the label associated with this parameterImplem.  A label is 
			 * a description of the parameterImplem, used
			 * both to describe the variable and identify it to any 
			 * prospective Trilinos solver component.
			 */
			const string& getLabel() const {return label_;}

			/** return a string indicating the type */
			string type() const ;

		private:

			/** set default values */
			void setDefaults(const string& label);

			/** report a call with a mismatched type */
			void wrongType(const string& wrongType) const ;

			/** report an attempt to change type */
			void typeChangeError(const string& attemptedType) const ;

			/* methods to check that a getType() method has been called on the
			 * correct type of parameterImplem */
			inline void mustBeChar() const {if (!isChar()) wrongType("char");}
			inline void mustBeCharString() const {if (!isCharString()) wrongType("string");}
			inline void mustBeInt() const {if (!isInt()) wrongType("int");}
			inline void mustBeDouble() const {if (!isDouble()) wrongType("double");}
			inline void mustBeIntArray() const {if (!isIntArray()) wrongType("int array");}
			inline void mustBeDoubleArray() const 
				{if (!isDoubleArray()) wrongType("double array");}
	

			/* the label */
			string label_;

			/* values. Only one of these will be valid. */
			char charValue_;
			string charStringValue_;
			int intValue_;
			TSFArray<int> intArrayValue_;
			double doubleValue_;
			TSFArray<double> doubleArrayValue_;

			/* flag indicating whether the type is defined. This is in the logic that prevents
			 * changing the type of a parameter after its type has been set 
			 */
			bool typeIsDefined_;

			/* flags that tell us the type */
			bool isChar_;
			bool isCharString_;
			bool isInt_; 
			bool isIntArray_;
			bool isDouble_;
			bool isDoubleArray_;
		};

} 

#endif /* _TSF_PARAMETERIMPLEM_H_ */
