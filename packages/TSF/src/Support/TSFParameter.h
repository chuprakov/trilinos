#ifndef TSFPARAMETER_H
#define TSFPARAMETER_H

#include "TSFDefs.h"
#include "TSFObject.h"
#include "TSFParameterBase.h"
#include "TSFUnsharedPtr.h"
#include "TSFHashUtils.h"
#include <string>
#include "TSFArray.h"

namespace TSF
{


  using std::ostream;
  using std::string;
  /** \ingroup Support
   *  TSFParameter:  The Trilinos Solver Framework Parameter class.
   *  The Parameter class encapsulates information about solver parameters.
   *  Any common type of data can be used to construct an Parameter object.
   *  There is also a constructor that accepts a string that can be parsed into
   *  parameter object.
   *
   * While you can freely change the <I> value </I> of a parameter,
   * it is an error to try to change the <I> type </I> of a parameter. An error
   * will result if you call a set() method having a different type
   * than the type with which the parameter was constructed.
   */
  class TSFParameter : public TSFObject
    {

    public:
      //@{ \name Constructors

      /* Default copy ctor and dtor are used. */

      /** TSFParameter default constructor.
       * Needed for use in containers.
       */
      TSFParameter(void);

      /** Virtual dtor */
      virtual ~TSFParameter(){;}

      /** construct with a label */
      TSFParameter(const string& label);

      /** TSFParameter character constructor. */
      TSFParameter(const string& label, char parameter);

      /** TSFParameter character string constructor. */
      TSFParameter(const string& label, const string& parameter);

      /** TSFParameter int constructor. */
      TSFParameter(const string& label, int parameter);

      /** TSFParameter int array constructor. */
      TSFParameter(const string& label, const TSFArray<int>& parameter);

      /** TSFParameter double constructor. */
      TSFParameter(const string& label, double parameter);

      /** TSFParameter double array constructor. */
      TSFParameter(const string& label, const TSFArray<double>& parameter);
      //@}

      //@{ \name Methods to set parameter values
      /**
       * It is an error to try to change the type of a parameter. An error will
       * result if you call a set() method having a different type than the type with
       * which the parameter was constructed.
       */

      /*
       * NOTE: The original version had label arguments for set methods. This is
       * dangerous, since if we could change the labels it would be possible to
       * put these guys into a hashtable keyed by labels, and then change the labels
       * later resulting in an inconsistency between label and key. A better way is
       * to make the label immutable, setting it only at ctor time.
       * KL Jan 16, 2002.
       */

      /** set character value */
      void set(char parameter) {ptr_->set(parameter);}

      /** set string value */
      void set(const string& parameter) {ptr_->set(parameter);}

      /** set integer value */
      void set(int parameter) {ptr_->set(parameter);}

      /** set int array value */
      void set(const TSFArray<int>& parameter) {ptr_->set(parameter);}

      /** set double value */
      void set(double parameter) {ptr_->set(parameter);}

      /** set double array value */
      void set(const TSFArray<double>& parameter) {ptr_->set(parameter);}

      //@}

      //@{ \name Methods to test parameter type
      /** True if parameter is character */
      bool isChar() const {return ptr_->isChar();}
      /** True if parameter is character string */
      bool isCharString() const {return ptr_->isCharString();}
      /** True if parameter is integer */
      bool isInt() const {return ptr_->isInt();}
      /** True if parameter is integer array */
      bool isIntArray() const {return ptr_->isIntArray();}
      /** True if parameter is double */
      bool isDouble() const {return ptr_->isDouble();}
      /** True if parameter is double array */
      bool isDoubleArray() const {return ptr_->isDoubleArray();}
      //@}

      //@{ \name Methods to extract parameter values

      /** Return character parameter. */
      char getChar() const {return ptr_->getChar();}

      /** Return character string parameter. */
      string getCharString() const {return ptr_->getCharString();}

      /** Return integer parameter. */
      int  getInt() const {return ptr_->getInt();}

      /** Return integer array parameter. */
      TSFArray<int>  getIntArray() const {return ptr_->getIntArray();}

      /** Return double parameter. */
      double getDouble() const {return ptr_->getDouble();}

      /** Return double array parameter. */
      TSFArray<double> getDoubleArray() const {return ptr_->getDoubleArray();}

      //@}

      /** Returns the label associated with this parameter.  A label is
       * a description of the parameter, used
       * both to describe the variable and identify it to any
       * prospective Trilinos solver component.
       */
      const string& getLabel() const {return ptr_->getLabel();}

      /** return a string indicating the type */
      string getType() const {return ptr_->type();}

      /** */
      string getValueString() const {return ptr_->getValueString();}

      /** */
      virtual void print(ostream& os) const ;
    private:

      TSFUnsharedPtr<TSFParameterBase, TSFCloner<TSFParameterBase> > ptr_;
    };

}

#endif /* _TSF_PARAMETER_H_ */
