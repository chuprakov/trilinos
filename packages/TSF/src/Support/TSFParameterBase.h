#ifndef TSFPARAMETERBASE_H
#define TSFPARAMETERBASE_H

#include "TSFConfig.h"
#include <string>
#include "TSFArray.h"

namespace TSF
{

  using namespace std;
  /**
   *  TSFParameterImplem:  low-level implementation of TSF parameter class.
   */
  class TSFParameterBase
    {
    public:
      /** construct with a label */
      TSFParameterBase(const string& label);

      /** virtual dtor */
      virtual ~TSFParameterBase(){}

      /** set character value */
      virtual void set(char parameter);

      /** set string value */
      virtual void set(const string& parameter);

      /** set integer value */
      virtual void set(int parameterImplem);

      /** set int array value */
      virtual void set(const TSFArray<int>& parameter);

      /** set double value */
      virtual void set(double parameter);

      /** set double array value */
      virtual void set(const TSFArray<double>& parameter);

      /** True if parameter is character */
      virtual bool isChar() const {return false;}
      /** True if parameter is character string */
      virtual bool isCharString() const {return false;}
      /** True if parameter is integer */
      virtual bool isInt() const {return false;}
      /** True if parameter is integer array */
      virtual bool isIntArray() const {return false;}
      /** True if parameter is double */
      virtual bool isDouble() const {return false;}
      /** True if parameter is double array */
      virtual bool isDoubleArray() const {return false;}

      /** Return character parameter. */
      virtual char getChar() const ;

      /** Return character string parameter. */
      virtual string getCharString() const ;

      /** Return integer parameter. */
      virtual int  getInt() const ;

      /** Return integer array parameter. */
      virtual TSFArray<int>  getIntArray() const ;

      /** Return double parameter. */
      virtual double getDouble() const ;

      /** Return double array parameter. */
      virtual TSFArray<double> getDoubleArray() const ;

      /** Return string representation of value */
      virtual string getValueString() const = 0 ;


      /** Returns the label associated with this parameterImplem.  A label is
       * a description of the parameterImplem, used
       * both to describe the variable and identify it to any
       * prospective Trilinos solver component.
       */
      const string& getLabel() const {return label_;}

      /** return a string indicating the type */
      virtual string type() const = 0 ;

      /** virtual copy ctor */
      virtual TSFParameterBase* clone() const = 0 ;

    protected:
      /** report a call to a get method with a mismatched type */
      void typeAccessError(const string& wrongType) const ;

      /** report a call to a set method with a mismatched type */
      void typeChangeError(const string& wrongType) const ;

      /* the label */
      string label_;

    };





  /** \ingroup Support
   * Character parameter type
   */
  class TSFCharParameter : public TSFParameterBase
    {
    public:
      /** construct with a label and value */
      TSFCharParameter(const string& label, char value);

      /** virtual dtor */
      virtual ~TSFCharParameter(){;}

      /** set the value */
      virtual void set(char value) {value_ = value;}

      /** inform the world that I am a character */
      virtual bool isChar() const {return true;}

      /** return my value */
      virtual char getChar() const {return value_;}

      /** return type name */
      virtual string type() const {return "char";}

      /** Return string representation of value */
      virtual string getValueString() const {return string(1, value_);}

      /** virtual copy ctor */
      virtual TSFParameterBase* clone() const {return new TSFCharParameter(*this);}

    private:
      char value_;
    };



  /** \ingroup Support
   * Character string parameter type
   */
  class TSFCharStringParameter : public TSFParameterBase
    {
    public:
      /** construct with a label and value */
      TSFCharStringParameter(const string& label, const string& value);

      /** virtual dtor */
      virtual ~TSFCharStringParameter(){;}

      /** set the value */
      virtual void set(const string& value) {value_ = value;}

      /** inform the world that I am a character string */
      virtual bool isCharString() const {return true;}

      /** return my value */
      virtual string getCharString() const {return value_;}

      /** return type name */
      virtual string type() const {return "char string";}

      /** Return string representation of value */
      virtual string getValueString() const {return value_;}

      /** virtual copy ctor */
      virtual TSFParameterBase* clone() const {return new TSFCharStringParameter(*this);}

    private:
      string value_;
    };


  /** \ingroup Support
   * Integer parameter type
   */
  class TSFIntParameter : public TSFParameterBase
    {
    public:
      /** construct with a label and value */
      TSFIntParameter(const string& label, int value);

      /** virtual dtor */
      virtual ~TSFIntParameter(){;}

      /** set the value */
      virtual void set(int value) {value_ = value;}

      /** inform the world that I am an integer */
      virtual bool isInt() const {return true;}

      /** return my value */
      virtual int getInt() const {return value_;}

      /** return type name */
      virtual string type() const {return "int";}

      /** Return string representation of value */
      virtual string getValueString() const {return TSF::toString(value_);}

      /** virtual copy ctor */
      virtual TSFParameterBase* clone() const {return new TSFIntParameter(*this);}

    private:
      int value_;
    };


  /** \ingroup Support
   * Integer array parameter type
   */
  class TSFIntArrayParameter : public TSFParameterBase
    {
    public:
      /** construct with a label and value */
      TSFIntArrayParameter(const string& label, const TSFArray<int>& value);

      /** virtual dtor */
      virtual ~TSFIntArrayParameter(){;}

      /** set the value */
      virtual void set(const TSFArray<int>& value) {value_ = value;}

      /** inform the world that I am an integer array */
      virtual bool isIntArray() const {return true;}

      /** return my value */
      virtual TSFArray<int> getIntArray() const {return value_;}

      /** return type name */
      virtual string type() const {return "int array";}

      /** Return string representation of value */
      virtual string getValueString() const {return TSF::toString(value_);}

      /** virtual copy ctor */
      virtual TSFParameterBase* clone() const {return new TSFIntArrayParameter(*this);}

    private:
      TSFArray<int> value_;
    };


  /** \ingroup Support
   * Double parameter type
   */
  class TSFDoubleParameter : public TSFParameterBase
    {
    public:
      /** construct with a label and value */
      TSFDoubleParameter(const string& label, const double& value);

      /** virtual dtor */
      virtual ~TSFDoubleParameter(){;}

      /** set the value */
      virtual void set(double value) {value_ = value;}

      /** inform the world that I am a double */
      virtual bool isDouble() const {return true;}

      /** return my value */
      virtual double getDouble() const {return value_;}

      /** return type name */
      virtual string type() const {return "double";}

      /** Return string representation of value */
      virtual string getValueString() const {return TSF::toString(value_);}

      /** virtual copy ctor */
      virtual TSFParameterBase* clone() const {return new TSFDoubleParameter(*this);}

    private:
      double value_;
    };

  /** \ingroup Support
   * Double array parameter type
   */
  class TSFDoubleArrayParameter : public TSFParameterBase
    {
    public:
      /** construct with a label and value */
      TSFDoubleArrayParameter(const string& label, const TSFArray<double>& value);

      /** virtual dtor */
      virtual ~TSFDoubleArrayParameter(){;}

      /** set the value */
      virtual void set(const TSFArray<double>& value) {value_ = value;}

      /** inform the world that I am an double array */
      virtual bool isDoubleArray() const {return true;}

      /** return my value */
      virtual TSFArray<double> getDoubleArray() const {return value_;}

      /** return type name */
      virtual string type() const {return "double array";}

      /** Return string representation of value */
      virtual string getValueString() const {return TSF::toString(value_);}

      /** virtual copy ctor */
      virtual TSFParameterBase* clone() const {return new TSFDoubleArrayParameter(*this);}

    private:
      TSFArray<double> value_;
    };
}

#endif /* TSFPARAMETERBASE_H */
