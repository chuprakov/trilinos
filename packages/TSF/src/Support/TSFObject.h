#ifndef TSFOBJECT_H
#define TSFOBJECT_H

#include "TSFDefs.h"

namespace TSF
{

  using std::string;
  using std::ostream;

  /**
   * TSFObject establishes the baseline interface for most classes in TSF. All classes
   * should define print() and hashCode() methods.
   */
  class TSFObject
    {
    public:
      /** Empty ctor */
      TSFObject();

      /** Virtual dtor. All classes derived from this should define virtual dtors. */
      virtual ~TSFObject();

      /** Print a short identifier to a stream */
      virtual void print(ostream& os) const  = 0;

      /** Print a short identifier to a stream, indented by indentDepth */
      virtual void printIndented(ostream& os, int indentDepth) const ;

      /** Print a verbose description to a stream. Default implementation is to call
       * print() */
      virtual void describe(ostream& os) const  ;

      /** Print to a string */
      virtual string toString() const ;

      /** Return a string identifying this object's type */
      virtual string typeName() const ;

      /** */
      static string spaces(int indentDepth);
    protected:
    private:
    };

  /** \relates TSFObject
   * Print to a stream
   */
  inline ostream& operator<<(ostream& os, const TSFObject& obj)
    {
      obj.print(os);
      return os;
    }

  /** \relates TSFObject
   * Write to a string
   */
  inline string toString(const TSFObject& obj)
    {
      return obj.toString();
    }

}


#endif
