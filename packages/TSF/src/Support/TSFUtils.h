#ifndef TSFUTILS_H
#define TSFUTILS_H

#include "TSFDefs.h"
#include <math.h>
#include <string>

namespace TSF
{
  using std::string;

  /**\ingroup Utilities
   *
   * Numerical constants, etc.
   */
  class TSFUtils
    {
    public:
      /**
       * print a description of the current build
       */
      static void aboutBuild();

      /**
       * Set a number to zero if it is less than chopVal.
       */
      static TSFReal chop(const TSFReal& x);

      /**
       * write a real as a string
       */
      static string toString(const TSFReal& x);

      /**
       * write an int as a string
       */
      static string toString(const int& x);

      /**
       * IEEE positive infinity
       */
      static TSFReal infinity() {return HUGE_VAL;}

      /**
       * IEEE negative infinity
       */
      static TSFReal negativeInfinity() {return -HUGE_VAL;}

      /**
       * pi.
       */
      static TSFReal pi() {return M_PI;}

      /**
       * Get the chopping value, below which numbers are considered to be zero
       */
      static TSFReal getChopVal() {return chopVal_;}
      /**
       * Set the chopping value, below which numbers are considered to be zero
       */
      static void setChopVal(TSFReal chopVal) {chopVal_ = chopVal;}

    private:
      static TSFReal chopVal_;
    };

  /** \relates TSFUtils */
  inline string toString(const int& x) {return TSFUtils::toString(x);}

  /** \relates TSFUtils */
  inline string toString(const TSFReal& x) {return TSFUtils::toString(x);}

  /** \relates TSFUtils */
  inline string toString(const string& x) {return x;}

}

#endif


