/* @HEADER@ */
/* @HEADER@ */

#ifndef TSFDOUBLE_HPP
#define TSFDOUBLE_HPP

#include "TSFConfigDefs.hpp"
#include "TSFVector.hpp"
#include "TSFEpetraVector.hpp"

namespace TSFDouble
{
  typedef TSFExtended::Vector<double> Vector;
  typedef TSFExtended::EpetraVector EpetraVector;
}

#endif
