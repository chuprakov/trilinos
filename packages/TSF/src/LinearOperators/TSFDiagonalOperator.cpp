#include <iostream>

#include "TSFDiagonalOperator.h"



using namespace TSF;



TSFDiagonalOperator::TSFDiagonalOperator(const TSFVector& diagonalValues)
  :
  TSFLinearOperatorBase(diagonalValues.space(), diagonalValues.space()),
  diagonalValues_(diagonalValues)
{
  ;
}


void TSFDiagonalOperator::apply(const TSFVector& in, TSFVector& out) const 
{
  out.dotStar(in, diagonalValues_);
}


void TSFDiagonalOperator::applyAdjoint(const TSFVector& in, 
                                       TSFVector& out) const 
{
  out.dotStar(in, diagonalValues_);
}

void TSFDiagonalOperator::applyInverse(const TSFVector& in, 
                                       TSFVector& out) const 
{
  out.dotSlash(in, diagonalValues_);
}


void TSFDiagonalOperator::print(ostream& os) const 
{
  diagonalValues_.print(os);
}



void TSFDiagonalOperator::getRow(int row, TSFArray<int>& indices, 
                                 TSFArray<TSFReal>& values) const
{
  indices.append(row);
  values.append(diagonalValues_[row]);
}
