/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// **********************************************************************/
/* @HEADER@ */

#include "TSFBlockOperator.hpp"
// #include "TSFProductVectorSpace.hpp"
// #include "TSFVectorSpace.hpp"
// //#include "TSFCoreVectorSpace.hpp"
// #include "TSFOpDescribableByTypeID.hpp"
 
#include "TSFZeroOperator.hpp"

using namespace TSFExtended;
using namespace Teuchos;
using std::ostream;


/*==================================================================*/
template <class Scalar>
BlockOperator<Scalar>::BlockOperator()
  : isUnspecified_(true),
    isFinal_(false)
{
  domain_ = new ProductVectorSpace<Scalar>();
  range_ = new ProductVectorSpace<Scalar>();
  sub_.resize(0);
  isSet_.resize(0);
}



/*==================================================================*/
template <class Scalar>
BlockOperator<Scalar>::BlockOperator(const VectorSpace<Scalar>& domain,
				     const VectorSpace<Scalar>& range)
  :isUnspecified_(false),
   isFinal_(false),
   domain_(domain),
   range_(range),
   nBlockRows_(range.numBlocks()), //put in check for prod vect sp
   nBlockCols_(domain.numBlocks()),
   sub_(range.numBlocks())
{
  for (int i=0; i<nBlockRows_; i++)
    {
      sub_[i].resize(nBlockCols_);
      isSet_[i].resize(nBlockCols_);
    }
}







/*==================================================================*/
template <class Scalar>
LinearOperator<Scalar> BlockOperator<Scalar>::getBlock(const int &i, 
						       const int &j) const
{
  TEST_FOR_EXCEPTION(i < 0 || i >=nBlockRows_, std::out_of_range,
		     "i is out of range in setBlock: i = " << i 
		     << "and j = " << j << endl);
  TEST_FOR_EXCEPTION(j < 0 || j >=nBlockCols_, std::out_of_range,
		     "j is out of range in setBlock: i = " << i 
		     << "and j = " << j << endl);

  TEST_FOR_EXCEPTION(!isSet[i][j], runtime_error,
		     "Block (" << i << ", " << j << ") is not set." << endl);

  return sub_[i][j];
}




/*==================================================================*/
template <class Scalar>
void BlockOperator<Scalar>::setBlock(const int &i, const int &j, 
				     const LinearOperator<Scalar>& sub)
{
  if (isUnspecified_)
    {
      buildSpaces(i, j, sub);   
      sub_[i][j] = sub;
      isSet_[i][j] = true;
    }
  else
    {
      TEST_FOR_EXCEPTION(i < 0 || i >=nBlockRows_, std::out_of_range,
			 "i is out of range in setBlock: i = " << i 
			 << "and j = " << j << endl);
      TEST_FOR_EXCEPTION(j < 0 || j >=nBlockCols_, std::out_of_range,
			 "j is out of range in setBlock: i = " << i 
			 << "and j = " << j << endl);

      chkSpaces(i, j, sub);
      sub_[i][j] = sub;
      isSet_[i][j] = true;
      isFinal_ = chkfinal();
    }
}  


/*==================================================================*/
template <class Scalar>
void BlockOperator<Scalar>::apply(
				  const TSFCore::ETransp            M_trans
				  ,const TSFCore::Vector<Scalar>    &x
				  ,TSFCore::Vector<Scalar>          *y
				  ,const Scalar            alpha = 1.0
				  ,const Scalar            beta  = 0.0
				  ) const 
{
//   TEST_FOR_EXCEPTION(dynamic_cast<ZeroOperator<Scalar>* >(op_) != 0, runtime_error,
// 		     "InverseOperator<Scalar>::apply() called on a ZeroOperator.");
//   TEST_FOR_EXCEPTION(op_->domain().dim() != op_->range().dim(), runtime_error,
// 		     "InverseOperator<Scalar>::apply() called on a non-square operator.");
//   SolverState<Scalar> haveSoln;
//   LinearOperator<Scalar> applyOp;      
//   if (M_trans == NOTRANS)
//     {
//       applyOp = op_;
//     }
//   else
//     {
//       applyOp = op_.transpose();
//     }
      
//   if (beta == 0.0)
//     {
//       Vt_S(&x, alpha);
//     }
//   else
//     {
//       applyOp.ptr()->apply(NOTRANS, *y, &x, beta, alpha);
//     }
      
//   DiagonalOperator<Scalar>* val1 = dynamic_cast<DiagonalOperator<Scalar>* >(op_);
//   IdentityOperator<Scalar>* val2 = dynamic_cast<IdentityOperator<Scalar>* >(op_);
//   if (val1 != 0 or val2 != 0)
//     {
//       haveSoln = applyOp_.ptr()->applyInverse(x, y);
//     }      
//   else
//     {
//       haveSoln = solver_.solve(applyOp, x, *y);
//     }
//   TEST_FOR_EXCEPTION(haveSoln.finalState() != SolveConverged, runtime_error,
// 		     "InverseOperator<Scalar>::apply() " << haveSoln.stateDescription());
}



// /*==================================================================*/
// template <class Scalar>
// void BlockOperator<Scalar>::apply(const Vector<Scalar>& arg,
// 				  Vector<Scalar>& out) const
// {
//   TEST_FOR_EXCEPTION(isFinal_, runtime_error, "Operator not finalized");
//   for (int i=0; i<nBlockRows_; i++)
//     {
//       Vector<Scalar> tmpRow = range().getBlock(i).createMember();
//       tmpRow.zero();
//       for (int j=0; j<nBlockCols_; j++)
// 	{
// 	  tmpRow = tmpRow + sub_[i][j] * arg.getBlock(j);
// 	}
//       out.setBlock(i, tmpRow);
//     }
// }


// /*==================================================================*/
// template <class Scalar>
// void BlockOperator<Scalar>::applyAdjoint(const Vector<Scalar>& arg,
// 					 Vector<Scalar>& out) const
// {
//   TEST_FOR_EXCEPTION(isFinal_, runtime_error, "Operator not finalized");
//   for (int i=0; i<nBlockCols_; i++)
//     {
//       Vector<Scalar> tmpRow = domain().getBlock(i).createMember();
//       for (int j=0; j<nBlockRows_; j++)
// 	{
// 	  Vector<Scalar> tmp = domain().getBlock(i).createMember();
// 	  sub_[j][i].applyAdjoint(arg.getBlock(j), tmp);
// 	  tmpRow.add(tmp, tmpRow);
// 	}
//       out.setBlock(i, tmpRow);
//     }
// }

/*==================================================================*/
template <class Scalar>
LinearOperator<Scalar>* BlockOperator<Scalar>::formTranspose()
{
  TEST_FOR_EXCEPTION(isFinal_, runtime_error, "Operator not finalized");
  opTrp_ = new BlockOperator(range(), domain());
  for (int i = 0; i < numBlockRows(); i++)
    {
      for (int j = 0; j < numBlockCols(); j++)
	{
	  //cerr << "getting Block i,j = " << i << " " << j << endl;
	  LinearOperator<Scalar> B = sub_[i][j];
	  LinearOperator<Scalar>* Btrp = &(B.getTranspose());
	  opTrp_.setBlock(j, i, *Btrp);
	}
    }
  return &opTrp_;
}




/*==================================================================*/
template <class Scalar>
void BlockOperator<Scalar>::print(ostream& os) const 
{
  TEST_FOR_EXCEPTION(isFinal_, runtime_error, "Operator not finalized");
  os << "<BlockOperator nRows=\"" << nBlockRows_ 
     << "\" nCols=\"" << nBlockCols_ << "\">" << endl;
  for (int i=0; i<nBlockRows_; i++)
    {
      os << "<BlockRow i=\"" << i << "\">" << endl;
      for (int j=0; j<nBlockCols_; j++)
	{
	  os << "<BlockCol j=\"" << j << "\">" << endl;
	  os << sub_[i][j] << endl;
	  os << "</BlockCol>" << endl;
	}
      os << "</BlockRow>" << endl;
    }
  os << "</BlockOperator>" << endl;
}




/*==================================================================*/
template <class Scalar>
void BlockOperator<Scalar>::finalize(const bool &zerofill)
{
  try
    {
      domain_.finalize();
    }
  catch(runtime_error)
    {
      TEST_FOR_EXCEPTION(true, runtime_error, "Domain not complete" << 
			 endl);
    }
  try
    {
      range_.finalize();
    }
  catch(runtime_error)
    {
      TEST_FOR_EXCEPTION(true, runtime_error, "Range not complete" << 
			 endl);
    }
  for (int i = 0; i < nBlockRows_; i++)
    {
      for (int j = 0; j < nBlockCols_; j++)
	{
	  if (!isSet_[i][j]) 
	    {
	      if (zerofill) 
		{
		  zeroFill(i, j);
		}
	      else
		{
		  TEST_FOR_EXCEPTION(true, runtime_error,
				     "Block (" << i << ", " << j << 
				     ") not set and zerofill is false" 
				     << endl);
		}
	    }
	}
    }
}





/*==================================================================*/
template <class Scalar>
bool BlockOperator<Scalar>::buildSpaces(const int &i, const int &j, 
					const LinearOperator<Scalar>& sub)
{
  try
    {
      domain_.setBlock(j, sub.domain());  //errors here caught in ProdSpace
    }
  catch(runtime_error)
    {
      TEST_FOR_EXCEPTION(true, runtime_error,
			 "Domain not compatible while setting block (" << i
			 << ", " << j << ")" << endl);
    }
  numCols_ = nBlockCols_;
  nBlockCols_ = domain_.getNumBlocks();

  try
    {
      range_.setBlock(i, sub.range());
    }
  catch(runtime_error)
    {
      TEST_FOR_EXCEPTION(true, runtime_error,
			 "Range not compatible while setting block (" << i
			 << ", " << j << ")" << endl);
    }
  numRows = nBlockRows_;
  nBlockRows_ = range_.getNumBlocks();

  sub_.resize(nBlockRows_);
  isSet_.resize(nBlockRows_);
  for (int i = 0; i < nBlockRows_; i++)
    {
      sub_[i].resize(nBlockRows_);
      isSet_[i].resize(nBlockRows_);
    }

  // Set new rows of isSet_ to false
  for (int i = nBlockRows_ - 1; i >= numRows; i--)
    {
      for (int j = 0; j < nBlockCols_; j++)
	{
	  isSet_[i][j] = false;
	}
    }
  //Set new cols of isSet_ to false
  for (int j = nBlockCols_ -1; j >= numCols; j--)
    {
      for (int i = 0; i < nBlockRows_; i++)
	{
	  isSet[i][j] = false;
	}
    }
}


/*==================================================================*/
template <class Scalar>
void BlockOperator<Scalar>::chkSpaces(const int &i, const int &j, 
				      const LinearOperator<Scalar>& sub) const
{
  if (sub.domain() != domain_.getBlock(j))
    {
      TEST_FOR_EXCEPTION(true, runtime_error, 
			 "Domain not compatible in block (" << i 
			 << ", " << j << ")" << endl);
    }
  if (sub.range() != range_.getBlock(i))
    {
      TEST_FOR_EXCEPTION(true, runtime_error, 
			 "Range not compatible in block (" << i 
			 << ", " << j << ")" << endl);
    }
}


/*==================================================================*/
template <class Scalar>
bool BlockOperator<Scalar>:: chkFinal() const
{
  for (int i = 0; i < nBlockCols_; i++)
    {
      for (int j = 0; j < nBlockRows_; j++)
	{
	  if (!isSet[i][j]) return false;
	}
    }
  return true;
}



/*==================================================================*/
template <class Scalar>
void BlockOperator<Scalar>:: zeroFill(const int &i, const int &j)
{
  VectorSpace<Scalar> d = domain_.getBlock(j);
  VectorSpace<Scalar> r = range_.getBlock(i);
  LinearOperator<Scalar> z = new ZeroOperator<Scalar>(d, r);
  setBlock(i, j, z);
}



/*==================================================================*/
template <class Scalar>
string BlockOperator<Scalar>::describe(const int &depth) const
{
  TEST_FOR_EXCEPTION(isFinal_, runtime_error, "Operator not finalized");
  string spaces = "";
  for (int i = 0; i < depth; i++)
    {
      spaces.append("   ");
    }

  string ret = "";
  ret.append(spaces + "Block Operator with "
	     + toString(numBlockRows())
	     + " rows and  "  + toString(numBlockCols())
	     + " columns" + endl);
  for (int i = 0; i < numBlockRows(); i++)
    {
      for (int j = 0; j < numBlockCols(); j++)
	{
	  ret.append(spaces + "   Block " + toString(i) + ", "
		     + toString(j) + " is" + endl);
	  OpDescribableByTypeID<Scalar> mat = 
	    dynamic_cast<OpDescribableByTypeID<Scalar> > 
	    (&(*(getBlock(i,j).getPtr())));
	  if (mat == 0)
	    {
	      ret.append(spaces + "   Unknown type of dimension " 
			 + toString(getBlock(i,j).range().dim()) 
			 + toString(getBlock(i,j).domain()) + endl);
	    }
	  else
	    {
	      ret.append(mat.describe(depth+1));
	    }

	}
    }
}





// /*==================================================================*/
// template <class Scalar>
// void BlockOperator<Scalar>::getBlock(int i, int j, 
// 			     LinearOperator<Scalar>& sub) const 
// {
//   TEST_FOR_EXCEPTION(i < 0 || i >=nBlockRows_, std::out_of_range,
// 		     "i is out of range in setBlock: i = " << i 
// 		     << "and j = " << j << endl);
//   TEST_FOR_EXCEPTION(j < 0 || j >=nBlockCols_, std::out_of_range,
// 		     "j is out of range in setBlock: i = " << i 
// 		     << "and j = " << j << endl);

//   TEST_FOR_EXCEPTION(!isSet[i][j], runtime_error,
// 		     "Block (" << i << ", " << j << ") is not set." << endl);

//   sub = sub_[i][j];
// }
