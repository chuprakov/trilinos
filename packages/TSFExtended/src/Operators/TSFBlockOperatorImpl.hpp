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


#ifndef TSFBLOCKOPERATORIMPL_HPP
#define TSFBLOCKOPERATORIMPL_HPP


#include "TSFProductVectorSpaceDecl.hpp"
#include "TSFLinearCombination.hpp"
#include "TSFZeroOperator.hpp"

using namespace TSFExtended;
using namespace TSFExtendedOps;
using namespace Teuchos;
using std::ostream;





/*==================================================================*/
template <class Scalar>
BlockOperator<Scalar>::BlockOperator(const VectorSpace<Scalar>& domain,
				     const  VectorSpace<Scalar>& range)
  :isFinal_(false),
   domain_(domain),
   range_(range),
   nBlockRows_(range.numBlocks()), //put in check for prod vect sp
   
   nBlockCols_(domain.numBlocks()),
   myBlockRow_(5),
   sub_(range.numBlocks()),
   isSet_(range.numBlocks())
{
  for (int i=0; i<nBlockRows_; i++)
    {
      sub_[i].resize(nBlockCols_);
      isSet_[i].resize(nBlockCols_);
      for (int j = 0; j < nBlockCols_; j++)
	{
	  isSet_[i][j] = false;
	}
    }
}



/*==================================================================*/
template <class Scalar>
RefCountPtr<const TSFCore::VectorSpace<Scalar> >  
BlockOperator<Scalar>::domain() const
{
  return domain_.ptr();
}

/*==================================================================*/
template <class Scalar>
RefCountPtr<const TSFCore::VectorSpace<Scalar> >  
BlockOperator<Scalar>::range() const
{
  return range_.ptr();
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

  TEST_FOR_EXCEPTION(!isSet_[i][j], runtime_error,
		     "Block (" << i << ", " << j << ") is not set." << endl);

  return sub_[i][j];
}




/*==================================================================*/
template <class Scalar>
void BlockOperator<Scalar>::setBlock( int i,  int j, 
				      const LinearOperator<Scalar>& sub)
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
  isFinal_ = chkFinal();
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
  if (M_trans == TSFCore::NOTRANS)
    {
      applyReg(x, alpha, y, beta);  
    }
  else
    {
      applyTrans(x, alpha, y, beta);
    }
      
}





/*==================================================================*/
template <class Scalar>
void BlockOperator<Scalar>::applyReg(const TSFCore::Vector<Scalar>& arg, 
				     const Scalar alpha,
				     TSFCore::Vector<Scalar>* out, 
				     const Scalar beta) const
{
  TEST_FOR_EXCEPTION(!isFinal_, runtime_error, "Operator not finalized");

  for (int i=0; i<nBlockRows_; i++)
    {
      Vector<Scalar> tmpE = range_.getBlock(i).createMember();
      TSFCore::Vector<Scalar>* tmp = tmpE.ptr().get();

      TSFCore::assign(tmp, 0.0);
      const TSFCore::ProductVector<Scalar>* argPV = 
	dynamic_cast<const TSFCore::ProductVector<Scalar>* > (&arg);
      TSFCore::ProductVector<Scalar>* outPV = 
	dynamic_cast< TSFCore::ProductVector<Scalar>* > (out);
      for (int j=0; j<nBlockCols_; j++)
	{
	  (sub_[i][j].ptr())->apply(TSFCore::NOTRANS, 
				    *(argPV->getBlock(j).get()), 
				    tmp, 1.0, 1.0);
	}
      if (beta == 0.0)
	{
	  TSFCore::Vt_S(tmp, alpha);
	}
      else
	{
	  TSFCore::Vt_S(tmp, alpha);
	  TSFCore::Vp_StV(tmp, beta, *(outPV->getBlock(i).get()));
	}
      Teuchos::RefCountPtr<TSFCore::Vector<Scalar> > block = outPV->getBlock(i);
      TSFCore::assign(block.get(), *tmp);
    }
}

/*==================================================================*/
template <class Scalar>
void BlockOperator<Scalar>::applyTrans(const TSFCore::Vector<Scalar>& arg,
				       const Scalar alpha,
				       TSFCore::Vector<Scalar>* out,
				       const Scalar beta) const
{
  TEST_FOR_EXCEPTION(!isFinal_, runtime_error, "Operator not finalized");

  for (int i=0; i<nBlockCols_; i++)
    {
      Vector<Scalar> tmpE = domain_.getBlock(i).createMember();
      TSFCore::Vector<Scalar>* tmp = tmpE.ptr().get();

      TSFCore::assign(tmp, 0.0);
      const TSFCore::ProductVector<Scalar>* argPV = 
	dynamic_cast<const TSFCore::ProductVector<Scalar>* > (&arg);
      TSFCore::ProductVector<Scalar>* outPV = 
	dynamic_cast<TSFCore::ProductVector<Scalar>* > (out);
      for (int j=0; j<nBlockRows_; j++)
	{
	  (sub_[j][i].ptr())->apply(TSFCore::TRANS, 
				    *(argPV->getBlock(j).get()), 
				    tmp, 1.0, 1.0);
	}
      if (beta == 0.0)
	{
	  TSFCore::Vt_S(tmp, alpha);
	}
      else
	{
	  TSFCore::Vt_S(tmp, alpha);
	  TSFCore::Vp_StV(tmp, beta, *(outPV->getBlock(i).get()));
	}
      Teuchos::RefCountPtr<TSFCore::Vector<Scalar> > block = outPV->getBlock(i);
      TSFCore::assign(block.get(), *tmp);
    }
}





/*==================================================================*/
template <class Scalar>
LinearOperator<Scalar> BlockOperator<Scalar>::formTranspose() const
{
  TEST_FOR_EXCEPTION(!isFinal_, runtime_error, "Operator not finalized");

  BlockOperator<Scalar> trans(range_, domain_);
  
  for (int i = 0; i < numBlockRows(); i++)
    {
      for (int j = 0; j < numBlockCols(); j++)
	{
	  LinearOperator<Scalar> B = sub_[i][j];
	  LinearOperator<Scalar> Btrp = (B.transpose()).form(); 
	  trans.setBlock(j, i, Btrp);
	}
    }
  return LinearOperator<Scalar>(&trans);
}





/*==================================================================*/
template <class Scalar>
void BlockOperator<Scalar>::getRow(const int& row, 
				   Teuchos::Array<int>& indices,
				   Teuchos::Array<Scalar>& values) const
{
  /* find col block in which row exists and row within the block */
  int K = 0;
  int blockRow = 0;
  int rowInBlock = 0;
  for (int i = 0; i < nBlockRows_; i++)
    {
      int numR = getBlock(i, 0).range().dim();
      K += numR;
      if (row < K)
	{
	  blockRow = i;
	  rowInBlock = row - (K - numR);
	  break;
	}
    }

  /* get the row elements for each block in the row.  */
  int offset = 0;
  Teuchos::Array<int> localInd;
  Teuchos::Array<Scalar> localVal;
  for (int i = 0; i < nBlockCols_; i++)
    {
      sub_[blockRow][i].getRow(rowInBlock,localInd, localVal);
      for (int j = 0; j < localInd.size(); j++)
	{
	  indices.append(localInd[j] + offset);
	  values.append(localVal[j]);
	}
      offset += sub_[blockRow][i].domain().dim();
    }

}




/*==================================================================*/
template <class Scalar>
void BlockOperator<Scalar>::print(ostream& os) const 
{
  TEST_FOR_EXCEPTION(!isFinal_, runtime_error, "Operator not finalized");
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
  isFinal_ = true;
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
	  if (!isSet_[i][j]) return false;
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
string BlockOperator<Scalar>::describe(int depth) const
{
  string spaces = "";
  for (int i = 0; i < depth; i++)
    {
      spaces.append("   ");
    }

  string ret = "";
  ret.append(spaces + "Block Operator with "
	     + Teuchos::toString(nBlockRows_)
	     + " rows and  "  + Teuchos::toString(nBlockCols_)
	     + " columns\n");
  if (!isFinal_) 
    {
      ret.append(spaces + "    Operator not finalized");
      return ret;
    }

  for (int i = 0; i < nBlockRows_; i++)
    {
      for (int j = 0; j < nBlockCols_; j++)
	{
	  ret.append(spaces + "   Block " + Teuchos::toString(i) + ", "
		     + Teuchos::toString(j) + " is");
	  OpDescribableByTypeID<Scalar>* mat = 
	    dynamic_cast<OpDescribableByTypeID<Scalar>* > (sub_[i][j].ptr().get());
	  if (mat == 0)
	    {
	      ret.append(spaces + "   Unknown type of dimension " 
			 + Teuchos::toString(getBlock(i,j).range().dim()) 
			 + Teuchos::toString(getBlock(i,j).domain().dim()) + "\n");
	    }
	  else
	    {
	      ret.append(mat->describe(depth+1) + "\n");
	    }

	}
    }
  return ret; 
}



#endif
