#include "TSFError.h"
#include "TSFLinearOperator.h"
#include "TSFBlockLinearOperator.h"
#include "TSFZeroOperator.h"
#include "TSFProductSpace.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"

#include <iostream.h>

using namespace TSF;
using std::ostream;


TSFBlockLinearOperator::TSFBlockLinearOperator(const TSFVectorSpace& domain,
																							 const TSFVectorSpace& range)
	: TSFLinearOperatorBase(domain, range),
		nBlockRows_(range.numBlocks()), 
		nBlockCols_(domain.numBlocks()),
		sub_(range.numBlocks())
{
	for (int i=0; i<nBlockRows_; i++)
		{
			sub_[i].resize(nBlockCols_);
			TSFVectorSpace r = range.getBlock(i);
			for (int j=0; j<nBlockCols_; j++)
				{
					TSFVectorSpace d = domain.getBlock(j);
					setBlock(i, j, new TSFZeroOperator(d, r));
				}
		}
}

void TSFBlockLinearOperator::getBlock(int i, int j, 
																			TSFLinearOperator& sub) const 
{
	if (i < 0 || i >= nBlockRows_ || j<0 || j>=nBlockCols_)
		{
			TSFError::raise("TSFBlockLinearOperator::getBlock index out of range");
		}
	sub = sub_[i][j];
}

void TSFBlockLinearOperator::setBlock(int i, int j, 
																			const TSFLinearOperator& sub)
{
	if (i < 0 || i >= nBlockRows_ || j<0 || j>=nBlockCols_)
		{
			TSFError::raise("TSFBlockLinearOperator::setBlock index out of range");
		}
/*     if (sub.domain() != domain_.getBlock(j)) */
/*       { */
/*         cerr << "error in BlockLinearOp: setBlcok: domain\n"; */
/*       } */
/*     if (sub.range() != range_.getBlock(i)) */
/*       { */
/*         cerr << "error in BlockLinearOp: setBlcok: range\n"; */
/*       } */
	sub_[i][j] = sub;
}



void TSFBlockLinearOperator::apply(const TSFVector& arg,
																	 TSFVector& out) const
{
	for (int i=0; i<nBlockRows_; i++)
		{
			TSFVector tmpRow = range().getBlock(i).createMember();
            tmpRow.zero();
			for (int j=0; j<nBlockCols_; j++)
				{
					TSFVector tmp = range().getBlock(i).createMember();
					sub_[i][j].apply(arg.getBlock(j), tmp);
                    tmpRow.add(tmp, tmpRow);
				}
			out.setBlock(i, tmpRow);
		}
}


void TSFBlockLinearOperator::applyAdjoint(const TSFVector& arg,
																					TSFVector& out) const
{
	for (int i=0; i<nBlockCols_; i++)
		{
			TSFVector tmpRow = domain().getBlock(i).createMember();
			for (int j=0; j<nBlockRows_; j++)
				{
					TSFVector tmp = domain().getBlock(i).createMember();
					sub_[j][i].applyAdjoint(arg.getBlock(j), tmp);
					tmpRow.add(tmp, tmpRow);
				}
			out.setBlock(i, tmpRow);
		}
}


TSFLinearOperator* TSFBlockLinearOperator::getTranspose()
{
  opTrp_ = new TSFBlockLinearOperator(range(), domain());
  for (int i = 0; i < numBlockRows(); i++)
    {
      for (int j = 0; j < numBlockCols(); j++)
        {
          TSFLinearOperator B = sub_[i][j];
          TSFLinearOperator* Btrp = &(B.getTranspose());
          opTrp_.setBlock(j, i, *Btrp);
        }
    }
  return &opTrp_;
}




void TSFBlockLinearOperator::print(ostream& os) const 
{
	os << "<BlockLinearOperator nRows=\"" << nBlockRows_ 
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
	os << "</BlockLinearOperator>" << endl;
}
