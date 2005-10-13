#include "TSFProductSpace.h"
#include "TSFVectorSpace.h"
#include "TSFBlockVector.h"
#include "TSFVectorBase.h"


using namespace TSF;
using namespace std;

TSFProductSpace::TSFProductSpace(const TSFVectorSpace& space)
	: blocks_(1)
{
	blocks_[0] = space;
}

TSFProductSpace::TSFProductSpace(const TSFVectorSpace& space0,
																 const TSFVectorSpace& space1)
	: blocks_(2)
{
	blocks_[0] = space0;
	blocks_[1] = space1;
}

TSFProductSpace::TSFProductSpace(const TSFVectorSpace& space0,
																 const TSFVectorSpace& space1,
																 const TSFVectorSpace& space2)
	: blocks_(3)
{
	blocks_[0] = space0;
	blocks_[1] = space1;
	blocks_[2] = space2;
}

TSFProductSpace::TSFProductSpace(const TSFVectorSpace& space0,
																 const TSFVectorSpace& space1,
																 const TSFVectorSpace& space2,
																 const TSFVectorSpace& space3)
	: blocks_(4)
{
	blocks_[0] = space0;
	blocks_[1] = space1;
	blocks_[2] = space2;
	blocks_[3] = space3;
}

TSFProductSpace::TSFProductSpace(const TSFArray<TSFVectorSpace>& blocks)
	: blocks_(blocks)
{}



TSFVectorSpaceBase* TSFProductSpace::deepCopy() const 
{
	TSFArray<TSFVectorSpace> copy(blocks_.size());
	for (int i=0; i<copy.size(); i++)
		{
			copy[i] = blocks_[i].deepCopy();
		}
	return new TSFProductSpace(copy);
}

int TSFProductSpace::dim() const 
{
	int rtn = 0;

	for (int i=0; i<numBlocks(); i++)
		{
			rtn += blocks_[i].dim();
		}
	
	return rtn;
}

TSFVectorBase* TSFProductSpace::createMember(const TSFVectorSpace& handle) const 
{
	return new TSFBlockVector(handle);
}

bool TSFProductSpace::checkEquality(const TSFVectorSpaceBase* other) const
{
	/* the type of the argument has been checked at the handle level, so
	 * we can safely cast to a TSFProductSpace */
	const TSFProductSpace* otherGuy 
		= dynamic_cast<const TSFProductSpace*>(other);

	/* first check is that they have to have the same number of blocks */
	if (numBlocks() != otherGuy->numBlocks()) return false;

	/* now compare block-by-block. If any pair of blocks is mismatched, the
	 * spaces are not equal */
	for (int i=0; i<numBlocks(); i++)
		{
			if (blocks_[i] != otherGuy->blocks_[i]) return false;
		}

	/* if we've gotten this far, the spaces are equal */
	return true;
}

int TSFProductSpace::numBlocks() const 
{
	return blocks_.size();
}

void TSFProductSpace::getBlock(int i, const TSFVectorSpace& /* self */,
															 TSFVectorSpace& sub) const
{
	sub = blocks_[i];
}

void TSFProductSpace::print(ostream& os) const 
{
  os << "ProductSpace";
}
