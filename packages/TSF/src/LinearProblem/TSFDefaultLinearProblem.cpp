#include "TSFDefaultLinearProblem.h"


TSFDefaultLinearProblem::TSFDefaultLinearProblem()
	: TSFLinearProblemBase()
{}

TSFVector TSFDefaultLinearProblem::getRHS(const TSFVectorSpace& space) const 
{
	TSFVector rhs = space.createMember();

	int low = space.firstLocalIndex();
	int high = space.nLocalIndices();

	for (int i=low; i<high; i++)
		{
			rhs.setElement(i, getRHSRowValue(i, space));
		}
	return rhs;
}

TSFVector TSFDefaultLinearProblem::getKnownSoln() const 
{
	TSFVector soln = space.createMember();

	int low = space.firstLocalIndex();
	int high = space.nLocalIndices();

	for (int i=low; i<high; i++)
		{
			soln.setElement(i, getSolutionRowValue(i, space));
		}
	return soln;
}




