// /////////////////////////////////////////////////////////////////////////
// TSFCoreTestingTools.hpp

#ifndef TSFCORE_TESTING_TOOLS_HPP
#define TSFCORE_TESTING_TOOLS_HPP

#include "TSFCoreTestingToolsDecl.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "TSFCoreLinearOp.hpp"

template<class Scalar>
std::ostream& TSFCore::operator<<( std::ostream& o, const Vector<Scalar>& v )
{
	RTOpPack::SubVectorT<Scalar> sv;
	v.getSubVector(Range1D(),&sv);
	for( Index i = 1; i <= sv.subDim(); ++i )
		o << " " << sv(i);
	o << std::endl;
	v.freeSubVector(&sv);
	return o;
}

template<class Scalar>
std::ostream& TSFCore::operator<<( std::ostream& o, const LinearOp<Scalar>& M )
{
	namespace mmp = MemMngPack;
	const Index dimDomain = M.domain()->dim(), dimRange = M.range()->dim();
	// We will extract by column if op==NOTRANS is supported and by row otherwise
	const ETransp opM = M.opSupported(NOTRANS) ? NOTRANS : TRANS;
	// Copy into dense matrix (by column or row)
	Teuchos::RefCountPtr<Vector<Scalar> >
		e_j = ( opM==NOTRANS ? M.domain() : M.range()  )->createMember(),
		t   = ( opM==NOTRANS ? M.range()  : M.domain() )->createMember(); // temp column or row
	const Index
		dimOpMDomain = ( opM==NOTRANS ? dimDomain : dimRange  ),
		dimOpMRange  = ( opM==NOTRANS ? dimRange  : dimDomain );
	RTOpPack::SubVectorT<Scalar> sv;
	std::vector<Scalar>  Md( dimOpMRange * dimOpMDomain ); // Column major
	const Index
		cs = ( opM==NOTRANS ? 1         : dimRange ),  // stride for columns or rows 
		rs = ( opM==NOTRANS ? dimRange  : 1        );  // stride for rows or columns
	Index i, j;
	for( j = 1; j <= dimOpMDomain; ++j ) {
		TSFCore::assign( e_j.get(), 0.0 );
		TSFCore::set_ele( j, 1.0, e_j.get() );
		M.apply(opM,*e_j,t.get());  // extract the ith column or row
		t->getSubVector(Range1D(),&sv);
		for( i = 1; i <= dimOpMRange; ++i ) Md[ (i-1)*cs + (j-1)*rs ] = sv(i);
		t->freeSubVector(&sv);
	}
	// Print the matrix
	for( i = 1; i <= dimRange; ++i ) {
		for( j = 1; j <= dimDomain; ++j )
			o << " " << Md[ (i-1) + (j-1)*dimRange ];
		o << std::endl;
	}
	return o;
}

#endif // TSFCORE_TESTING_TOOLS_HPP
