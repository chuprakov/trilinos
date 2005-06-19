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
// ***********************************************************************
// @HEADER

#ifndef TSFCORE_VECTOR_STD_OPS_TESTER_HPP
#define TSFCORE_VECTOR_STD_OPS_TESTER_HPP

#include "TSFCoreVectorStdOpsTesterDecl.hpp"
#include "TSFCoreTestingTools.hpp"
#include "Teuchos_arrayArg.hpp"

namespace TSFCore {

// VectorStdOpsTesterComparable (using partial specialization to only do tests in some cases)

template <bool isComparable, class Scalar>
class VectorStdOpsTesterComparable {
public:
	static bool checkComparableStdOps(
		const VectorSpace<Scalar>                                     &vecSpc
		,Vector<Scalar>                                               *z
		,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &maxRelErr
		,std::ostream                                                 *out
		,const bool                                                   &dumpAll
		)
		{
			//TEST_FOR_EXCEPTION(true,std::logic_error,"Error, This should not even compile?");
			//return false;
			return Teuchos::ScalarTraits<Scalar>::ThisShouldNotCompile();
		}
};

template <class Scalar>
class VectorStdOpsTesterComparable<false,Scalar> {
public:
	static bool checkComparableStdOps(
		const VectorSpace<Scalar>                                     &vecSpc
		,Vector<Scalar>                                               *z
		,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &maxRelErr
		,std::ostream                                                 *out
		,const bool                                                   &dumpAll
		)
		{
			if(out) *out << "\nThis scalar type does not support comparable operations so we can not test min(), max() and other such functions.\n";
			return true;
		}
};

template <class Scalar>
class VectorStdOpsTesterComparable<true,Scalar> {
public:
	static bool checkComparableStdOps(
		const VectorSpace<Scalar>                                     &vecSpc
		,Vector<Scalar>                                               *z
		,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &maxRelErr
		,std::ostream                                                 *out
		,const bool                                                   &dumpAll
		)
		{
			typedef Teuchos::ScalarTraits<Scalar> ST;

			bool success = true, result;
			
			if(out) *out << "\nTesting comparable operations ...\n";
			
			const Scalar scalarSmall(1e-5), scalarMedium(2.0), scalarLarge(100.0);
			if(out) *out << "\nassign(&*z,"<<scalarMedium<<");\n";
			assign(&*z,Scalar(scalarMedium));
			if(out && dumpAll) *out << "\nz =\n" << *z;
			if(out) *out << "\nset_ele(1,"<<scalarSmall<<",&*z);\n";
			set_ele(1,scalarSmall,&*z);
			if(out && dumpAll) *out << "\nz =\n" << *z;
			if(out) *out << "\nset_ele(2,"<<scalarLarge<<",&*z);\n";
			set_ele(2,scalarLarge,&*z);
			if(out && dumpAll) *out << "\nz =\n" << *z;
			if(out) *out << "\nset_ele(vecSpc.dim()-1,"<<scalarSmall<<",&*z);\n";
			set_ele(vecSpc.dim()-1,scalarSmall,&*z);
			if(out && dumpAll) *out << "\nz =\n" << *z;
			if(out) *out << "\nset_ele(vecSpc.dim(),"<<scalarLarge<<",&*z);\n";
			set_ele(vecSpc.dim(),scalarLarge,&*z);
			if(out && dumpAll) *out << "\nz =\n" << *z;

			Scalar minEle; Index minIndex;
			Scalar maxEle; Index maxIndex;

			if(!testRelErr<Scalar>(
					 "min(*z)",min(*z)
					 ,"scalarSmall",scalarSmall
					 ,"maxRelErr",maxRelErr,out
					 )
				) success=false;

			if(out) *out << "\nmin(*z,&minEle,&minIndex);\n";
			minEle = ST::zero(); minIndex = 0;
			min(*z,&minEle,&minIndex);
			if(!testRelErr<Scalar>(
					 "minEle",minEle
					 ,"scalarSmall",scalarSmall
					 ,"maxRelErr",maxRelErr,out
					 )
				) success=false;
			result = minIndex == 1;
			if(out) *out << "\nminIndex = " << minIndex << " == 1 ? " << passfail(result) << std::endl;
			if(!result) success = false;

			if(out) *out << "\nminGreaterThanBound(*z,"<<scalarMedium<<",&minEle,&minIndex);\n";
			minEle = ST::zero(); minIndex = 0;
			minGreaterThanBound(*z,scalarMedium,&minEle,&minIndex);
			if(!testRelErr<Scalar>(
					 "minEle",minEle
					 ,"scalarLarge",scalarLarge
					 ,"maxRelErr",maxRelErr,out
					 )
				) success=false;
			result = minIndex == 2;
			if(out) *out << "\nminIndex = " << minIndex << " == 2 ? " << passfail(result) << std::endl;
			if(!result) success = false;

			if(out) *out << "\nminGreaterThanBound(*z,"<<scalarLarge<<",&minEle,&minIndex);\n";
			minEle = ST::zero(); minIndex = 0;
			minGreaterThanBound(*z,scalarLarge,&minEle,&minIndex);
			result = minIndex < 0;
			if(out) *out << "\nminIndex = " << minIndex << " < 0 ? " << passfail(result) << std::endl;
			if(!result) success = false;
		
			if(!testRelErr<Scalar>(
					 "max(*z)",max(*z)
					 ,"scalarLarge",scalarLarge
					 ,"maxRelErr",maxRelErr,out
					 )
				) success=false;

			if(out) *out << "\nmax(*z,&maxEle,&maxIndex);\n";
			maxEle = ST::zero(); maxIndex = 0;
			max(*z,&maxEle,&maxIndex);
			if(!testRelErr<Scalar>(
					 "maxEle",maxEle
					 ,"scalarLarge",scalarLarge
					 ,"maxRelErr",maxRelErr,out
					 )
				) success=false;
			result = maxIndex == 2;
			if(out) *out << "\nmaxIndex = " << maxIndex << " == 2 ? " << passfail(result) << std::endl;
			if(!result) success = false;

			if(out) *out << "\nmaxLessThanBound(*z,"<<scalarMedium<<",&maxEle,&maxIndex);\n";
			maxEle = ST::zero(); maxIndex = 0;
			maxLessThanBound(*z,scalarMedium,&maxEle,&maxIndex);
			if(!testRelErr<Scalar>(
					 "maxEle",maxEle
					 ,"scalarSmall",scalarSmall
					 ,"maxRelErr",maxRelErr,out
					 )
				) success=false;
			result = maxIndex == 1;
			if(out) *out << "\nmaxIndex = " << maxIndex << " == 1 ? " << passfail(result) << std::endl;
			if(!result) success = false;

			if(out) *out << "\nmaxLessThanBound(*z,"<<scalarSmall<<",&maxEle,&maxIndex);\n";
			maxEle = ST::zero(); maxIndex = 0;
			maxLessThanBound(*z,scalarSmall,&maxEle,&maxIndex);
			result = ( maxIndex < 0 );
			if(out) *out << "\nmaxIndex = " << maxIndex << " < 0 ? " << passfail(result) << std::endl;
			if(!result) success = false;
			
			return success;
		}
};

// VectorStdOpsTester

template <class Scalar>
VectorStdOpsTester<Scalar>::VectorStdOpsTester(
	const ScalarMag   &maxRelErr
	)
	:maxRelErr_(maxRelErr)
{}

template <class Scalar>
bool VectorStdOpsTester<Scalar>::checkStdOps(
	const VectorSpace<Scalar>    &vecSpc
	,std::ostream                *out
	,const bool                  &dumpAll
	)
{
	using Teuchos::arrayArg;
	typedef Teuchos::ScalarTraits<Scalar> ST;

	if(out)
		*out << "\n*** Entering VectorStdOpsTester<"<<ST::name()<<">::checkStdOps(...) ...\n"
				 << "using a \'" << vecSpc.description() << "\' object ...\n";

	bool success = true;
	//bool result;
	//Scalar sresult1, sresult2;
	if(out) *out << "\nvecSpc.dim() = " << vecSpc.dim() << std::endl;

	if(out) *out << "\nCreating vectors v1, v2, v3 and z ...\n";
	Teuchos::RefCountPtr<Vector<Scalar> >
		v1 = vecSpc.createMember(),
		v2 = vecSpc.createMember(),
		v3 = vecSpc.createMember(),
		z  = vecSpc.createMember();

	if(out) *out << "\nassign(&*v1,-2.0);\n";
	assign(&*v1,Scalar(-2.0));
	if(out) *out << "\nassign(&*v2,-3.0);\n";
	assign(&*v2,Scalar(-3.0));
	if(out) *out << "\nassign(&*v3,-4.0);\n";
	assign(&*v3,Scalar(-4.0));
	
	if(out) *out << "\nabs(&*z,*v1);\n";
	abs(&*z,*v1);
	if(!testRelErr<Scalar>("sum(*z)",sum(*z),"2.0*vecSpc.dim()",Scalar(2.0)*Scalar(vecSpc.dim()),"maxRelErr",maxRelErr(),out)) success=false;
	
	if(out) *out << "\nreciprocal(&*z,*v1);\n";
	reciprocal(&*z,*v1);
	if(!testRelErr<Scalar>("sum(*z)",sum(*z),"-0.5*vecSpc.dim()",Scalar(-0.5)*Scalar(vecSpc.dim()),"maxRelErr",maxRelErr(),out)) success=false;

	if(out) *out << "\nlinear_combination(2,{0.5,0.25},{&*v1,&*v2},0.0,&*z);\n";
	linear_combination(2,arrayArg<Scalar>(0.5,0.25)(),arrayArg<const Vector<Scalar>*>(&*v1,&*v2)(),Scalar(0.0),&*z);
	if(!testRelErr<Scalar>("sum(*z)",sum(*z),"(-0.5*2.0-0.25*3.0)*vecSpc.dim()",Scalar(-0.5*2.0-0.25*3.0)*Scalar(vecSpc.dim()),"maxRelErr",maxRelErr(),out)) success=false;

	if(out) *out << "\nassign(&*z,2.0);\n";
	assign(&*z,Scalar(2.0));
	if(out) *out << "\nlinear_combination(3,{0.5,0.25,0.125},{&*v1,&*v2,&*v2},0.5,&*z);\n";
	linear_combination(3,arrayArg<Scalar>(0.5,0.25,0.125)(),arrayArg<const Vector<Scalar>*>(&*v1,&*v2,&*v3)(),Scalar(0.5),&*z);
	if(!testRelErr<Scalar>(
			 "sum(*z)",sum(*z)
			 ,"(0.5*2.0-0.5*2.0-0.25*3.0-0.125*4.0)*vecSpc.dim()",Scalar(0.5*2.0-0.5*2.0-0.25*3.0-0.125*4.0)*Scalar(vecSpc.dim())
			 ,"maxRelErr",maxRelErr(),out
			 )
		) success=false;

	if(out) *out << "\nassign(&*z,2.0);\n";
	assign(&*z,Scalar(2.0));
	if(!testRelErr<Scalar>(
			 "norm_2(*z,*v2)",norm_2(*z,*v2)
			 ,"sqrt(2.0*3.0*3.0*vecSpc.dim())",ST::magnitude(ST::squareroot(Scalar(2.0*3.0*3.0)*Scalar(vecSpc.dim())))
			 ,"maxRelErr",maxRelErr(),out
			 )
		) success=false;

	if(!VectorStdOpsTesterComparable<ST::isComparable,Scalar>::checkComparableStdOps(vecSpc,&*z,maxRelErr(),out,dumpAll)) success=false;

	// ToDo: Add tests for *all* standard operators!

  if(out) *out
		<< "\n*** Leaving VectorStdOpsTester<"<<ST::name()<<">::checkStdOps(...) ...\n";

	return success;

}

} // namespace TSFCore

#endif // TSFCORE_VECTOR_STD_OPS_TESTER_HPP
