#include "DenseSerialVector.h"
#include "TSFTimeMonitor.h"
#include "TSFUtils.h"
#include "TSFBlas.h"

namespace TSF
{

	int DenseSerialVector::one_ = 1;
	TSFReal DenseSerialVector::onePointZero_ = 1.0;
	TSFReal DenseSerialVector::negativeOnePointZero_ = -1.0;
}

using namespace TSF;

TSFTimer DenseSerialVector::blasTimer_("DSF BLAS");

DenseSerialVector& DenseSerialVector::operator=(const DenseSerialVector& other)
{
	if (this != &other) 
		{
			resize(other.n_);
			TSFBlas<TSFReal>::copy(&n_, other.x_, &one_, x_, &one_);
		}
	return *this;
}
	


void DenseSerialVector::setScalar(const TSFReal& a)
{
	TSFReal* yy = x_;
	for (int i=0; i<n_; i++, yy++)
		{
			*yy = a;
		}
}



void DenseSerialVector::negate()
{
	TSFReal* yy = x_;
	for (int i=0; i<n_; i++, yy++)
		{
			*yy = -(*yy);
		}
}

void DenseSerialVector::add(const DenseSerialVector& other)
{
	if (n_ != other.n_) 
		{
			TSFError::raise("DenseSerialVector::add length mismatch");
		}
	TSFTimeMonitor timer(blasTimer_);
	TSFBlas<TSFReal>::axpy(&n_, &onePointZero_, other.x_, &one_, x_, &one_);
}



void DenseSerialVector::subtract(const DenseSerialVector& other)
{
	if (n_ != other.n_) 
		{
			TSFError::raise("DenseSerialVector::subtract length mismatch");
		}
	TSFTimeMonitor timer(blasTimer_);
	TSFBlas<TSFReal>::axpy(&n_, &negativeOnePointZero_, 
												 other.x_, &one_, x_, &one_);
}


void DenseSerialVector::daxpy(const DenseSerialVector& other, const TSFReal& a)
{
	if (n_ != other.n_) 
		{
			TSFError::raise("DenseSerialVector::daxpy length mismatch");
		}
	TSFTimeMonitor timer(blasTimer_);
	TSFBlas<TSFReal>::axpy(&n_, &a, other.x_, &one_, x_, &one_);
}


void DenseSerialVector::eMult(const DenseSerialVector& other)
{
	if (n_ != other.n_) 
		{
			TSFError::raise("DenseSerialVector::eMult length mismatch");
		}
	TSFReal* yy = x_;
	TSFReal* xx = other.x_;	
	for (int i=0; i<n_; i++)
		{
			*yy++ *= *xx++;
		}
}

void DenseSerialVector::scalarMult(const TSFReal& a)
{
	TSFBlas<TSFReal>::scal(&n_, &a, x_, &one_);
}

void DenseSerialVector::scalarPow(const TSFReal& a)
{
	if (TSFUtils::chop(a-1.0)==0) return;

	TSFReal* yy = x_;

	if (TSFUtils::chop(a+1.0)==0)
		{ 
			for (int i=0; i<n_; i++)
				{
					*yy = 1.0/(*yy);
					yy++;
				}
		}
	else if (TSFUtils::chop(a-2.0)==0)
		{ 
			for (int i=0; i<n_; i++)
				{
					*yy = *yy*(*yy);
					yy++;
				}
		}
	else if (TSFUtils::chop(a-3.0)==0)
		{ 
			for (int i=0; i<n_; i++)
				{
					TSFReal tmp = *yy*(*yy);
					*yy++ *= tmp;
				}
		}
	else if (TSFUtils::chop(a-4.0)==0)
		{ 
			for (int i=0; i<n_; i++)
				{
					TSFReal tmp = *yy*(*yy);
					*yy++ = tmp*tmp;
				}
		}
	else
		{
			for (int i=0; i<n_; i++)
				{
					*yy++ = pow(*yy, a);
				}
		}
}


TSFReal DenseSerialVector::dot(const DenseSerialVector& other) const 
{
	if (n_ != other.n_) 
		{
			TSFError::raise("DenseSerialVector::dot length mismatch");
		}
	TSFReal rtn = 0.0;
	TSFReal* xx = x_;
	TSFReal* yy = other.x_;
	
	rtn = TSFBlas<TSFReal>::dot(&n_, xx, &one_, yy, &one_);
	return rtn;
}

TSFReal DenseSerialVector::norm2Squared() const
{
	TSFReal rtn = 0.0;

	rtn = TSFBlas<TSFReal>::nrm2(&n_, x_, &one_);

	return rtn*rtn;
}


TSFReal DenseSerialVector::sumElements() const
{
	TSFReal rtn = 0.0;
	TSFReal* xx = x_;


	for (int i=0; i<n_; i++)
		{
			rtn += *(xx++);
		}
	return rtn;
}

TSFReal DenseSerialVector::maxNorm() const
{
	TSFReal mx = -1.0e30;
	TSFReal* xx = x_;

	for (int i=0; i<n_; i++)
		{
			if (*xx > mx) mx = *xx;
			xx++;
		}
	return mx;
}

string DenseSerialVector::toString() const
{
	/* BVBW changed from strstream to std::ostringstream */
#if HAVE_STRSTREAM
	std::ostrstream str;
#else
	std::ostringstream str;
#endif
	str << *this << ends;
	return str.str();
}

string DenseSerialVector::summary() const
{
	return "DenseSerialVector[dim=" + TSFUtils::toString(n_) + "]";
}

namespace TSF
{

	ostream& operator<<(ostream& os, const DenseSerialVector& x)
	{
		os << "[";
		for (int i=0; i<x.length(); i++)
			{
				os << x[i];
				if (i<(x.length()-1)) os << ", ";
			}
		os << "]";
		return os;
	}


}

			









