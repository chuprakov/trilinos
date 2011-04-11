/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_VECTOROPSIMPL_HPP
#define PLAYA_VECTOROPSIMPL_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorImpl.hpp"
#include "PlayaVectorOpsDecl.hpp"
#include "Teuchos_MPIComm.hpp"
#include "Teuchos_RCP.hpp"

namespace PlayaFunctors
{
template <class Scalar> class BoundedMinLocFunctor;

template <class Scalar> class BoundedMaxLocFunctor;

template <class Scalar> class Norm2Dist;

template <class Scalar> class Norm1Dist;

template <class Scalar> class NormInfDist;
}

namespace Playa
{

using namespace PlayaFunctors;

/* */
template <class Scalar> inline
Scalar minloc(const Vector<Scalar>& x, int& gni)
{
  return minlocWithBound(-HUGE_VAL, x, gni);
}

/* */
template <class Scalar> inline
Scalar maxloc(const Vector<Scalar>& x, int& gni)
{
  return maxlocWithBound(-HUGE_VAL, x, gni);
}

/* */
template <class Scalar> inline
Scalar minlocWithBound(const Scalar& lowerBound, 
  const Vector<Scalar>& x, int& gni)
{
  IndexedValue<Scalar> y = 
    x.applyUnaryReductionFunctor(BoundedMinLocFunctor<Scalar>(x.comm(), lowerBound, x.space().baseGlobalNaturalIndex()));
  gni = y.where;
  return y.what;
}

/* */
template <class Scalar> inline
Scalar maxlocWithBound(const Scalar& upperBound, 
  const Vector<Scalar>& x, int& gni)
{
  IndexedValue<Scalar> y = 
    x.applyUnaryReductionFunctor(BoundedMaxLocFunctor<Scalar>(x.comm(), upperBound, x.space().baseGlobalNaturalIndex()));
  gni = y.where;
  return y.what;
}


/* */
template <class Scalar>
Scalar norm2Dist(const Vector<Scalar>& x, const Vector<Scalar>& y)
{
  return x.applyBinaryReductionFunctor(
    PlayaFunctors::Norm2Dist<Scalar>(x.comm()), y);
}

/* */
template <class Scalar>
Scalar norm1Dist(const Vector<Scalar>& x, const Vector<Scalar>& y)
{
  return x.applyBinaryReductionFunctor(
    PlayaFunctors::Norm1Dist<Scalar>(x.comm()), y);
}

/* */
template <class Scalar>
Scalar normInfDist(const Vector<Scalar>& x, const Vector<Scalar>& y)
{
  return x.applyBinaryReductionFunctor(
    PlayaFunctors::NormInfDist<Scalar>(x.comm()), y);
}

/** \relates Vector \brief Compute the Euclidean norm of a vector */
template <class Scalar>
Scalar norm2(const Vector<Scalar>& x) {return x.norm2();}

/** \relates Vector \brief Compute the one-norm of a vector */
template <class Scalar>
Scalar norm1(const Vector<Scalar>& x) {return x.norm1();}

/** \relates Vector \brief Compute the infinity norm of a vector */
template <class Scalar>
Scalar normInf(const Vector<Scalar>& x) {return x.normInf();}

} // end namespace Playa

namespace PlayaFunctors
{


using namespace Playa;

/** \brief Find minimum element above a lower bound, returning value and
 * location */
template <class Scalar>
class BoundedMinLocFunctor : public ReductionFunctorBase<Scalar>
{
public:
  /** */
  BoundedMinLocFunctor(const MPIComm& comm, const Scalar& bound,
    int baseGNI)
    : ReductionFunctorBase<Scalar>(comm), min_(), 
      bound_(bound), baseGNI_(baseGNI)
    {
      min_.what = HUGE_VAL;
      min_.where = -1;
    }

  /** */
  void step(int i, const Scalar& x) const 
    {
      if (x < min_.what && x > bound_) 
      {
        min_.what = x;
        min_.where = i;
      }
    }

  /** */
  void postProc() const 
    {
      min_.where += baseGNI_;

      IndexedValue<Scalar> out = min_;

      this->comm().allReduce(&min_, &out, 1, MPIComm::DOUBLE_INT, 
        MPIComm::MINLOC);
      min_ = out;
    }

  /** */
  IndexedValue<Scalar> result() const 
    {
      return min_;
    }

  /** */
  std::string description() const {return "MinLoc()";}

private:
  MPIComm comm_;
  mutable IndexedValue<Scalar> min_;
  Scalar bound_;
  int baseGNI_;
};


/** \brief Find maximum element below an upper bound, returning value and
 * location */
template <class Scalar>
class BoundedMaxLocFunctor : public ReductionFunctorBase<Scalar>
{
public:
  /** */
  BoundedMaxLocFunctor(const MPIComm& comm, const Scalar& bound,
    int baseGNI)
    : ReductionFunctorBase<Scalar>(comm), max_(), 
      bound_(bound), baseGNI_(baseGNI)
    {
      max_.what = -HUGE_VAL;
      max_.where = -1;
    }

  /** */
  void step(int i, const Scalar& x) const 
    {
      if (x > max_.what && x < bound_) 
      {
        max_.what = x;
        max_.where = i;
      }
    }

  /** */
  void postProc() const 
    {
      max_.where += baseGNI_;

      IndexedValue<Scalar> out = max_;

      this->comm().allReduce(&max_, &out, 1, MPIComm::DOUBLE_INT, 
        MPIComm::MAXLOC);
      max_ = out;
    }

  /** */
  IndexedValue<Scalar> result() const 
    {
      return max_;
    }

  /** */
  std::string description() const {return "MinLoc()";}

private:
  MPIComm comm_;
  mutable IndexedValue<Scalar> max_;
  Scalar bound_;
  int baseGNI_;
};



/** \brief Euclidean distance between two vectors */
template <class Scalar>
class Norm2Dist : public ReductionFunctorBase<Scalar>
{
public:
  Norm2Dist(const MPIComm& comm)
    : ReductionFunctorBase<Scalar>(comm), val_(0.0) {}

  void step(int i, const Scalar& x, const Scalar& y) const 
    {
      Scalar d = x-y;
      val_ += d*d;
    }

  void postProc() const 
    {
      Scalar final = val_;
      this->comm().allReduce(&val_, &final, 1, MPIComm::DOUBLE, MPIComm::SUM);
      val_ = final;
    }

  Scalar result() const 
    {
      return ::sqrt(val_);
    }

  /** */
  std::string description() const {return "Norm2Dist()";}

private:
  mutable Scalar val_;
};



/** \brief One-norm distance between two vectors */
template <class Scalar>
class Norm1Dist : public ReductionFunctorBase<Scalar>
{
public:
  Norm1Dist(const MPIComm& comm)
    : ReductionFunctorBase<Scalar>(comm), val_(0.0) {}

  void step(int i, const Scalar& x, const Scalar& y) const 
    {
      Scalar d = x-y;
      val_ += ::fabs(d);
    }

  void postProc() const 
    {
      Scalar final = val_;
      this->comm().allReduce(&val_, &final, 1, MPIComm::DOUBLE, MPIComm::SUM);
      val_ = final;
    }

  Scalar result() const 
    {
      return val_;
    }

  /** */
  std::string description() const {return "Norm1Dist()";}

private:
  mutable Scalar val_;
};

/** \brief Infinity-norm distance between two vectors */
template <class Scalar>
class NormInfDist : public ReductionFunctorBase<Scalar>
{
public:
  NormInfDist(const MPIComm& comm)
    : ReductionFunctorBase<Scalar>(comm), val_(0.0) {}

  void step(int i, const Scalar& x, const Scalar& y) const 
    {
      Scalar d = ::fabs(x-y);
      if (d > val_) val_ = d;
    }

  void postProc() const 
    {
      Scalar final = val_;
      this->comm().allReduce(&val_, &final, 1, MPIComm::DOUBLE, MPIComm::MAX);
      val_ = final;
    }

  Scalar result() const 
    {
      return val_;
    }

  /** */
  std::string description() const {return "NormInfDist()";}

private:
  mutable Scalar val_;
};


/** \brief Specify return type of BoundedMinLocFunctor */
template <class Scalar>
class VectorFunctorTraits<Scalar, BoundedMinLocFunctor<Scalar> >
{
public:
  typedef IndexedValue<Scalar> ReturnType ;
};


/** \brief Specify return type of BoundedMaxLocFunctor */
template <class Scalar>
class VectorFunctorTraits<Scalar, BoundedMaxLocFunctor<Scalar> >
{
public:
  typedef IndexedValue<Scalar> ReturnType ;
};

  
}

#endif
