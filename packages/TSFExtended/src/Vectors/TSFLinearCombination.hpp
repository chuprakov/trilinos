/* @HEADER@ */
/* @HEADER@ */

#ifndef TSFLINEARCOMBINATION_HPP
#define TSFLINEARCOMBINATION_HPP

#include "TSFConfigDefs.hpp"
#include "TSFVector.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace TSFExtendedOps
{
  /** 
   * Class LC1 holds the information required to perform
   * a 1-term linear combination, i.e., a scalar times a vector.
   */
  template <class Scalar>
  class LC1
  {
    public:
    /** */
    LC1(const TSFExtended::Vector<Scalar>& x) : alpha_(1.0), x_(x) {;}

    /** */
    LC1(const Scalar& alpha, const TSFExtended::Vector<Scalar>& x)
      : alpha_(alpha), x_(x) {;}

    /** 
     * Evaluate the term into the argument vector, overwriting 
     * the previous value of the argument. */
    void evalInto(TSFExtended::Vector<Scalar>& result) const ;

    /** Add the term into the argument vector */
    void addInto(TSFExtended::Vector<Scalar>& result) const ;

    /** */
    TSFExtended::Vector<Scalar> eval() const ;

    /** */
    bool containsVector(const TSFCore::Vector<Scalar>* vec) const 
    {return vec == x_.ptr().get();}

    private:
    Scalar alpha_;

    TSFExtended::Vector<Scalar> x_;
  };

  /**
   *
   */
  template <class Scalar, class Node1, class Node2>
  class LCN
  {
    public:
    /** */
    LCN(const Node1& x1, const Node2& x2)
      : x1_(x1), x2_(x2) {;}
      
    /** */
    void evalInto(TSFExtended::Vector<Scalar>& result) const ;

    /** */
    void addInto(TSFExtended::Vector<Scalar>& result) const ;

    /** */
    TSFExtended::Vector<Scalar> eval() const ;

    /** */
    bool containsVector(const TSFCore::Vector<Scalar>* vec) const
    {return x1_.containsVector(vec) || x2_.containsVector(vec);}
    
    private:
    Node1 x1_;

    Node2 x2_;
  };


  template <class Scalar> inline
  LC1<Scalar> operator*(const Scalar& alpha, const TSFExtended::Vector<Scalar>& x)
  {
    return LC1<Scalar>(alpha, x);
  }

  template <class Scalar> inline
  LCN<Scalar, LC1<Scalar>, LC1<Scalar> >
  operator+(const TSFExtended::Vector<Scalar>& x1, const TSFExtended::Vector<Scalar>& x2)
  {
    return LCN<Scalar, LC1<Scalar>, LC1<Scalar> >(x1, x2);
  }


  template <class Scalar> inline
  LCN<Scalar, LC1<Scalar>, LC1<Scalar> >
  operator+(const LC1<Scalar>& x1, const TSFExtended::Vector<Scalar>& x2)
  {
    return LCN<Scalar, LC1<Scalar>, LC1<Scalar> >(x1, x2);
  }

  template <class Scalar> inline
  LCN<Scalar, LC1<Scalar>, LC1<Scalar> >
  operator+(const TSFExtended::Vector<Scalar>& x1, const LC1<Scalar>& x2)
  {
    return LCN<Scalar, LC1<Scalar>, LC1<Scalar> >(x1, x2);
  }

  
  
  template <class Scalar, class Node1, class Node2> inline
  LCN<Scalar, LC1<Scalar>, LCN<Scalar, Node1, Node2> >
  operator+(const TSFExtended::Vector<Scalar>& x1, const LCN<Scalar, Node1, Node2>& x2)
  {
    return LCN<Scalar, LC1<Scalar>, LCN<Scalar, Node1, Node2> >(x1, x2);
  }

  template <class Scalar, class Node1, class Node2> inline
  LCN<Scalar, LCN<Scalar, Node1, Node2>, LC1<Scalar> >
  operator+(const LCN<Scalar, Node1, Node2>& x1, const TSFExtended::Vector<Scalar>& x2)
  {
    return LCN<Scalar, LCN<Scalar, Node1, Node2>, LC1<Scalar> >(x1, x2);
  }

  template <class Scalar> inline
  LCN<Scalar, LC1<Scalar>, LC1<Scalar> > 
  operator+(const LC1<Scalar>& x1, const LC1<Scalar>& x2)
  {
    return LCN<Scalar, LC1<Scalar>, LC1<Scalar> >(x1, x2);
  }

  template <class Scalar, class Node1, class Node2> inline
  LCN<Scalar, LC1<Scalar>, LCN<Scalar, Node1, Node2> > 
  operator+(const LC1<Scalar>& x1, const LCN<Scalar, Node1, Node2>& x2)
  {
    return LCN<Scalar, LC1<Scalar>, LCN<Scalar, Node1, Node2> >(x1, x2);
  }
  
  template <class Scalar, class Node1, class Node2> inline
  LCN<Scalar, LCN<Scalar, Node1, Node2>, LC1<Scalar> > 
  operator+(const LCN<Scalar, Node1, Node2>& x1, const LC1<Scalar>& x2)
  {
    return LCN<Scalar, LCN<Scalar, Node1, Node2>, LC1<Scalar> >(x1, x2);
  }

  
  template <class Scalar, class Node1, class Node2, 
            class Node3, class Node4> inline
  LCN<Scalar, LCN<Scalar, Node1, Node2>, LCN<Scalar, Node3, Node4> >
  operator+(const LCN<Scalar, Node1, Node2>& x1, 
            const LCN<Scalar, Node3, Node4>& x2)
  {
    return LCN<Scalar, LCN<Scalar, Node1, Node2>, 
      LCN<Scalar, Node3, Node4> >(x1, x2);
  }

  template <class Scalar, class Node1, class Node2> inline
  void LCN<Scalar, Node1, Node2>::evalInto(TSFExtended::Vector<Scalar>& result) const
  {
    x1_.evalInto(result);
    x2_.addInto(result);
  } 

  template <class Scalar, class Node1, class Node2> inline
  void LCN<Scalar, Node1, Node2>::addInto(TSFExtended::Vector<Scalar>& result) const
  {
    x1_.addInto(result);
    x2_.addInto(result);
  }

  template <class Scalar, class Node1, class Node2> inline
  TSFExtended::Vector<Scalar> LCN<Scalar, Node1, Node2>::eval() const
  {
    TSFExtended::Vector<Scalar> result = x1_.eval();
    x2_.addInto(result);
    return result;
  }
  
  
  template <class Scalar> inline
  void LC1<Scalar>::evalInto(TSFExtended::Vector<Scalar>& result) const
  {
    result.acceptCopyOf(x_).scale(alpha_);
  }

  template <class Scalar> inline
  void LC1<Scalar>::addInto(TSFExtended::Vector<Scalar>& result) const
  {
    result.update(alpha_, x_);
  } 

  template <class Scalar> inline
  TSFExtended::Vector<Scalar> LC1<Scalar>::eval() const 
  {
    TSFExtended::Vector<Scalar> result = x_.copy();
    result.scale(alpha_);
    
    return result;
  }
}

namespace TSFExtended
{
  /* definition of assignment from 1-term linear combination to a vector */
  template <class Scalar>
  Vector<Scalar>& Vector<Scalar>::operator=(const TSFExtendedOps::LC1<Scalar>& x)
  {
    if (ptr().get()==0)
      {
        *this = x.eval();
      }
    else if (x.containsVector(ptr().get()))
      {
        Vector<Scalar> rtn = x.eval();
        acceptCopyOf(rtn);
      }
    else
      {
        x.evalInto(*this);
      }
    return *this;
  }

 
  /* definition of assignment from N-term linear combination to a vector */
  template <class Scalar>
  template <class Node1, class Node2>
  Vector<Scalar>& Vector<Scalar>::operator=(const TSFExtendedOps::LCN<Scalar, Node1, Node2>& x)
  {
    if (ptr().get()==0)
      {
        *this = x.eval();
      }
    else if (x.containsVector(ptr().get()))
      {
        Vector<Scalar> rtn = x.eval();
        acceptCopyOf(rtn);
      }
    else
      {
        x.evalInto(*this);
      }
    return *this;
  }

}

#endif  /* DOXYGEN_DEVELOPER_ONLY */


#endif
