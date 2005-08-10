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

#ifndef TSFLINEARCOMBINATION_HPP
#define TSFLINEARCOMBINATION_HPP

#include "TSFConfigDefs.hpp"
#include "TSFVector.hpp"
#include "TSFLinearOperatorDecl.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace TSFExtendedOps
{
  /** 
   * 
   */
  enum LCSign {LCAdd = 1, LCSubtract = -1};

  template <class Scalar>
  class ConvertibleToVector
  {
  public:
    virtual ~ConvertibleToVector(){;}
    /** */
    virtual TSFExtended::Vector<Scalar> eval() const = 0 ;

    /** */
    Scalar operator*(const ConvertibleToVector<Scalar>& other) const
    {return eval().dot(other.eval());}

    /** */
    Scalar norm1() const {return eval().norm1();}

    /** */
    Scalar norm2() const {return eval().norm2();}

    /** */
    Scalar normInf() const {return eval().normInf();}
  };

  /** 
   * Class LC1 holds the information required to perform
   * a 1-term linear combination, i.e., a scalar times a vector
   * or an operator times a vectir.
   */
  template <class Scalar>
  class LC1 : public ConvertibleToVector<Scalar>
  {
    public:
    /** */
    LC1(const TSFExtended::Vector<Scalar>& x) 
      : alpha_(1.0), op_(), x_(x) {;}

    /** */
    LC1(const Scalar& alpha, const TSFExtended::Vector<Scalar>& x)
      : alpha_(alpha), op_(), x_(x) {;}

    /** */
    LC1(const Scalar& alpha,
        const TSFExtended::LinearOperator<Scalar>& op, 
        const TSFExtended::Vector<Scalar>& x)
      : alpha_(alpha), op_(op), x_(x) {;}

    /** 
     * Evaluate the term into the argument vector, overwriting 
     * the previous value of the argument. */
    void evalInto(TSFExtended::Vector<Scalar>& result) const ;

    /** Add the term into the argument vector */
    void addInto(TSFExtended::Vector<Scalar>& result, 
                 LCSign sign = LCAdd) const ;

    /** Evaluate the term and return its value */
    virtual TSFExtended::Vector<Scalar> eval() const ;

    /** Determine whether this term contains the given vector */
    bool containsVector(const Thyra::VectorBase<Scalar>* vec) const 
    {return vec == x_.ptr().get();}

    /** */
    LC1<Scalar> operator*(const Scalar& beta) const ;
    
    private:
    Scalar alpha_;

    TSFExtended::LinearOperator<Scalar> op_;

    TSFExtended::Vector<Scalar> x_;
  };

  /**
   * Class LCN contains an n-term linear combination. 
   */
  template <class Scalar, class Node1, class Node2>
  class LCN : public ConvertibleToVector<Scalar>
  {
    public:
    /** */
    LCN(const Node1& x1, const Node2& x2, LCSign sign = LCAdd)
      : x1_(x1), x2_(x2), sign_(sign) {;}

    /** */
    operator TSFExtended::Vector<Scalar>() const {return eval();}
      
    /** */
    void evalInto(TSFExtended::Vector<Scalar>& result) const ;

    /** */
    void addInto(TSFExtended::Vector<Scalar>& result, 
                 LCSign sign = LCAdd) const ;

    /** */
    virtual TSFExtended::Vector<Scalar> eval() const ;

    /** */
    bool containsVector(const Thyra::VectorBase<Scalar>* vec) const
    {return x1_.containsVector(vec) || x2_.containsVector(vec);}
    
    private:
    Node1 x1_;

    Node2 x2_;

    LCSign sign_;
  };


  template <class Scalar> inline
  LC1<Scalar> operator*(const Scalar& alpha, 
                        const TSFExtended::Vector<Scalar>& x)
  {
    return LC1<Scalar>(alpha, x);
  } 

  template <class Scalar> inline
  LC1<Scalar> operator*(const TSFExtended::Vector<Scalar>& x, 
                        const Scalar& alpha)
  {
    return LC1<Scalar>(alpha, x);
  }

  template <class Scalar> inline
  LC1<Scalar> operator*(const TSFExtended::LinearOperator<Scalar>& op, 
                        const TSFExtended::Vector<Scalar>& x)
  {
    return LC1<Scalar>(1.0, op, x);
  }

  template <class Scalar> inline
  LC1<Scalar> LC1<Scalar>::operator*(const Scalar& beta) const
  {
    return LC1<Scalar>(alpha_*beta, op_, x_);
  }



  template <class Scalar> inline
  LCN<Scalar, LC1<Scalar>, LC1<Scalar> >
  operator+(const TSFExtended::Vector<Scalar>& x1, 
            const TSFExtended::Vector<Scalar>& x2)
  {
    return LCN<Scalar, LC1<Scalar>, LC1<Scalar> >(x1, x2);
  }
  
  template <class Scalar> inline
  LCN<Scalar, LC1<Scalar>, LC1<Scalar> >
  operator-(const TSFExtended::Vector<Scalar>& x1, 
            const TSFExtended::Vector<Scalar>& x2)
  {
    return LCN<Scalar, LC1<Scalar>, LC1<Scalar> >(x1, x2, LCSubtract);
  }





  template <class Scalar> inline
  LCN<Scalar, LC1<Scalar>, LC1<Scalar> >
  operator+(const LC1<Scalar>& x1, 
            const TSFExtended::Vector<Scalar>& x2)
  {
    return LCN<Scalar, LC1<Scalar>, LC1<Scalar> >(x1, x2);
  }

  template <class Scalar> inline
  LCN<Scalar, LC1<Scalar>, LC1<Scalar> >
  operator-(const LC1<Scalar>& x1, 
            const TSFExtended::Vector<Scalar>& x2)
  {
    return LCN<Scalar, LC1<Scalar>, LC1<Scalar> >(x1, x2, LCSubtract);
  }



  template <class Scalar> inline
  LCN<Scalar, LC1<Scalar>, LC1<Scalar> >
  operator+(const TSFExtended::Vector<Scalar>& x1, const LC1<Scalar>& x2)
  {
    return LCN<Scalar, LC1<Scalar>, LC1<Scalar> >(x1, x2);
  }

  template <class Scalar> inline
  LCN<Scalar, LC1<Scalar>, LC1<Scalar> >
  operator-(const TSFExtended::Vector<Scalar>& x1, const LC1<Scalar>& x2)
  {
    return LCN<Scalar, LC1<Scalar>, LC1<Scalar> >(x1, x2, LCSubtract);
  }

  
  
  template <class Scalar, class Node1, class Node2> inline
  LCN<Scalar, LC1<Scalar>, LCN<Scalar, Node1, Node2> >
  operator+(const TSFExtended::Vector<Scalar>& x1, 
            const LCN<Scalar, Node1, Node2>& x2)
  {
    return LCN<Scalar, LC1<Scalar>, LCN<Scalar, Node1, Node2> >(x1, x2);
  }

  template <class Scalar, class Node1, class Node2> inline
  LCN<Scalar, LC1<Scalar>, LCN<Scalar, Node1, Node2> >
  operator-(const TSFExtended::Vector<Scalar>& x1, 
            const LCN<Scalar, Node1, Node2>& x2)
  {
    return LCN<Scalar, LC1<Scalar>, LCN<Scalar, Node1, Node2> >(x1, x2, 
                                                                LCSubtract);
  }



  template <class Scalar, class Node1, class Node2> inline
  LCN<Scalar, LCN<Scalar, Node1, Node2>, LC1<Scalar> >
  operator+(const LCN<Scalar, Node1, Node2>& x1, 
            const TSFExtended::Vector<Scalar>& x2)
  {
    return LCN<Scalar, LCN<Scalar, Node1, Node2>, LC1<Scalar> >(x1, x2);
  }

  template <class Scalar, class Node1, class Node2> inline
  LCN<Scalar, LCN<Scalar, Node1, Node2>, LC1<Scalar> >
  operator-(const LCN<Scalar, Node1, Node2>& x1, 
            const TSFExtended::Vector<Scalar>& x2)
  {
    return LCN<Scalar, LCN<Scalar, Node1, Node2>, LC1<Scalar> >(x1, x2, 
                                                                LCSubtract);
  }




  template <class Scalar> inline
  LCN<Scalar, LC1<Scalar>, LC1<Scalar> > 
  operator+(const LC1<Scalar>& x1, const LC1<Scalar>& x2)
  {
    return LCN<Scalar, LC1<Scalar>, LC1<Scalar> >(x1, x2);
  }

  template <class Scalar> inline
  LCN<Scalar, LC1<Scalar>, LC1<Scalar> > 
  operator-(const LC1<Scalar>& x1, const LC1<Scalar>& x2)
  {
    return LCN<Scalar, LC1<Scalar>, LC1<Scalar> >(x1, x2, LCSubtract);
  }



  template <class Scalar, class Node1, class Node2> inline
  LCN<Scalar, LC1<Scalar>, LCN<Scalar, Node1, Node2> > 
  operator+(const LC1<Scalar>& x1, const LCN<Scalar, Node1, Node2>& x2)
  {
    return LCN<Scalar, LC1<Scalar>, LCN<Scalar, Node1, Node2> >(x1, x2);
  }

  template <class Scalar, class Node1, class Node2> inline
  LCN<Scalar, LC1<Scalar>, LCN<Scalar, Node1, Node2> > 
  operator-(const LC1<Scalar>& x1, const LCN<Scalar, Node1, Node2>& x2)
  {
    return LCN<Scalar, LC1<Scalar>, LCN<Scalar, Node1, Node2> >(x1, x2, 
                                                                LCSubtract);
  }
  


  template <class Scalar, class Node1, class Node2> inline
  LCN<Scalar, LCN<Scalar, Node1, Node2>, LC1<Scalar> > 
  operator+(const LCN<Scalar, Node1, Node2>& x1, const LC1<Scalar>& x2)
  {
    return LCN<Scalar, LCN<Scalar, Node1, Node2>, LC1<Scalar> >(x1, x2);
  }

  template <class Scalar, class Node1, class Node2> inline
  LCN<Scalar, LCN<Scalar, Node1, Node2>, LC1<Scalar> > 
  operator-(const LCN<Scalar, Node1, Node2>& x1, const LC1<Scalar>& x2)
  {
    return LCN<Scalar, LCN<Scalar, Node1, Node2>, LC1<Scalar> >(x1, x2, 
                                                                LCSubtract);
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

  template <class Scalar, class Node1, class Node2, 
            class Node3, class Node4> inline
  LCN<Scalar, LCN<Scalar, Node1, Node2>, LCN<Scalar, Node3, Node4> >
  operator-(const LCN<Scalar, Node1, Node2>& x1, 
            const LCN<Scalar, Node3, Node4>& x2)
  {
    return LCN<Scalar, LCN<Scalar, Node1, Node2>, 
      LCN<Scalar, Node3, Node4> >(x1, x2, LCSubtract);
  }




  template <class Scalar, class Node1, class Node2> inline
  void LCN<Scalar, Node1, Node2>::evalInto(TSFExtended::Vector<Scalar>& result) const
  {
    x1_.evalInto(result);
    x2_.addInto(result, sign_);
  } 

  template <class Scalar, class Node1, class Node2> inline
  void LCN<Scalar, Node1, Node2>::addInto(TSFExtended::Vector<Scalar>& result,
                                          LCSign sign) const
  {
    x1_.addInto(result, sign);
    x2_.addInto(result, sign_ * sign);
  }

  template <class Scalar, class Node1, class Node2> inline
  TSFExtended::Vector<Scalar> LCN<Scalar, Node1, Node2>::eval() const
  {
    TSFExtended::Vector<Scalar> result = x1_.eval();
    x2_.addInto(result, sign_);
    return result;
  }
  
  
  template <class Scalar> inline
  void LC1<Scalar>::evalInto(TSFExtended::Vector<Scalar>& result) const
  {
    if (op_.ptr().get() != 0)
      {
        op_.apply(x_, result, alpha_);
      }
    else
      {
        result.acceptCopyOf(x_);
        result.scale(alpha_);
      }
  }

  template <class Scalar> inline
  void LC1<Scalar>::addInto(TSFExtended::Vector<Scalar>& result,
                            LCSign sign) const
  {
    if (op_.ptr().get() != 0)
      {
        op_.apply(x_, result, alpha_, 1.0);
      }
    else
      {
        result.update(sign*alpha_, x_);
      }
  } 

  template <class Scalar> inline
  TSFExtended::Vector<Scalar> LC1<Scalar>::eval() const 
  {
    TSFExtended::Vector<Scalar> result;

    if (op_.ptr().get() != 0)
      {
        result = op_.range().createMember();
        op_.apply(x_, result, alpha_);
      }
    else
      {
        result = x_.copy();
        result.scale(alpha_);    
      }

    return result;
  }
}

namespace TSFExtended
{
  /* definition of assignment from 1-term linear combination to a vector */
  template <class Scalar> inline
  Vector<Scalar>& Vector<Scalar>::operator=(const TSFExtendedOps::LC1<Scalar>& x)
  {
    if (this->ptr().get()==0)
      {
        *this = x.eval();
      }
    else if (x.containsVector(this->ptr().get()))
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
  template <class Node1, class Node2> inline
  Vector<Scalar>& Vector<Scalar>::operator=(const TSFExtendedOps::LCN<Scalar, Node1, Node2>& x)
  {
    if (this->ptr().get()==0)
      {
        *this = x.eval();
      }
    else if (x.containsVector(this->ptr().get()))
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

  template <class Scalar>
  template <class Node1, class Node2> inline
  Vector<Scalar>::Vector(const TSFExtendedOps::LCN<Scalar, Node1, Node2>& x)
    : Handle<Thyra::VectorBase<Scalar> >(x.eval().ptr())
  {;}

  template <class Scalar> inline
  Vector<Scalar>::Vector(const TSFExtendedOps::LC1<Scalar>& x)
    : Handle<Thyra::VectorBase<Scalar> >(x.eval().ptr())
  {;}

}

#endif  /* DOXYGEN_DEVELOPER_ONLY */


#endif
