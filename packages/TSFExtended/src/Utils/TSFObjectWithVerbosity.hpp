/* @HEADER@ */
/* @HEADER@ */

#ifndef TSF_OBJECTWITHVERBOSITY_H
#define TSF_OBJECTWITHVERBOSITY_H

#include "TSFConfigDefs.hpp"


namespace TSFExtended
{
  /** \enum Verbosity settings */
  enum VerbositySetting {VerbSilent=0, VerbLow=1, VerbMedium=2, 
                         VerbHigh=3, VerbExtreme=4};
  
#ifndef DOXYGEN_DEVELOPER_ONLY
  /**
   * ObjectWithVerbosity and the related verbosity() method
   * provide an interface for getting/setting
   * verbosity flags for classes or instances. 
   *
   * All objects start out with a verbosity setting of VerbSilent.
   * 
   * You can set verbosity for a single instance of a class, or for
   * the whole class. To set for an instance, use the verbosity()
   * member function, for example,
   * \code
   * Mesh mesh1 = reader1.getMesh();
   * Mesh mesh2 = reader2.getMesh();
   * Mesh mesh3 = reader3.getMesh();
   * mesh1.verbosity() = VerbHigh;
   * \endcode
   * which sets the verbosity of <tt>mesh1</tt> to VerbHigh and leaves
   * those of <tt>mesh2</tt> and <tt>mesh3</tt> unchanged.
   *
   * Alternatively, you can set a default verbosity for an entire
   * class, for example,
   * \code
   * Mesh mesh1 = reader1.getMesh();
   * Mesh mesh2 = reader2.getMesh();
   * Mesh mesh3 = reader3.getMesh();
   * mesh1.verbosity() = VerbHigh;
   * verbosity<Mesh>() = VerbMedium;
   * \endcode
   * which sets the default verbosity to VerbMedium. Since <tt>mesh1</tt>
   * has its own verbosity setting of VerbHigh, 
   * it will use it rather than the
   * default, but <tt>mesh2</tt> and <tt>mesh3</tt> will use VerbMedium.
   * 
   */
  template <class X>
  class ObjectWithVerbosity
  {
  public:
    /** Construct, starting silent */
    ObjectWithVerbosity() : verbosity_(VerbSilent), setLocally_(false) {;}

    /** Read-only access to the verbosity */
    VerbositySetting verbosity() const 
    {
      if (setLocally_) return verbosity_;
      return classVerbosity();
    }

    /** Writeable access to the verbosity setting */
    VerbositySetting& verbosity() 
    {
      setLocally_ = true; 
      return verbosity_;
    }

    /** Writeable access to the default verbosity for the class */
    static VerbositySetting& classVerbosity() 
    {
      static VerbositySetting rtn = VerbSilent;
      return rtn;
    }
  private:
    /** */
    VerbositySetting verbosity_;

    /** */
    bool setLocally_;
  };

  /** 
   * \relates ObjectWithVerbosity
   * Global method for setting verbosity of a class
   */
  template <class X> VerbositySetting& verbosity() 
  {
    return X::classVerbosity();
  }
#endif  /* DOXYGEN_DEVELOPER_ONLY */   
}





#endif
