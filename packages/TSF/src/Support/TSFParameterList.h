#ifndef TSFPARAMETERLIST_H
#define TSFPARAMETERLIST_H

#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include "TSFParameter.h"
#include <string>
#include "TSFArray.h"

namespace TSF
{

  using namespace std;

  /** forward declare class TSFParameterListImplem */
  class TSFParameterListImplem;

  /** \ingroup Utilities
   *
   */
  class TSFParameterList : public TSFObject
    {
    public:
      /** empty ctor, for use in containers */
      TSFParameterList();

      /** Virtual dtor */
      virtual ~TSFParameterList();

      /** construct with a descriptive label */
      TSFParameterList(const string& label);

      /** return descriptive label */
      const string& getLabel() const ;

      //@{ \name Child lists
      /** Add a new child list. The child's name is expected to be unique. */
      void addChild(const TSFParameterList& list);

      /** */
      void setChild(const TSFParameterList& list);

      /** Look up a child list by name */
      TSFParameterList getChild(const string& name) const ;

      /** get list of all children */
      TSFArray<TSFParameterList> listChildren() const ;
      //@}


      //@{ \name Creating parameters

      /** Add a new parameter to the list. The parameter's name is checked against
       * the list of existing names to ensure that each name appears only once.
       * If the parameter's name already appears, an error is thrown. To change the
       * value of an existing parameter, use the setValue() method instead. */
      void addParameter(const TSFParameter& parameter);

      //@}


      //@{ \name Setting parameter values

      /** set the value of an existing character parameter.  The
       * parameter's name is expected to be in the list; if not, an
       * error is thrown. To add a new parameter to the list, use the
       * add() method instead. */
      void setValue(const string& name, char value);

      /** set the value of an existing string parameter.  The
       * parameter's name is expected to be in the list; if not, an
       * error is thrown. To add a new parameter to the list, use the
       * add() method instead. */
      void setValue(const string& name, const string& value);

      /** set the value of an existing integer parameter.  The
       * parameter's name is expected to be in the list; if not, an
       * error is thrown. To add a new parameter to the list, use the
       * add() method instead. */
      void setValue(const string& name, int value);


      /** set the value of an existing double parameter.  The
       * parameter's name is expected to be in the list; if not, an
       * error is thrown. To add a new parameter to the list, use the
       * add() method instead. */
      void setValue(const string& name, double value);

      /** set the value of an existing int array parameter.  The
       * parameter's name is expected to be in the list; if not, an
       * error is thrown. To add a new parameter to the list, use the
       * add() method instead. */
      void setValue(const string& name, const TSFArray<int>& value);

      /** set the value of an existing double parameter.  The
       * parameter's name is expected to be in the list; if not, an
       * error is thrown. To add a new parameter to the list, use the
       * add() method instead. */
      void setValue(const string& name, const TSFArray<double>& value);

      /** */
      void setParameter(const TSFParameter& param);

      //@}

      //@{ \name Getting parameters
      /** Look up a parameter by name */
      TSFParameter getParameter(const string& name) const ;

      /** Return a list of all parameters */
      TSFArray<TSFParameter> listParameters() const ;
      //@}

      //@{ \name Getting information about labels
      /** Return a list of all parameter names */
      TSFArray<string> listParameterNames() const ;

      /** Return a list of all child names */
      TSFArray<string> listChildNames() const ;
      //@}

      /** */
      virtual TSFParameterList overrideWith(const TSFParameterList& newParams) const ;

      /** */
      virtual void print(ostream& os) const ;

    private:
      TSFUnsharedPtr<TSFParameterListImplem> ptr_;
    };

}

#endif
