#ifndef TSFPARAMETERLISTIMPLEM_H
#define TSFPARAMETERLISTIMPLEM_H

#include "TSFConfig.h"
#include "TSFParameterList.h"
#include <string>
#include "TSFArray.h"
#include "TSFHashtable.h"

namespace TSF
{

  using std::string;
  using std::ostream;

  class TSFParameterListImplem
    {
    public:
      /** construct with a descriptive label */
      TSFParameterListImplem(const string& label);

      /** return descriptive label */
      const string& getLabel() const {return label_;}

      //@{ \name Child lists
      /** Add a new child list. The child's name is expected to be unique. */
      void addChild(const TSFParameterList& list);

      /** Add a new child list. The child's name is expected to be unique. */
      void setChild(const TSFParameterList& list);

      /** Look up a child list by name */
      void getChild(const string& name, TSFParameterList& list) const ;

      /** get list of all children */
      TSFArray<TSFParameterList> listChildren() const ;
      //@}


      //@{ \name Creating parameters

      /** Add a new parameter to the list. The parameter's name is checked against
       * the list of existing names to ensure that each name appears only once.
       * If the parameter's name already appears, an error is thrown. To change the
       * value of an existing parameter, use the setValue() method instead. */
      void addParameter(const TSFParameter& parameter);

      /** Add a new parameter to the list. The parameter's name is checked against
       * the list of existing names to ensure that each name appears only once.
       * If the parameter's name already appears, an error is thrown. To change the
       * value of an existing parameter, use the setValue() method instead. */
      void setParameter(const TSFParameter& parameter);

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
      void print(ostream& os, int indentDepth) const ;

    private:
      /* the label */
      string label_;

      /* child lists of this list, keyed by label. */
      TSFHashtable<string, TSFParameterList> children_;

      /* parameters, keyed by label. */
      TSFHashtable<string, TSFParameter> parameters_;
    };
}

#endif
