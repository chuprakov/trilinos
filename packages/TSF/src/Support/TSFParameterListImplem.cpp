#include "TSFParameterList.h"
#include "TSFParameterListImplem.h"



#include "TSFError.h"

using namespace TSF;
using namespace std;

TSFParameterListImplem::TSFParameterListImplem(const string& label)
	: label_(label), children_(), parameters_()
{}

void TSFParameterListImplem::addChild(const TSFParameterList& list) 
{
	/* check that the key does not already exist in the list */
	if (!children_.containsKey(list.getLabel()))
		{
			children_.put(list.getLabel(), list);
		}
	else 
		{
			TSFError::raise("key=" + list.getLabel() + " previously defined");
		}
}

void TSFParameterListImplem::setChild(const TSFParameterList& list) 
{
  children_.put(list.getLabel(), list);
}

void TSFParameterListImplem::getChild(const string& name, TSFParameterList& list) const 
{
	if (!children_.containsKey(list.getLabel()))
		{
			TSFError::raise("ParameterList child label=" + name + " does not exist");
		}
	else 
		{
			list = children_.get(list.getLabel());
		}
}

TSFArray<TSFParameterList> TSFParameterListImplem::listChildren() const
{
	TSFArray<TSFParameterList> values;
	TSFArray<string> keys;

	children_.arrayify(keys, values);

	return values;
}

TSFArray<string> TSFParameterListImplem::listChildNames() const
{
	TSFArray<TSFParameterList> values;
	TSFArray<string> keys;

	children_.arrayify(keys, values);

	return keys;
}

TSFArray<TSFParameter> TSFParameterListImplem::listParameters() const
{
	TSFArray<TSFParameter> values;
	TSFArray<string> keys;

	parameters_.arrayify(keys, values);

	return values;
}

TSFArray<string> TSFParameterListImplem::listParameterNames() const
{
	TSFArray<TSFParameter> values;
	TSFArray<string> keys;

	parameters_.arrayify(keys, values);

	return keys;
}

void TSFParameterListImplem::addParameter(const TSFParameter& parameter)
{
		/* check that the key does not already exist in the list */
	if (!parameters_.containsKey(parameter.getLabel()))
		{
			parameters_.put(parameter.getLabel(), parameter);
		}
	else 
		{
			TSFError::raise("key=" + parameter.getLabel() + " previously defined in parameter list");
		}
}

void TSFParameterListImplem::setParameter(const TSFParameter& parameter)
{
  parameters_.put(parameter.getLabel(), parameter);
}


TSFParameter TSFParameterListImplem::getParameter(const string& name) const 
{
	if (!parameters_.containsKey(name))
		{
			TSFError::raise("ParameterList parameter label=" + name + " does not exist");
		}
	return parameters_.get(name);
}

void TSFParameterListImplem::setValue(const string& name, char value) 
{
	getParameter(name).set(value);
}

void TSFParameterListImplem::setValue(const string& name, const string& value) 
{
	getParameter(name).set(value);
}

void TSFParameterListImplem::setValue(const string& name, int value) 
{
	getParameter(name).set(value);
}

void TSFParameterListImplem::setValue(const string& name, const TSFArray<int>& value) 
{
	getParameter(name).set(value);
}

void TSFParameterListImplem::setValue(const string& name, double value) 
{
	getParameter(name).set(value);
}
void TSFParameterListImplem::setValue(const string& name, const TSFArray<double>& value) 
{
	getParameter(name).set(value);
}

void TSFParameterListImplem::print(ostream& os, int indentDepth) const 
{
  os << TSFObject::spaces(indentDepth) << label_ << endl;
  TSFArray<TSFParameter> p = listParameters();
  for (int i=0; i<p.length(); i++)
    {
      p[i].printIndented(os, indentDepth+1);
    }
  TSFArray<TSFParameterList> c = listChildren();
  
  for (int i=0; i<c.length(); i++)
    {
      c[i].printIndented(os, indentDepth+1);
    }
}
