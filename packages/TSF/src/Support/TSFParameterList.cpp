#include "TSFParameterList.h"
#include "TSFParameterListImplem.h"

#include "TSFError.h"

using namespace TSF;
using namespace std;

TSFParameterList::TSFParameterList()
	: ptr_(0)
{}

TSFParameterList::TSFParameterList(const string& label)
	: ptr_(new TSFParameterListImplem(label))
{}

TSFParameterList::~TSFParameterList() {;}

const string& TSFParameterList::getLabel() const
{
	return ptr_->getLabel();
}


/* ------------- setting and getting children  ------------------------------*/

void TSFParameterList::addChild(const TSFParameterList& list)
{
	ptr_->addChild(list);
}

void TSFParameterList::setChild(const TSFParameterList& list)
{
	ptr_->setChild(list);
}

TSFParameterList TSFParameterList::getChild(const string& name) const
{
	TSFParameterList rtn;
	ptr_->getChild(name, rtn);
	return rtn;
}

/* ------------- adding and setting parameters ------------------------------*/

void TSFParameterList::addParameter(const TSFParameter& parameter)
{
	ptr_->addParameter(parameter);
}

void TSFParameterList::setParameter(const TSFParameter& parameter)
{
	ptr_->setParameter(parameter);
}

void TSFParameterList::setValue(const string& name, char value)
{
	ptr_->setValue(name, value);
}

void TSFParameterList::setValue(const string& name, const string& value)
{
	ptr_->setValue(name, value);
}

void TSFParameterList::setValue(const string& name, int value)
{
	ptr_->setValue(name, value);
}

void TSFParameterList::setValue(const string& name, double value)
{
	ptr_->setValue(name, value);
}

void TSFParameterList::setValue(const string& name, const TSFArray<int>& value)
{
	ptr_->setValue(name, value);
}

void TSFParameterList::setValue(const string& name, const TSFArray<double>& value)
{
	ptr_->setValue(name, value);
}

/* ------------- getting and listing parameters ------------------------------*/

TSFParameter TSFParameterList::getParameter(const string& name) const
{
	return ptr_->getParameter(name);
}

TSFArray<TSFParameter> TSFParameterList::listParameters() const 
{
	return ptr_->listParameters();
}

TSFArray<TSFParameterList> TSFParameterList::listChildren() const 
{
	return ptr_->listChildren();
}

/* ------------- lists of labels ------------------------------*/

TSFArray<string> TSFParameterList::listParameterNames() const 
{
	return ptr_->listParameterNames();
}

TSFArray<string> TSFParameterList::listChildNames() const 
{
	return ptr_->listChildNames();
}

void TSFParameterList::print(ostream& os) const 
{
  ptr_->print(os, 0);
}

TSFParameterList TSFParameterList::overrideWith(const TSFParameterList& newParams) const
{
  TSFParameterList rtn(*this);

  TSFArray<TSFParameter> p = newParams.listParameters();
  TSFArray<TSFParameterList> c = newParams.listChildren();

  for (int i=0; i<p.length(); i++)
    {
      rtn.setParameter(p[i]);
    }
  for (int i=0; i<c.length(); i++)
    {
      rtn.setChild(c[i]);
    }
  
  return rtn;
}


