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

const string& TSFParameterList::getLabel() const
{
	return ptr_->getLabel();
}

TSFParameterList TSFParameterList::deepCopy() const 
{
	return ptr_->deepCopy();
}

/* ------------- setting and getting children  ------------------------------*/

void TSFParameterList::addChild(const TSFParameterList& list)
{
	ptr_->addChild(list);
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


