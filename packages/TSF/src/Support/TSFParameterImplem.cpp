#include "TSFParameterImplem.h"
#include "TSFError.h"

using namespace TSF;
using namespace std;


//========================================= ctors =================================

TSFParameterImplem::TSFParameterImplem(const string& label) 
{
  setDefaults(label);
}

TSFParameterImplem::TSFParameterImplem(const string& label, char parameter)  
{
	setDefaults(label);
  set(parameter);
}

TSFParameterImplem::TSFParameterImplem(const string& label, const string& parameter)  
{
	setDefaults(label);
  set(parameter);
}

TSFParameterImplem::TSFParameterImplem(const string& label, int parameter)  
{
	setDefaults(label);
  set(parameter);
}

TSFParameterImplem::TSFParameterImplem(const string& label, const TSFArray<int>& parameter)  
{
	setDefaults(label);
  set(parameter);
}

TSFParameterImplem::TSFParameterImplem(const string& label, double parameter)  
{
	setDefaults(label);
  set(parameter);
}

TSFParameterImplem::TSFParameterImplem(const string& label, 
																			 const TSFArray<double>& parameter)  
{
	setDefaults(label);
  set(parameter);
}

//========== set methods ==========================================

void TSFParameterImplem::setDefaults(const string& label)
{
	label_ = label;

	typeIsDefined_ = false;

	isChar_ = false;
	isCharString_ = false;
	isInt_ = false;
	isIntArray_ = false;
	isDouble_ = false;
	isDoubleArray_ = false;
}

void TSFParameterImplem::set(char parameter)  
{
	if (typeIsDefined_ && !isChar()) typeChangeError("char");

  charValue_ = parameter;
  isChar_ = true;
	typeIsDefined_ = true;
}

void TSFParameterImplem::set(const string& parameter)  
{
	if (typeIsDefined_ && !isCharString()) typeChangeError("string");

	charStringValue_ = parameter;
  isCharString_ = true;
	typeIsDefined_ = true;
}

void TSFParameterImplem::set(int parameter)  
{
	if (typeIsDefined_ && !isInt()) typeChangeError("int");

  intValue_ = parameter;
  isInt_ = true;
	typeIsDefined_ = true;
}

void TSFParameterImplem::set(const TSFArray<int>& parameter)
{
	if (typeIsDefined_ && !isIntArray()) typeChangeError("int array");

	intArrayValue_ = parameter;
  isIntArray_ = true;
	typeIsDefined_ = true;
}

void TSFParameterImplem::set(double parameter)  
{
	if (typeIsDefined_ && !isDouble()) typeChangeError("double");

  doubleValue_ = parameter;
  isDouble_ = true;
	typeIsDefined_ = true;
}

void TSFParameterImplem::set(const TSFArray<double>& parameter)
{
	if (typeIsDefined_ && !isDoubleArray()) typeChangeError("double array");

	doubleArrayValue_ = parameter;
  isDoubleArray_ = true;
	typeIsDefined_ = true;
}

//========== internal utilities ==========================================

void TSFParameterImplem::wrongType(const string& wrongType) const
{
	TSFError::raise("tried to access " + type() + "-valued parameter "
									+ getLabel() + " as if it were a " + wrongType);
}

void TSFParameterImplem::typeChangeError(const string& attemptedType) const 
{
	TSFError::raise("tried to change parameter " + getLabel() + " from type " + type()
									+ " to type " + attemptedType);
}

string TSFParameterImplem::type() const 
{
	 // Call each get function.  One of them will eventually work!
  if (isChar()) return "char";
  if (isCharString()) return "char string";
  if (isInt()) return "int";
  if (isIntArray()) return "int array";
  if (isDouble()) return "double";
  if (isDoubleArray()) return "double array";

	TSFError::raise("parameter " + getLabel() + " has no type defined");
	return "undefined"; 
}

