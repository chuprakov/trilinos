#include "TSFParameterBase.h"
#include "TSFError.h"

using namespace TSF;
using namespace std;



TSFParameterBase::TSFParameterBase(const string& label)
	: label_(label)
{}


void TSFParameterBase::set(char /* value */) 
{
	typeChangeError("char");
}

void TSFParameterBase::set(const string& /* value */) 
{
	typeChangeError("char string");
}

void TSFParameterBase::set(int /* value */) 
{
	typeChangeError("int");
}

void TSFParameterBase::set(const TSFArray<int>& /* value */) 
{
	typeChangeError("int array");
}

void TSFParameterBase::set(double /* value */) 
{
	typeChangeError("double");
}

void TSFParameterBase::set(const TSFArray<double>& /* value */) 
{
	typeChangeError("double array");
}

char TSFParameterBase::getChar() const 
{
	typeAccessError("char");
	return 0 ; // -Wall
}

string TSFParameterBase::getCharString() const 
{
	typeAccessError("char string");
	return ""; // -Wall
}

int TSFParameterBase::getInt() const 
{
	typeAccessError("int");
	return 0 ; // -Wall
}

double TSFParameterBase::getDouble() const 
{
	typeAccessError("double");
	return 0.0; // -Wall
}

TSFArray<int> TSFParameterBase::getIntArray() const
{
	typeAccessError("int array");
	return TSFArray<int>(0); // -Wall
}

TSFArray<double> TSFParameterBase::getDoubleArray() const
{
	typeAccessError("double array");
	return TSFArray<double>(0); // -Wall
}







void TSFParameterBase::typeAccessError(const string& wrongType) const
{
	TSFError::raise("tried to access a " + type() + "-valued parameter "
									+ getLabel() + " as if it were a " + wrongType);
}

void TSFParameterBase::typeChangeError(const string& wrongType) const 
{
	TSFError::raise("tried to change parameter [" + getLabel() + "] from type " + type()
									+ " to type " + wrongType);
}

/* --------------- ctors for derived types ------------------- */

TSFCharParameter::TSFCharParameter(const string& label, char value)
	: TSFParameterBase(label), value_(value)
{}

TSFCharStringParameter::TSFCharStringParameter(const string& label, 
																							 const string& value)
	: TSFParameterBase(label), value_(value)
{}

TSFIntParameter::TSFIntParameter(const string& label, int value)
	: TSFParameterBase(label), value_(value)
{}

TSFDoubleParameter::TSFDoubleParameter(const string& label, 
																			 const double& value)
	: TSFParameterBase(label), value_(value)
{}

TSFIntArrayParameter::TSFIntArrayParameter(const string& label, 
																			 const TSFArray<int>& value)
	: TSFParameterBase(label), value_(value)
{}

TSFDoubleArrayParameter::TSFDoubleArrayParameter(const string& label, 
																								 const TSFArray<double>& value)
	: TSFParameterBase(label), value_(value)
{}

