#include "TSFObject.h"
#include "TSFArray.h"

using namespace TSF;
using std::string;

TSFObject::TSFObject()
{}

TSFObject::~TSFObject()
{}

void TSFObject::describe(ostream& os) const 
{
  print(os);
}

string TSFObject::spaces(int indentDepth) 
{
  static TSFArray<string> rtn;
  static int firstCall = true;

  if (firstCall)
    {
      rtn.append("");
      firstCall = false;
    }

  if (indentDepth >= rtn.length())
    {
      int oldLen = rtn.length();
      for (int i=oldLen; i<indentDepth+1; i++)
        {
          string tmp="";
          for (int s=0; s<i; s++)
            {
              tmp += " ";
            }
          rtn.append(tmp);
        }
    }

  return rtn[indentDepth];
}

void TSFObject::printIndented(ostream& os, int indentDepth) const 
{
  os << spaces(indentDepth);
  print(os);
}


string TSFObject::toString() const 
{
  std::ostringstream oss;
  
  print(oss);

  oss << ends;

  return oss.str();
}

string TSFObject::typeName() const 
{
  return typeid(*this).name();
}
