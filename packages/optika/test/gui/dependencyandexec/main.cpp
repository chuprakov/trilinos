// @HEADER
// ***********************************************************************
// 
//         Optika: A Tool For Developing Parameter Obtaining GUIs
//                Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, with Sandia Corporation, the 
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Kurtis Nusbaum (klnusbaum@gmail.com) 
// 
// ***********************************************************************
// @HEADER
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Optika_GUI.hpp"
#include "Teuchos_StandardDependencies.hpp"

void print(Teuchos::RCP<const Teuchos::ParameterList> theList){
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::writeParameterListToXmlOStream(*theList, *out);
}

class IntVisualTester : public Teuchos::SingleArguementFunctionObject<int, int>{

public:
  int runFunction() const;
};

int IntVisualTester::runFunction() const{
	if(getParameterValue() <= 32){
		return 1;
	}
	else{
		return 0;
	}
}

class IntFuncTester : public Teuchos::SingleArguementFunctionObject<int, int>{

public:
  int runFunction() const;
};

int IntFuncTester::runFunction() const{
  return getParameterValue() + 10;
}

class FondueTempTester : public Teuchos::SingleArguementFunctionObject<double, double>{

public:
  double runFunction() const;
};

double FondueTempTester::runFunction() const{
	return getParameterValue()-100;
}

int main(){
 Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
 Teuchos::RCP<Teuchos::ParameterList> My_deplist = rcp(new Teuchos::ParameterList);
 Teuchos::RCP<Teuchos::DependencySheet> depSheet1 = rcp(new Teuchos::DependencySheet);

  Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
    stringFoodTypeValidator = Teuchos::rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<int>(
        Teuchos::tuple<std::string>( "Cheese", "Soda", "Chips" )
        ,"Food Type"
        )
      );

  Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
    cheeseValidator = Teuchos::rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<int>(
        Teuchos::tuple<std::string>( "Swiss", "American", "Super Awesome Cheese" )
        ,"Food Selector"
        )
      );
  Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
    sodaValidator = Teuchos::rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<int>(
        Teuchos::tuple<std::string>( "Pepsi", "Coke", "Kurtis Cola", "Bad Cola" )
        ,"Food Selector"
        )
      );
  Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
    chipsValidator = Teuchos::rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<int>(
        Teuchos::tuple<std::string>( "Lays", "Doritos", "Kurtis Super Awesome Brand" )
        ,"Food Selector"
        )
      );

  Teuchos::StringValidatorDependency::ValueToValidatorMap testValidatorMap1;
  testValidatorMap1["Cheese"] = cheeseValidator;
  testValidatorMap1["Soda"] = sodaValidator;
  testValidatorMap1["Chips"] = chipsValidator;


  Teuchos::ParameterList&
    stringValiDepList = My_deplist->sublist("String Validator Dependency", false, "String Validator Dependency testing list.\nWorking June 27th 2009");
  stringValiDepList.set("Food Selector", "Swiss", "select the food you want", cheeseValidator);
  stringValiDepList.set("Food Type", "Cheese", "String Validator Dependency Tester", stringFoodTypeValidator);
  depSheet1->addDependency(Teuchos::RCP<Teuchos::StringValidatorDependency>(
  	new Teuchos::StringValidatorDependency(
    My_deplist->getEntryRCP("Food Type"),
	  My_deplist->getEntryRCP("Food Selector"),
	  testValidatorMap1 )));
 

  Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
    stringRangeValidator = Teuchos::rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<int>(
        Teuchos::tuple<std::string>( "1-10", "10-33", "50-60" )
        ,"Range selector"
        )
      );
  Teuchos::ParameterList&
    stringValiDepList2 = My_deplist->sublist("String Validator Dependency (other validators)", false, "String Validator Dependency testing list for EnhancedNumber Validators.");
  stringValiDepList2.set("Range selector", "1-10", "selects the range to validate", stringRangeValidator);
  Teuchos::RCP<Teuchos::EnhancedNumberValidator<int> > range110Vali = 
  	Teuchos::rcp(new Teuchos::EnhancedNumberValidator<int>(1,10));
  Teuchos::RCP<Teuchos::EnhancedNumberValidator<int> > range1033Vali = 
  	Teuchos::rcp(new Teuchos::EnhancedNumberValidator<int>(10,33));
  Teuchos::RCP<Teuchos::EnhancedNumberValidator<int> > range5060Vali = 
  	Teuchos::rcp(new Teuchos::EnhancedNumberValidator<int>(50,60));
  Teuchos::StringValidatorDependency::ValueToValidatorMap rangeValidatorMap1;
  rangeValidatorMap1["1-10"] = range110Vali;
  rangeValidatorMap1["10-33"] = range1033Vali;
  rangeValidatorMap1["50-60"] = range5060Vali;
  stringValiDepList2.set("RangeValue", 3, "the value of the range", range110Vali);
  depSheet1->addDependency(Teuchos::RCP<Teuchos::StringValidatorDependency>(
  	new Teuchos::StringValidatorDependency(
	My_deplist->getEntryRCP("Range selector"),
	My_deplist->getEntryRCP("RangeValue"),
	rangeValidatorMap1)));


  Teuchos::ParameterList&
    boolValidatorDepList = My_deplist->sublist("Bool Validator Dependency List", false, "Bool Validator Dependency testing list.\nWorking June 27th 2009");
  boolValidatorDepList.set("Use Validator?", true, "truns the validator on and off");
  Teuchos::RCP<Teuchos::EnhancedNumberValidator<int> > basicVali = 
  	Teuchos::rcp(new Teuchos::EnhancedNumberValidator<int>(1,10));
  boolValidatorDepList.set("do I have a validator?", 4, "does it have a validator?", basicVali);
  depSheet1->addDependency(Teuchos::RCP<Teuchos::BoolValidatorDependency>(
  	new Teuchos::BoolValidatorDependency(
	    My_deplist->getEntryRCP("Use Validator?"),
	    My_deplist->getEntryRCP("do I have a validator?"),
	    basicVali, 
	    Teuchos::RCP<Teuchos::ParameterEntryValidator>())));


  Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
    lowTempCheeseValidator = Teuchos::rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<int>(
        Teuchos::tuple<std::string>( "PepperJack", "Swiss", "American" )
        ,"Cheese to Fondue"
        )
      );

  Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
    highTempCheeseValidator = Teuchos::rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<int>(
        Teuchos::tuple<std::string>( "Munster", "Provalone", "Kurtis Super Awesome Cheese")
        ,"Cheese to Fondue"
        )
      );
  Teuchos::ParameterList&
    rangeValidatorDepList = My_deplist->sublist("Range Validator and NumberVisual Dependency List", false, "Range Validator and Number Visual Dependency testing list.");
  rangeValidatorDepList.set("Temperature",101.0, "The temperature of the fondue");
  rangeValidatorDepList.set("Cheese to Fondue", "Swiss", "The cheese we'll be using in our fondue pot.", lowTempCheeseValidator);
  Teuchos::RangeValidatorDependency<double>::RangeToValidatorMap tempranges;
  tempranges[std::pair<double,double>(100,200)] = lowTempCheeseValidator;
  tempranges[std::pair<double,double>(200,300)] = highTempCheeseValidator;
  Teuchos::RCP<Teuchos::RangeValidatorDependency<double> > cheeseTempDep = Teuchos::rcp(
  	new Teuchos::RangeValidatorDependency<double>(
	    My_deplist->getEntryRCP("Temperature"),
	    My_deplist->getEntryRCP("Cheese to Fondue"),
	    tempranges));
 
  depSheet1->addDependency(cheeseTempDep);
  
  Teuchos::RCP<FondueTempTester> fondueFunc = Teuchos::rcp(new FondueTempTester());
  Teuchos::RCP<Teuchos::NumberVisualDependency<double> > fondueDep = 
      Teuchos::RCP<Teuchos::NumberVisualDependency<double> >(new Teuchos::NumberVisualDependency<double>(
      My_deplist->getEntryRCP("Temperature"),
      My_deplist->getEntryRCP("Cheese to Fondue"),
      fondueFunc));
  depSheet1->addDependency(fondueDep);

  Teuchos::ParameterList&
    numberArrayLengthDepList = My_deplist->sublist("Number Array Length Dependency List", false, "Number Array Length ependecy testing list.");
  numberArrayLengthDepList.set("Array Length", 8, "array length setter");
  Teuchos::Array<double> variableLengthArray(10,23.0);
  Teuchos::RCP<Teuchos::EnhancedNumberValidator<double> > varLengthArrayVali = Teuchos::RCP<Teuchos::EnhancedNumberValidator<double> >(
  	new Teuchos::EnhancedNumberValidator<double>(10,50,4) );
  numberArrayLengthDepList.set("Variable Length Array", variableLengthArray, "variable length array",
  Teuchos::RCP<Teuchos::ArrayNumberValidator<double> >(new Teuchos::ArrayNumberValidator<double>(varLengthArrayVali)));

  Teuchos::RCP<Teuchos::NumberArrayLengthDependency<int, double> > arrayLengthDep(
    new Teuchos::NumberArrayLengthDependency<int, double>(My_deplist->getEntryRCP("Array Length"), 
    My_deplist->getEntryRCP("Variable Length Array")));
  depSheet1->addDependency(arrayLengthDep);



  Teuchos::ParameterList&
    numberValiAspDepList = My_deplist->sublist("Number Validator Aspect Dependency List", false, "Number Validator Aspect Dependency testing list.");
  Teuchos::RCP<Teuchos::EnhancedNumberValidator<int> > intVali2 = 
  	Teuchos::rcp(new Teuchos::EnhancedNumberValidator<int>(0,20));
  numberValiAspDepList.set("Int", 8, "Int tester", intVali2);
  numberValiAspDepList.set("Int2", 8, "int2 tester", intVali2);
  numberValiAspDepList.set("Int dependee", 1, "Int dependee");

  Teuchos::ParameterList&
    boolVisDepList = My_deplist->sublist("Bool Visual Dependency List", false, "Bool Visual Dependency testing list.");
  boolVisDepList.set("ShowPrecs", true, "Whether or not to should the Preciondtioner list");
  Teuchos::ParameterList&
    Prec_List0 = boolVisDepList.sublist("Preconditioner",false,"Sublist that defines the preconditioner.");
  Prec_List0.set("Type", "ILU", "The tpye of preconditioner to use");
  Teuchos::RCP<Teuchos::EnhancedNumberValidator<double> > droptolValidator = Teuchos::rcp(new Teuchos::EnhancedNumberValidator<double>(0,10,1e-3));
  Prec_List0.set("Drop Tolerance", 1e-3
                ,"The tolerance below which entries from the\n""factorization are left out of the factors.", droptolValidator);
  Teuchos::RCP<Teuchos::BoolVisualDependency> precDep1 = Teuchos::rcp(new Teuchos::BoolVisualDependency(
  My_deplist->getEntryRCP("ShowPrecs"),
  My_deplist->getEntryRCP("Preconditioner"),
  true));
  depSheet1->addDependency(precDep1);




Teuchos::ParameterList&
    stringVisDepList = My_deplist->sublist("String Visual Dependency List", false, "String Visual Dependency testing list.\nWorking June 29 2009");
  Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
    favCheeseValidator = Teuchos::rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<int>(
        Teuchos::tuple<std::string>( "Swiss", "American", "Cheder" )
        ,"Favorite Cheese"
        )
      );
   
   stringVisDepList.set(
   	"Favorite Cheese", "American", "Your favorite type of cheese", favCheeseValidator);
   Teuchos::RCP<Teuchos::EnhancedNumberValidator<int> > swissValidator = Teuchos::rcp(new Teuchos::EnhancedNumberValidator<int>(0,10));
   stringVisDepList.set("Swiss rating", 0, "How you rate swiss on a scale of 1 to 10", swissValidator);
   Teuchos::RCP<Teuchos::StringVisualDependency> swissDep1 = 
      Teuchos::RCP<Teuchos::StringVisualDependency>(new Teuchos::StringVisualDependency(
      My_deplist->getEntryRCP("Favorite Cheese"),
      My_deplist->getEntryRCP("Swiss rating"),
      "Swiss", 
      true));
   depSheet1->addDependency(swissDep1);





  Teuchos::RCP<IntVisualTester> intVisTester = 
    Teuchos::rcp(new IntVisualTester());
Teuchos::ParameterList&
    numberVisDepList = My_deplist->sublist("Number Visual Dependency List", false, "Number Visual Dependency testing list.");
  numberVisDepList.set("Ice", 50, "Ice stuff");
  numberVisDepList.set("Room Temp", 10, "Room temperature");
  Teuchos::RCP<Teuchos::NumberVisualDependency<int> > iceDep = 
      Teuchos::RCP<Teuchos::NumberVisualDependency<int> >(new Teuchos::NumberVisualDependency<int>(
      My_deplist->getEntryRCP("Room Temp"),
      My_deplist->getEntryRCP("Ice"),
      intVisTester));
  depSheet1->addDependency(iceDep);




 void (*myFunc)(Teuchos::RCP<const Teuchos::ParameterList>);
  myFunc = print;
  Optika::getInput(My_deplist, depSheet1, myFunc);

  std::cout << "Dep List: \n";
  Teuchos::writeParameterListToXmlOStream(*My_deplist, *out);

  std::cout << "Deps: \n";
  depSheet1->printDeps();


  return 0;
}


