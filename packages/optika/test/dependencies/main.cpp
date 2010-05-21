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
#include "Teuchos_LocalTestingHelpers.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Optika_StandardDependencies.hpp"
#include "Optika_DependencySheet.hpp"

int intFuncTester(int argument){
	return argument+10;
}

int intVisualTester(int argument){
	if(argument <= 32){
		return 1;
	}
	else{
		return 0;
	}
}

double fondueTempTester(double argument){
	return argument-100;
}


int testValiDeps(Teuchos::FancyOStream &out){
	bool success = true;
	Teuchos::RCP<Teuchos::ParameterList> My_deplist = Teuchos::RCP<Teuchos::ParameterList>(new Teuchos::ParameterList);
	Teuchos::RCP<Optika::DependencySheet> depSheet1 = Teuchos::RCP<Optika::DependencySheet>(new Optika::DependencySheet(My_deplist));

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

	Optika::StringValidatorDependency::ValueToValidatorMap testValidatorMap1;
	testValidatorMap1["Cheese"] = cheeseValidator;
	testValidatorMap1["Soda"] = sodaValidator;
	testValidatorMap1["Chips"] = chipsValidator;

	Teuchos::ParameterList& stringValiDepList = My_deplist->sublist("String Validator Dependency", false, "String Validator Dependency testing list.\nWorking June 27th 2009");
	stringValiDepList.set("Food Selector", "Swiss", "select the food you want", cheeseValidator);
	stringValiDepList.set("Food Type", "Cheese", "String Validator Dependency Tester", stringFoodTypeValidator);

	Teuchos::RCP<Optika::StringValidatorDependency> 
	stringValiDep = Teuchos::RCP<Optika::StringValidatorDependency>(
		new Optika::StringValidatorDependency("Food Type", 
			Teuchos::sublist(My_deplist,"String Validator Dependency"),
			"Food Selector", 
			Teuchos::sublist(My_deplist,"String Validator Dependency"),
			testValidatorMap1, 
			cheeseValidator
		)
	);

	depSheet1->addDependency(stringValiDep);
	
	TEST_NOTHROW(stringValiDepList.validateParameters(stringValiDepList));
	TEST_ASSERT(depSheet1->hasDependents(stringValiDepList.getEntryPtr("Food Type")));
	Optika::DependencySheet::DepSet stringValiDepSet = depSheet1->getDependenciesForParameter(stringValiDepList.getEntryPtr("Food Type"));
	TEST_ASSERT(stringValiDepSet.size() == 1);
	stringValiDepList.set("Food Type","Soda");
	stringValiDep->evaluate();
	TEST_ASSERT(stringValiDepList.getEntry("Food Selector").validator().shares_resource(sodaValidator));
	TEST_THROW(stringValiDepList.validateParameters(stringValiDepList), Teuchos::Exceptions::InvalidParameterValue);
	stringValiDepList.set("Food Selector", "Pepsi");
	TEST_NOTHROW(stringValiDepList.validateParameters(stringValiDepList));
 

	Teuchos::ParameterList& 
	stringValiDepList2 = My_deplist->sublist(
		"String Validator Dependency (other validators)",
		false,
		"String Validator Dependency testing list for EnhancedNumber Validators.\nWorking June 27th 2009"
	);

	Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
	stringRangeValidator = Teuchos::rcp(
		new Teuchos::StringToIntegralParameterEntryValidator<int>(
		Teuchos::tuple<std::string>( "1-10", "10-33", "50-60" ),
		"Range selector"
		)
	);

	Teuchos::RCP<Optika::EnhancedNumberValidator<int> > range110Vali = 
	Teuchos::rcp(new Optika::EnhancedNumberValidator<int>(1,10));
	Teuchos::RCP<Optika::EnhancedNumberValidator<int> > range1033Vali = 
	Teuchos::rcp(new Optika::EnhancedNumberValidator<int>(10,33));
	Teuchos::RCP<Optika::EnhancedNumberValidator<int> > range5060Vali = 
	Teuchos::rcp(new Optika::EnhancedNumberValidator<int>(50,60));

	stringValiDepList2.set("Range selector", "1-10", "selects the range to validate", stringRangeValidator);

	Optika::StringValidatorDependency::ValueToValidatorMap rangeValidatorMap1;
	rangeValidatorMap1["1-10"] = range110Vali;
	rangeValidatorMap1["10-33"] = range1033Vali;
	rangeValidatorMap1["50-60"] = range5060Vali;
	stringValiDepList2.set("RangeValue", 3, "the value of the range", range110Vali);

	Teuchos::RCP<Optika::StringValidatorDependency> 
	stringValiDep2 = Teuchos::RCP<Optika::StringValidatorDependency>(
		new Optika::StringValidatorDependency(
			"Range selector", 
			Teuchos::sublist(My_deplist,"String Validator Dependency (other validators)"),
			"RangeValue", 
			Teuchos::sublist(My_deplist,"String Validator Dependency (other validators)"),
			rangeValidatorMap1, 
			range110Vali
		)
	);

	depSheet1->addDependency(stringValiDep2);

	TEST_NOTHROW(stringValiDepList2.validateParameters(stringValiDepList2));
	TEST_ASSERT(depSheet1->hasDependents(stringValiDepList2.getEntryPtr("Range selector")));
	Optika::DependencySheet::DepSet stringValiDepSet2 = depSheet1->getDependenciesForParameter(stringValiDepList2.getEntryPtr("Range selector"));
	TEST_ASSERT(stringValiDepSet2.size() == 1);
	stringValiDepList2.set("Range selector","50-60");
	stringValiDep2->evaluate();
	TEST_ASSERT(stringValiDepList2.getEntry("RangeValue").validator().shares_resource(range5060Vali));
	TEST_THROW(stringValiDepList2.validateParameters(stringValiDepList2), Teuchos::Exceptions::InvalidParameterValue);
	stringValiDepList2.set("RangeValue", 55);
	TEST_NOTHROW(stringValiDepList2.validateParameters(stringValiDepList2));


	Teuchos::ParameterList&
	boolValidatorDepList = My_deplist->sublist(
		"Bool Validator Dependency List", 
		false,
		"Bool Validator Dependency testing list.\nWorking June 27th 2009"
	);

	boolValidatorDepList.set("Use Validator?", true, "truns the validator on and off");
	Teuchos::RCP<Optika::EnhancedNumberValidator<int> > basicVali = Teuchos::rcp(new Optika::EnhancedNumberValidator<int>(1,10));
	boolValidatorDepList.set("do I have a validator?", 4, "does it have a validator?", basicVali);

	Teuchos::RCP<Optika::BoolValidatorDependency> 
	boolValiDep = Teuchos::RCP<Optika::BoolValidatorDependency>(
		new Optika::BoolValidatorDependency(
			"Use Validator?", 
			Teuchos::sublist(My_deplist,"Bool Validator Dependency List"),
			"do I have a validator?", 
			Teuchos::sublist(My_deplist,"Bool Validator Dependency List"),
			basicVali, 
			Teuchos::RCP<Teuchos::ParameterEntryValidator>()
		)
	);

	depSheet1->addDependency(boolValiDep);

	TEST_ASSERT(depSheet1->hasDependents(boolValidatorDepList.getEntryPtr("Use Validator?")));
	TEST_ASSERT(boolValidatorDepList.getEntry("do I have a validator?").validator().shares_resource(basicVali));
	TEST_NOTHROW(boolValidatorDepList.validateParameters(boolValidatorDepList));
	Optika::DependencySheet::DepSet boolValiDepSet = depSheet1->getDependenciesForParameter(boolValidatorDepList.getEntryPtr("Use Validator?"));
	TEST_ASSERT(boolValiDepSet.size() == 1);
	boolValidatorDepList.set("Use Validator?",false);
	boolValiDep->evaluate();
	TEST_ASSERT(boolValidatorDepList.getEntry("do I have a validator?").validator().is_null());


	Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
	lowTempCheeseValidator = Teuchos::rcp(
		new Teuchos::StringToIntegralParameterEntryValidator<int>(
			Teuchos::tuple<std::string>( "PepperJack", "Swiss", "American" ),
			"Cheese to Fondue"
		)
	);

	Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
	highTempCheeseValidator = Teuchos::rcp(
		new Teuchos::StringToIntegralParameterEntryValidator<int>(
			Teuchos::tuple<std::string>( "Munster", "Provalone", "Kurtis Super Awesome Cheese"),
			"Cheese to Fondue"
		)
	);

	Teuchos::ParameterList& 
	rangeValidatorDepList = My_deplist->sublist(
		"Range Validator and NumberVisual Dependency List",
		false,
		"Range Validator and Number Visual Dependency testing list.\nWorking June 27th 2009"
	);
	rangeValidatorDepList.set("Temperature",101.0, "The temperature of the fondue");
	rangeValidatorDepList.set("Cheese to Fondue", "Swiss", "The cheese we'll be using in our fondue pot.", lowTempCheeseValidator);
	Optika::RangeValidatorDependency<double>::RangeToValidatorMap tempranges;
	tempranges[std::pair<double,double>(100,200)] = lowTempCheeseValidator;
	tempranges[std::pair<double,double>(200,300)] = highTempCheeseValidator;
	Teuchos::RCP<Optika::RangeValidatorDependency<double> > 
	cheeseTempDep = Teuchos::RCP<Optika::RangeValidatorDependency<double> >(
		new Optika::RangeValidatorDependency<double>(
			"Temperature", 
			Teuchos::sublist(My_deplist,"Range Validator and NumberVisual Dependency List"),
			"Cheese to Fondue", 
			Teuchos::sublist(My_deplist,"Range Validator and NumberVisual Dependency List"),
			tempranges, 
			lowTempCheeseValidator
		)
	);
	depSheet1->addDependency(cheeseTempDep);

	TEST_ASSERT(depSheet1->hasDependents(rangeValidatorDepList.getEntryPtr("Temperature")));
	Optika::DependencySheet::DepSet rangeValiDepSet = depSheet1->getDependenciesForParameter(rangeValidatorDepList.getEntryPtr("Temperature"));
	TEST_ASSERT(rangeValiDepSet.size() == 1);
	rangeValidatorDepList.set("Temperature",250.0);
	cheeseTempDep->evaluate();
	TEST_ASSERT(rangeValidatorDepList.getEntry("Cheese to Fondue").validator().shares_resource(highTempCheeseValidator));
	TEST_THROW(rangeValidatorDepList.validateParameters(rangeValidatorDepList), Teuchos::Exceptions::InvalidParameterValue);
	rangeValidatorDepList.set("Cheese to Fondue", "Provalone");
	TEST_NOTHROW(rangeValidatorDepList.validateParameters(rangeValidatorDepList));

	return (success ? 0:1);
}
  
int testVisualDeps(Teuchos::FancyOStream &out){
	bool success = false;
	Teuchos::RCP<Teuchos::ParameterList> My_deplist = Teuchos::RCP<Teuchos::ParameterList>(new Teuchos::ParameterList);
	Teuchos::RCP<Optika::DependencySheet> depSheet1 = Teuchos::RCP<Optika::DependencySheet>(new Optika::DependencySheet(My_deplist));
	double (*fondueFunc)(double);
	fondueFunc = fondueTempTester;

	Teuchos::RCP<Optika::NumberVisualDependency<double> > fondueDep = 
	Teuchos::RCP<Optika::NumberVisualDependency<double> >(
		new Optika::NumberVisualDependency<double>(
			"Temperature", 
			Teuchos::sublist(My_deplist,"Range Validator and NumberVisual Dependency List"),
			"Cheese to Fondue", 
			Teuchos::sublist(My_deplist,"Range Validator and NumberVisual Dependency List"),
			fondueFunc
		)
	);
	depSheet1->addDependency(fondueDep);




  Teuchos::ParameterList&
    numberArrayLengthDepList = My_deplist->sublist("Number Array Length Dependency List", false, "Number Array Length ependecy testing list.\nWorking June 27th 2009");
  numberArrayLengthDepList.set("Array Length", 8, "array length setter");
  Teuchos::Array<double> variableLengthArray(10,23.0);
  Teuchos::RCP<Optika::EnhancedNumberValidator<double> > varLengthArrayVali = Teuchos::RCP<Optika::EnhancedNumberValidator<double> >(
  	new Optika::EnhancedNumberValidator<double>(10,50,4) );
  numberArrayLengthDepList.set("Variable Length Array", variableLengthArray, "variable length array",
  Teuchos::RCP<Optika::ArrayNumberValidator<double> >(new Optika::ArrayNumberValidator<double>(varLengthArrayVali)));

  Teuchos::RCP<Optika::NumberArrayLengthDependency> arrayLengthDep = Teuchos::RCP<Optika::NumberArrayLengthDependency>(
  	new Optika::NumberArrayLengthDependency(
	"Array Length", 
	Teuchos::sublist(My_deplist,"Number Array Length Dependency List"),
	"Variable Length Array", 
	Teuchos::sublist(My_deplist,"Number Array Length Dependency List")));
  depSheet1->addDependency(arrayLengthDep);



  Teuchos::ParameterList&
    numberValiAspDepList = My_deplist->sublist("Number Validator Aspect Dependency List", false, "Number Validator Aspect Dependency testing list.\nWorking June 27th 2009");
  Teuchos::RCP<Optika::EnhancedNumberValidator<int> > intVali2 = 
  	Teuchos::rcp(new Optika::EnhancedNumberValidator<int>(0,20));
  numberValiAspDepList.set("Int", 8, "Int tester", intVali2);
  numberValiAspDepList.set("Int2", 8, "int2 tester", intVali2);
  numberValiAspDepList.set("Int dependee", 1, "Int dependee");

  int (*func)(int);
  func = intFuncTester;

  Teuchos::RCP<Optika::NumberValidatorAspectDependency<int> > intDep1 =
     Teuchos::RCP<Optika::NumberValidatorAspectDependency<int> >(
       new Optika::NumberValidatorAspectDependency<int>(
	 "Int dependee",
	Teuchos::sublist(My_deplist,"Number Validator Aspect Dependency List"),
         "Int",
	Teuchos::sublist(My_deplist,"Number Validator Aspect Dependency List"),
	 intVali2,
	 Optika::NumberValidatorAspectDependency<int>::Max,
	 func
	 ));
  Teuchos::RCP<Optika::NumberValidatorAspectDependency<int> > intDep2 =
     Teuchos::RCP<Optika::NumberValidatorAspectDependency<int> >(
       new Optika::NumberValidatorAspectDependency<int>(
	 "Int dependee",
	Teuchos::sublist(My_deplist,"Number Validator Aspect Dependency List"),
         "Int2",
	Teuchos::sublist(My_deplist,"Number Validator Aspect Dependency List"),
	 intVali2,
	 Optika::NumberValidatorAspectDependency<int>::Max,
	 func
	 ));

  depSheet1->addDependency(intDep1);
  depSheet1->addDependency(intDep2);

  



  Teuchos::ParameterList&
    boolVisDepList = My_deplist->sublist("Bool Visual Dependency List", false, "Bool Visual Dependency testing list.\nWorking June 29 2009");
  boolVisDepList.set("ShowPrecs", true, "Whether or not to should the Preciondtioner list");
  Teuchos::ParameterList&
    Prec_List0 = boolVisDepList.sublist("Preconditioner",false,"Sublist that defines the preconditioner.");
  Prec_List0.set("Type", "ILU", "The tpye of preconditioner to use");
  Teuchos::RCP<Optika::EnhancedNumberValidator<double> > droptolValidator = Teuchos::rcp(new Optika::EnhancedNumberValidator<double>(0,10,1e-3));
  Prec_List0.set("Drop Tolerance", 1e-3
                ,"The tolerance below which entries from the\n""factorization are left out of the factors.", droptolValidator);
  Teuchos::RCP<Optika::BoolVisualDependency> precDep1 = Teuchos::RCP<Optika::BoolVisualDependency>(new Optika::BoolVisualDependency(
  "ShowPrecs", 
	Teuchos::sublist(My_deplist,"Bool Visual Dependency List"),
  "Preconditioner", 
	Teuchos::sublist(My_deplist,"Bool Visual Dependency List"),
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
   Teuchos::RCP<Optika::EnhancedNumberValidator<int> > swissValidator = Teuchos::rcp(new Optika::EnhancedNumberValidator<int>(0,10));
   stringVisDepList.set("Swiss rating", 0, "How you rate swiss on a scale of 1 to 10", swissValidator);
   Teuchos::RCP<Optika::StringVisualDependency> swissDep1 = 
      Teuchos::RCP<Optika::StringVisualDependency>(new Optika::StringVisualDependency(
      "Favorite Cheese", 
	Teuchos::sublist(My_deplist,"String Visual Dependency List"),
      "Swiss rating", 
	Teuchos::sublist(My_deplist,"String Visual Dependency List"),
      "Swiss", 
      true));
   depSheet1->addDependency(swissDep1);





  int (*visfunc)(int);
  visfunc = intVisualTester;
Teuchos::ParameterList&
    numberVisDepList = My_deplist->sublist("Number Visual Dependency List", false, "Number Visual Dependency testing list.\nWorking June 27th 2009");
  numberVisDepList.set("Ice", 50, "Ice stuff");
  numberVisDepList.set("Room Temp", 10, "Room temperature");
  Teuchos::RCP<Optika::NumberVisualDependency<int> > iceDep = 
      Teuchos::RCP<Optika::NumberVisualDependency<int> >(new Optika::NumberVisualDependency<int>(
      "Room Temp", 
	Teuchos::sublist(My_deplist,"Number Visual Dependency List"),
      "Ice", 
	Teuchos::sublist(My_deplist,"Number Visual Dependency List"),
      visfunc));
  depSheet1->addDependency(iceDep);






  return (success ? 0:1);
}

int main(int argc, char* argv[]){
	bool success = true;
	Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
	if(testValiDeps(*out) == 1){
		success = false;
	}
/*	if(testVisualDeps(*out) ==1){
		success = false;
	}*/
	return (success ? 0:1);
}

