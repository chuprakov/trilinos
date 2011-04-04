// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardConditions.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_StandardConditions.hpp"

namespace Teuchos{

class DoubleTesterFunc : public SingleArguementFunctionObject<double, double>{
public:
  double runFunction() const{
    return getArguementValue() - 100.0;
  }
};

/**
 * Test all the conditions
 */
TEUCHOS_UNIT_TEST(Teuchos_Conditions, testConditions){
	//Settin up initial list
	RCP<ParameterList> testingList = rcp(new ParameterList("Condition Testing List"));

	/*
	 * Testing for string condition
	 */
	Array<std::string> validValues(tuple<std::string>("mountain dew", "pepsi", "coke", "fanta"));
	RCP<StringValidator> stringVali1 = rcp(new StringValidator(validValues));

	testingList->set("string param", "fanta", "parameter for testing string conditions", stringVali1);

	StringCondition::ValueList conValues1(tuple<std::string>("pepsi", "coke"));
	RCP<StringCondition> stringCon1 = rcp( new StringCondition(testingList->getEntryRCP("string param"), conValues1));
	TEST_ASSERT(!stringCon1->isConditionTrue());
	testingList->set("string param", "coke");
	TEST_ASSERT(stringCon1->isConditionTrue());
	RCP<StringCondition> stringCon2 = rcp( 
    new StringCondition(
      testingList->getEntryRCP("string param"), conValues1, false));
	testingList->set("string param", "fanta");
	TEST_ASSERT(stringCon2->isConditionTrue());
	testingList->set("string param", "coke");
	TEST_ASSERT(!stringCon2->isConditionTrue());

	/*
	 * Testing for number condition
	 */
	testingList->set("double param", 5.0, "parameter for testing number conditions");

	RCP<NumberCondition<double> > numberCon1 = 
    rcp( new NumberCondition<double>(testingList->getEntryRCP("double param"), null, true));
	TEST_ASSERT(numberCon1->isConditionTrue());
	testingList->set("double param", -1.0);
	TEST_ASSERT(!numberCon1->isConditionTrue());

  RCP<DoubleTesterFunc> doubleTesterFunc = rcp( new DoubleTesterFunc());
	RCP<NumberCondition<double> > numberCon2 = 
    rcp( new NumberCondition<double>(testingList->getEntryRCP("double param"), doubleTesterFunc, false));
	TEST_ASSERT(numberCon2->isConditionTrue());
	testingList->set("double param", 101.0);
	TEST_ASSERT(!numberCon2->isConditionTrue());

	/*
	 * Testing bool conditions
	 */
	testingList->set("bool param", true, "parameter for testing bool conditions");

	RCP<BoolCondition> boolCon1 = rcp( new BoolCondition(testingList->getEntryRCP("bool param")));
	TEST_ASSERT(boolCon1->isConditionTrue());
	testingList->set("bool param", false);
	TEST_ASSERT(!boolCon1->isConditionTrue());

	RCP<BoolCondition> boolCon2 = rcp( new BoolCondition(testingList->getEntryRCP("bool param"), false));
	TEST_ASSERT(boolCon2->isConditionTrue());
	testingList->set("bool param", true);
	TEST_ASSERT(!boolCon2->isConditionTrue());

	/*
	 * Test Not condition
	 */
	RCP<NotCondition> notCon1 = rcp(new NotCondition(numberCon1));
	TEST_ASSERT(!notCon1->isConditionTrue());
	testingList->set("double param", -1.0);
	TEST_ASSERT(notCon1->isConditionTrue());

	/*
	 * Test And condition
	 */
	Condition::ConstConditionList conList1(tuple<RCP<const Condition> >(stringCon1, boolCon1));
	RCP<AndCondition> andCon1 = rcp(new AndCondition(conList1));
	TEST_ASSERT(andCon1->isConditionTrue());
	Condition::ConstConditionList conList2(tuple<RCP<const Condition> >(stringCon1, boolCon2));
	RCP<AndCondition> andCon2 = rcp(new AndCondition(conList2));
	TEST_ASSERT(!andCon2->isConditionTrue());
	Condition::ConstConditionList conList3(tuple<RCP<const Condition> >(stringCon2, boolCon2));
	RCP<AndCondition> andCon3 = rcp(new AndCondition(conList3));
	TEST_ASSERT(!andCon3->isConditionTrue());

	/*
	 * Testing or condition
	 */
	RCP<OrCondition> orCon1 = rcp(new OrCondition(conList1));
	TEST_ASSERT(orCon1->isConditionTrue());
	RCP<OrCondition> orCon2 = rcp(new OrCondition(conList2));
	TEST_ASSERT(orCon2->isConditionTrue());
	RCP<OrCondition> orCon3 = rcp(new OrCondition(conList3));
	TEST_ASSERT(!orCon3->isConditionTrue());

	/*
	 * Testing equal condition
	 */
	RCP<EqualsCondition> equalsCon1 = rcp(new EqualsCondition(conList1));
	TEST_ASSERT(equalsCon1->isConditionTrue());
	RCP<EqualsCondition> equalsCon2 = rcp(new EqualsCondition(conList2));
	TEST_ASSERT(!equalsCon2->isConditionTrue());
	RCP<EqualsCondition> equalsCon3 = rcp(new EqualsCondition(conList3));
	TEST_ASSERT(equalsCon3->isConditionTrue());
}

//Test getters and setters
TEUCHOS_UNIT_TEST(Teuchos_Conditions, testConditionGetterAndSetters){
	//Settin up initial list
	RCP<ParameterList> testingList = rcp(new ParameterList("Condition Testing List"));

	Array<std::string> validValues(tuple<std::string>("mountain dew", "pepsi", "coke", "fanta"));
	RCP<StringValidator> stringVali1 = rcp(new StringValidator(validValues));

	testingList->set("string param", "fanta", "parameter for testing string conditions", stringVali1);

	StringCondition::ValueList conValues1(tuple<std::string>("pepsi", "coke"));
	RCP<StringCondition> stringCon1 = rcp( new StringCondition(testingList->getEntryRCP("string param"), conValues1));
	Dependency::ConstParameterEntryList stringParameters = stringCon1->getAllParameters();
	TEST_ASSERT(stringParameters.size() == 1);
	TEST_ASSERT(stringParameters.find(testingList->getEntryRCP("string param")) != stringParameters.end());

	/*
	 * Testing for number condition
	 */
	testingList->set("double param", 5.0, "parameter for testing number conditions");

	RCP<NumberCondition<double> > numberCon1 = rcp( new NumberCondition<double>(testingList->getEntryRCP("double param")));
	Dependency::ConstParameterEntryList numberParameters = numberCon1->getAllParameters();
	TEST_ASSERT(numberParameters.size() == 1);
	TEST_ASSERT(numberParameters.find(testingList->getEntryRCP("double param")) != numberParameters.end());

	/*
	 * Testing bool conditions
	 */
	testingList->set("bool param", true, "parameter for testing bool conditions");

	RCP<BoolCondition> boolCon1 = rcp( new BoolCondition(testingList->getEntryRCP("bool param")));
	Dependency::ConstParameterEntryList boolParameters = boolCon1->getAllParameters();
	TEST_ASSERT(boolParameters.size() == 1);
	TEST_ASSERT(boolParameters.find(testingList->getEntryRCP("bool param")) != boolParameters.end());

	/*
	 * Test Not condition
	 */
	RCP<NotCondition> notCon1 = rcp(new NotCondition(numberCon1));
	Dependency::ConstParameterEntryList notParameters = notCon1->getAllParameters();
	TEST_ASSERT(notParameters.size() == 1);
	TEST_ASSERT(notParameters.find(testingList->getEntryRCP("double param")) != notParameters.end());

	/*
	 * Test And condition
	 */
	Condition::ConstConditionList conList1(tuple<RCP<const Condition> >(stringCon1, boolCon1));
	RCP<AndCondition> andCon1 = rcp(new AndCondition(conList1));
	Dependency::ConstParameterEntryList andParameters = andCon1->getAllParameters();
	TEST_ASSERT(andParameters.size() == 2);
	TEST_ASSERT(andParameters.find(testingList->getEntryRCP("string param")) != andParameters.end());
	TEST_ASSERT(andParameters.find(testingList->getEntryRCP("bool param")) != andParameters.end());

	/*
	 * Testing or condition
	 */
	RCP<OrCondition> orCon1 = rcp(new OrCondition(conList1));
	Dependency::ConstParameterEntryList orParameters = orCon1->getAllParameters();
	TEST_ASSERT(orParameters.size() == 2);
	TEST_ASSERT(orParameters.find(testingList->getEntryRCP("string param")) != orParameters.end());
	TEST_ASSERT(orParameters.find(testingList->getEntryRCP("bool param")) != orParameters.end());

	/*
	 * Testing Equsl condition
	 */
	Condition::ConstConditionList conList2(tuple<RCP<const Condition> >(numberCon1, boolCon1));
	RCP<EqualsCondition> equalsCon1 = rcp(new EqualsCondition(conList2));
	Dependency::ConstParameterEntryList equalsParameters = equalsCon1->getAllParameters();
	TEST_ASSERT(equalsParameters.size() == 2);
	TEST_ASSERT(equalsParameters.find(testingList->getEntryRCP("double param")) != equalsParameters.end());
	TEST_ASSERT(equalsParameters.find(testingList->getEntryRCP("bool param")) != equalsParameters.end());

	/*
	 * Testing BinaryLogicCondition add
	 */
	equalsCon1->addCondition(orCon1);
	Dependency::ConstParameterEntryList equalsParameters2 = equalsCon1->getAllParameters();
	TEST_ASSERT(equalsParameters2.size() == 3);
	TEST_ASSERT(equalsParameters2.find(testingList->getEntryRCP("string param")) != equalsParameters2.end());
	TEST_ASSERT(equalsParameters2.find(testingList->getEntryRCP("double param")) != equalsParameters2.end());
	TEST_ASSERT(equalsParameters2.find(testingList->getEntryRCP("bool param")) != equalsParameters2.end());

}

//Test that exceptions get thrown when they should.
TEUCHOS_UNIT_TEST(Teuchos_Conditions, testConditionException){
	//Settin up initial list
	RCP<ParameterList> testingList = rcp(new ParameterList("Condition Testing List"));
	testingList->set("double param",1.0);
	testingList->set("string param", "awesome");
	RCP<ParameterList> testingList2 = rcp(new ParameterList("Condition Testing List"));
	testingList2->set("bool param", true);

	TEST_THROW(BoolCondition boolCon1(testingList->getEntryRCP("bool param")), InvalidConditionException);
	TEST_THROW(StringCondition stringCon1(testingList->getEntryRCP("double param"), "coke"), InvalidConditionException);
	TEST_THROW(BoolCondition boolCon1(testingList->getEntryRCP("double param")), InvalidConditionException);
	Condition::ConstConditionList conList1;
	TEST_THROW(AndCondition andCon1(conList1), InvalidConditionException);
	RCP<const Condition> con1;
	TEST_THROW(NotCondition notCon1(con1), InvalidConditionException);
}

} //namespace Teuchos
