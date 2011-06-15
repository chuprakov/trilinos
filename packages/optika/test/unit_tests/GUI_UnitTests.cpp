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

#include <QtGui>
#include <QtTest/QtTest>
#include "Teuchos_UnitTestHarness.hpp"
#include "Optika_treemodel.hpp"
#include "Optika_delegate.hpp"
#include "Optika_treeview.hpp"
#include <QApplication>
#include <QSpinBox>
#include <iostream>

/**
 *  THIS IS JUST A PROOF OF CONCEPT!
 *  I just wanted to see if I could get the QTest
 *  framework working with our CTest Framework.
 *  I don't have time right now to actually
 *  do extensive GUI unit testing. I will come
 *  back and actually write the tests a.s.a.p.
 *
 *  kln 21/01/2011
 */
namespace Optika{

class OptikaGUITests: public QObject{
Q_OBJECT
private slots:
  //void basicTests();
  void dependencyTests();
private:
  /*static void verifyParameterAndType(
    RCP<ParameterList> pl,
    QString name, 
    QString type, 
    TreeModel* model);

  static QModelIndex getEntryIndex(
    RCP<ParameterList> pl, 
    std::string name,
    TreeModel* model);*/
};

//QModelIndex OptikaGUITests::getEntryIndex(

#define GET_ENTRY_INDEX(\
  PL, \
  NAME, \
  MODEL) \
  RCP<ParameterEntry> NAME##Entry = PL->getEntryRCP( #NAME  ); \
  QVERIFY(nonnull( NAME##Entry )); \
  QModelIndex NAME##Index = MODEL->findParameterEntryIndex( NAME##Entry ); \
  QVERIFY( NAME##Index.isValid());

/*
void OptikaGUITests::verifyParameterAndType(
  RCP<ParameterList> pl,
  QString name, 
  QString type, 
  TreeModel* model)
{
  QModelIndex itemIndex = getEntryIndex(pl, name.toStdString(), model);
  QCOMPARE(model->data(itemIndex, Qt::DisplayRole).toString(),name);
  QModelIndex typeIndex = itemIndex.sibling(itemIndex.row(),2);
  QVERIFY(typeIndex.isValid());
  QCOMPARE(model->data(typeIndex, Qt::DisplayRole).toString(),type);

}

void OptikaGUITests::basicTests(){
  RCP<ParameterList> My_List = 
    RCP<ParameterList>(new ParameterList);

  double *pointer = 0;
  My_List->set("Double pointer", pointer);
  My_List->set(
    "Max Iters", 
    1550, 
    "Determines the maximum number of iterations in the solver");
  My_List->set(
    "Tolerance", 1e-10, "The tolerance used for the convergence check");
  
  TreeModel* model = new TreeModel(My_List);

  verifyParameterAndType(My_List, "Max Iters", intId, model);
  verifyParameterAndType(My_List, "Double pointer", unrecognizedId, model);
  verifyParameterAndType(My_List, "Tolerance", doubleId, model);
}*/


void OptikaGUITests::dependencyTests(){
  RCP<DependencySheet> dependencySheet = rcp(new DependencySheet);
  RCP<ParameterList> validParameters = 
    getParametersFromXmlFile("deptests.xml", dependencySheet);
  TreeModel* model = new TreeModel(validParameters, dependencySheet);
  Delegate* delegate = new Delegate;
  TreeView* treeView = new TreeView(model, delegate);
  GET_ENTRY_INDEX(validParameters, Preconditioner, model)
  QVERIFY(treeView->isRowHidden(
    PreconditionerIndex.row(), PreconditionerIndex.parent()));

}
  


} //namespace Optika

/**
 *  Test the basics of Optika Functionality
 */
/*TEUCHOS_UNIT_TEST(Optika_GUITests, testBasics){
  int argc = 0;
  char* argv[1];
  char name [] = "optika_basic_test";
  argv[0]=name;
  QApplication app(argc, argv);
  QTEST_DISABLE_KEYPAD_NAVIGATION
  Optika::BasicTests bt;
  TEUCHOS_ASSERT(QTest::qExec(&bt, argc, argv) == 0);
}*/

QTEST_MAIN(Optika::OptikaGUITests)
#include "GUI_UnitTests.moc"

