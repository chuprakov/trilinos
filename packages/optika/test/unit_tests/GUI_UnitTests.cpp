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
  void typeTest();
  void dependencyTests();
  void cleanupTestCase();
private:
  static inline QModelIndex getWidgetIndex(const QModelIndex& index);
  QObjectCleanupHandler cleaner;
};

void OptikaGUITests::cleanupTestCase(){
  cleaner.clear();
}
  


//QModelIndex OptikaGUITests::getEntryIndex(

#define GET_ENTRY_INDEX(\
  PL, \
  NAME, \
  MODEL) \
  RCP<ParameterEntry> NAME##Entry = PL->getEntryRCP( #NAME  ); \
  QVERIFY(nonnull( NAME##Entry )); \
  QModelIndex NAME##Index = MODEL->findParameterEntryIndex( NAME##Entry ); \
  QVERIFY( NAME##Index.isValid());

#define VERIFY_PARAMETER_TYPE(PL, NAME, TYPE, MODEL) \
  GET_ENTRY_INDEX( PL , NAME , MODEL ) \
  QCOMPARE( MODEL->data( NAME##Index, Qt::DisplayRole).toString(), \
    QString::fromStdString( #NAME) );  \
  QModelIndex NAME##TypeIndex = NAME##Index.sibling(NAME##Index.row(),2); \
  QVERIFY( NAME##TypeIndex.isValid()); \
  std::cout <<  \
    MODEL->data( NAME##TypeIndex, Qt::DisplayRole).toString().toStdString() \
    <<std::endl; \
  QCOMPARE( MODEL->data( NAME##TypeIndex, Qt::DisplayRole).toString(), TYPE );


void OptikaGUITests::typeTest(){
  cleaner.clear();
  RCP<ParameterList> My_List = 
    RCP<ParameterList>(new ParameterList);

  double *pointer = 0;
  Array<double*> doubleStarArray;
  Array<int> intArray;
  My_List->set("Doublepointer", pointer);
  My_List->set(
    "MaxIters", 
    1550, 
    "Determines the maximum number of iterations in the solver");
  My_List->set(
    "Tolerance", 1e-10, "The tolerance used for the convergence check");
  My_List->set("DoublePointerArray", doubleStarArray);
  My_List->set("IntArray", intArray);
  
  TreeModel* model = new TreeModel(My_List);
  cleaner.add(model);

  VERIFY_PARAMETER_TYPE(My_List, MaxIters, intId, model)
  VERIFY_PARAMETER_TYPE(My_List, Doublepointer, unrecognizedId, model);
  VERIFY_PARAMETER_TYPE(My_List, Tolerance, doubleId, model);
  VERIFY_PARAMETER_TYPE(My_List, DoublePointerArray, unrecognizedId, model);
  VERIFY_PARAMETER_TYPE(My_List, IntArray, arrayId + " " + intId, model);
  cleaner.remove(model);
  delete model; 
}

inline QModelIndex OptikaGUITests::getWidgetIndex(const QModelIndex& index){
  return index.sibling(index.row(),1);
}
  
  
#define VERIFY_HIDDEN_ROW(INDEX) \
  QVERIFY(treeView->isRowHidden(  INDEX.row(), INDEX.parent()));

#define VERIFY_SHOWN_ROW(INDEX) \
  QVERIFY(!treeView->isRowHidden(  INDEX.row(), INDEX.parent()));

void OptikaGUITests::dependencyTests(){
  cleaner.clear();
  RCP<DependencySheet> dependencySheet = rcp(new DependencySheet);
  RCP<ParameterList> validParameters = 
    getParametersFromXmlFile("deptests.xml", dependencySheet);
  TreeModel* model = new TreeModel(validParameters, dependencySheet);
  Delegate* delegate = new Delegate;
  TreeView* treeView = new TreeView(model, delegate);
  cleaner.add(model);
  cleaner.add(delegate);
  cleaner.add(treeView);
  QStyleOptionViewItem genericStyleItem;
  
//Testing Bool visual dependency
  GET_ENTRY_INDEX(validParameters, Preconditioner, model)
  GET_ENTRY_INDEX(validParameters, ShowPrecs, model)

  VERIFY_HIDDEN_ROW(PreconditionerIndex)
  QModelIndex precWidgetIndex = getWidgetIndex(ShowPrecsIndex);
  QComboBox* precBox = (QComboBox*)delegate->createEditor(
    0, genericStyleItem, precWidgetIndex);
  precBox->setCurrentIndex(precBox->findText(Delegate::getBoolEditorTrue()));
  delegate->setModelData(precBox, model, precWidgetIndex);
  VERIFY_SHOWN_ROW(PreconditionerIndex)


//StringVisualDependency testing
  GET_ENTRY_INDEX(validParameters, Favorite_Cheese, model)
  GET_ENTRY_INDEX(validParameters, Swiss_rating, model)

  VERIFY_HIDDEN_ROW(Swiss_ratingIndex)
  QModelIndex cheeseWidgetIndex = getWidgetIndex(Favorite_CheeseIndex);
  QComboBox* cheeseBox = (QComboBox*)delegate->createEditor(
    0,genericStyleItem, cheeseWidgetIndex);
  cheeseBox->setCurrentIndex(cheeseBox->findText("Swiss"));
  delegate->setModelData(cheeseBox, model, cheeseWidgetIndex);
  VERIFY_SHOWN_ROW(Swiss_ratingIndex)
 
//Testing Number Visual Dependencies
  GET_ENTRY_INDEX(validParameters, Temp, model)
  GET_ENTRY_INDEX(validParameters, Num_ice_cubes, model)
  VERIFY_SHOWN_ROW(Num_ice_cubesIndex)
  QModelIndex tempWidgetIndex = getWidgetIndex(TempIndex);
  QDoubleSpinBox* tempSpinner = (QDoubleSpinBox*)delegate->createEditor(
    0,genericStyleItem, tempWidgetIndex);
  tempSpinner->setValue(33.0);
  delegate->setModelData(tempSpinner, model, tempWidgetIndex);
  VERIFY_HIDDEN_ROW(Num_ice_cubesIndex)

  //Test condition visual dependency
  GET_ENTRY_INDEX(validParameters, ParamA, model)
  GET_ENTRY_INDEX(validParameters, ParamB, model)
  GET_ENTRY_INDEX(validParameters, OptParam, model)
  VERIFY_SHOWN_ROW(OptParamIndex)
  QModelIndex paramAWidgetIndex = getWidgetIndex(ParamAIndex);
  QModelIndex paramBWidgetIndex = getWidgetIndex(ParamBIndex);
  QSpinBox* paramASpinner = (QSpinBox*)delegate->createEditor(
    0, genericStyleItem, paramAWidgetIndex);
  QSpinBox* paramBSpinner = (QSpinBox*)delegate->createEditor(
    0, genericStyleItem, paramBWidgetIndex);
  paramASpinner->setValue(0);
  delegate->setModelData(paramASpinner, model, paramAWidgetIndex);
  VERIFY_SHOWN_ROW(OptParamIndex)
  paramBSpinner->setValue(0);
  delegate->setModelData(paramBSpinner, model, paramBWidgetIndex);
  VERIFY_HIDDEN_ROW(OptParamIndex)
  paramBSpinner->setValue(1);
  delegate->setModelData(paramBSpinner, model, paramBWidgetIndex);
  VERIFY_SHOWN_ROW(OptParamIndex)

  




  //Test Number Array Length Dependency
  GET_ENTRY_INDEX(validParameters, NumBuckets, model)
  GET_ENTRY_INDEX(validParameters, AmtInBuckets, model)
  Array<double> bucketsArray = model->getArray<double>(AmtInBucketsIndex);
  QCOMPARE(bucketsArray.size(),(Array<double>::size_type)3);
  QModelIndex numBucketsWidgetIndex = getWidgetIndex(NumBucketsIndex);
  QSpinBox* numBucketsSpinner = (QSpinBox*)delegate->createEditor(
    0, genericStyleItem, numBucketsWidgetIndex);
  numBucketsSpinner->setValue(5);
  delegate->setModelData(numBucketsSpinner, model, numBucketsWidgetIndex);
  bucketsArray = model->getArray<double>(AmtInBucketsIndex);
  QCOMPARE(bucketsArray.size(),(Array<double>::size_type)5);
  numBucketsSpinner->setValue(0);
  delegate->setModelData(numBucketsSpinner, model, numBucketsWidgetIndex);
  bucketsArray = model->getArray<double>(AmtInBucketsIndex);
  VERIFY_HIDDEN_ROW(AmtInBucketsIndex)
  QCOMPARE(bucketsArray.size(),(Array<double>::size_type)0);
  numBucketsSpinner->setValue(2);
  delegate->setModelData(numBucketsSpinner, model, numBucketsWidgetIndex);
  bucketsArray = model->getArray<double>(AmtInBucketsIndex);
  VERIFY_SHOWN_ROW(AmtInBucketsIndex)
  QCOMPARE(bucketsArray.size(),(Array<double>::size_type)2);







  cleaner.remove(model);
  cleaner.remove(treeView);
  cleaner.remove(delegate);
  delete model;
  delete treeView;
  delete delegate;
}
  


} //namespace Optika

QTEST_MAIN(Optika::OptikaGUITests)
#include "GUI_UnitTests.moc"

