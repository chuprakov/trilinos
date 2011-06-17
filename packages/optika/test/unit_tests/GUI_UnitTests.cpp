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
private:
  static inline QModelIndex getWidgetIndex(const QModelIndex& index);
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
}*/

#define VERIFY_PARAMETER_TYPE(PL, NAME, TYPE, MODEL) \
  GET_ENTRY_INDEX( PL , NAME , MODEL ) \
  QCOMPARE( MODEL->data( NAME##Index, Qt::DisplayRole).toString(), \
    QString::fromStdString( #NAME) );  \
  QModelIndex NAME##TypeIndex = NAME##Index.sibling(NAME##Index.row(),2); \
  QVERIFY( NAME##TypeIndex.isValid()); \
  QCOMPARE( MODEL->data( NAME##TypeIndex, Qt::DisplayRole).toString(), TYPE );


void OptikaGUITests::typeTest(){
  RCP<ParameterList> My_List = 
    RCP<ParameterList>(new ParameterList);

  double *pointer = 0;
  My_List->set("Doublepointer", pointer);
  My_List->set(
    "MaxIters", 
    1550, 
    "Determines the maximum number of iterations in the solver");
  My_List->set(
    "Tolerance", 1e-10, "The tolerance used for the convergence check");
  
  TreeModel* model = new TreeModel(My_List);

  VERIFY_PARAMETER_TYPE(My_List, MaxIters, intId, model)
  VERIFY_PARAMETER_TYPE(My_List, Doublepointer, unrecognizedId, model);
  VERIFY_PARAMETER_TYPE(My_List, Tolerance, doubleId, model);
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
  RCP<DependencySheet> dependencySheet = rcp(new DependencySheet);
  RCP<ParameterList> validParameters = 
    getParametersFromXmlFile("deptests.xml", dependencySheet);
  TreeModel* model = new TreeModel(validParameters, dependencySheet);
  Delegate* delegate = new Delegate;
  TreeView* treeView = new TreeView(model, delegate);
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

  delete model;
  delete treeView;
  delete delegate;
}
  


} //namespace Optika

QTEST_MAIN(Optika::OptikaGUITests)
#include "GUI_UnitTests.moc"

