#include "Optika_FilenameWidget.hpp"
#include <QPushButton>
#include <QLabel>
#include <QVBoxLayout>
#include <QFileDialog>


namespace Optika{


FilenameWidget::FilenameWidget(QString currentFileName, bool mustAlreadyExist, QWidget *parent):QWidget(parent){
	this->mustAlreadyExist = mustAlreadyExist;
	this->currentFileName = currentFileName;
	QPushButton *changeButton = new QPushButton("Change Path",this);
	connect(changeButton, SIGNAL(clicked(bool)), this, SLOT(getNewFileName()));
	pathLabel = new QLabel(currentFileName,this);
	QVBoxLayout *layout = new QVBoxLayout(this);
	layout->addWidget(changeButton);
	layout->addWidget(pathLabel);
	setLayout(layout);
}

QString FilenameWidget::getCurrentFileName(){
	return currentFileName;
}

void FilenameWidget::setCurrentFileName(QString newName){
	currentFileName = newName;
	pathLabel->setText(newName);	
}

void FilenameWidget::getNewFileName(){
	QString defaultPath;
	if(currentFileName == ""){
		defaultPath = QDir::homePath();
	}
	else{
		defaultPath = currentFileName;
	}
	if(mustAlreadyExist){
		setCurrentFileName(QFileDialog::getOpenFileName(this, tr("File"), defaultPath));
	}
	else{
		setCurrentFileName(QFileDialog::getSaveFileName(this, tr("File"), defaultPath));
	}
}


}

