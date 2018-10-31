/****************************************************************************

	FloatEditDelegate Class
	Copyright (C) 2009 Andre Deperrois adeperrois@xflr5.com

	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*****************************************************************************/


#include "FloatEditDelegate.h"


FloatEditDelegate::FloatEditDelegate(QObject *parent)
: QItemDelegate(parent)
{
    QFEM = false;
}

QWidget *FloatEditDelegate::createEditor(QWidget *parent, const QStyleOptionViewItem &/*option*/, const QModelIndex & index ) const
{
	if(m_Precision[index.column()]>=0)
	{
		//we have a number
		NumberEdit *editor = new NumberEdit();
		editor->setParent(parent);
		editor->setAutomaticPrecision(m_Precision[index.column()]);
		double value = index.model()->data(index, Qt::EditRole).toDouble();
		editor->setValue(value);
        if (QFEM)
        {
            if(index.column()==0) editor->setEnabled(false);
        }
		return editor;
	}
	else
	{
		//we have a string
		QLineEdit *editor = new QLineEdit(parent);
		editor->setAlignment(Qt::AlignLeft);
		return editor;
	}


}

void FloatEditDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const
{
	if(m_Precision[index.column()]>=0)
	{
		double value = index.model()->data(index, Qt::EditRole).toDouble();
		NumberEdit *floatEdit = static_cast<NumberEdit*>(editor);
		floatEdit->setValue(value);
	}
	else
	{
		QLineEdit *pLine = static_cast<QLineEdit*>(editor);
		pLine->setText(index.model()->data(index, Qt::EditRole).toString());
	}
}



void FloatEditDelegate::setModelData(QWidget *editor, QAbstractItemModel *model,
								const QModelIndex &index) const
{
	if(m_Precision[index.column()]>=0)
	{
		NumberEdit *floatEdit = static_cast<NumberEdit*>(editor);
		double value = floatEdit->getValue(true);
		model->setData(index, value, Qt::EditRole);
	}
	else
	{
		QLineEdit *pLine = static_cast<QLineEdit*>(editor);
		model->setData(index, pLine->text(), Qt::EditRole);
	}
}


void FloatEditDelegate::updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &/* index */) const
{
	editor->setGeometry(option.rect);
}



void FloatEditDelegate::SetPrecision(int *PrecisionTable)
{
	m_Precision = PrecisionTable;
}



void FloatEditDelegate::paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
	QString strong;
	QStyleOptionViewItem myOption = option;
	if(m_Precision[index.column()]>=0)
	{
		myOption.displayAlignment = Qt::AlignRight | Qt::AlignVCenter;
		strong = QString("%1").arg(index.model()->data(index, Qt::DisplayRole).toDouble(),0,'f', m_Precision[index.column()]);

	}
	else
	{
		myOption.displayAlignment = Qt::AlignLeft | Qt::AlignVCenter;
		strong = index.model()->data(index, Qt::DisplayRole).toString();
	}
	drawDisplay(painter, myOption, myOption.rect, strong);
	drawFocus(painter, myOption, myOption.rect);
}















