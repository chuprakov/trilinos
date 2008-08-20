// @HEADER
// ************************************************************************
// 
//            Phalanx: A Partial Differential Equation Assembly 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#ifndef PHX_MD_FIELD_DEF_H
#define PHX_MD_FIELD_DEF_H

#include <sstream>
#include "Teuchos_TestForException.hpp"
#include "Phalanx_Print_Utilities.hpp"

//**********************************************************************
#ifdef PHX_DEBUG
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
const std::string PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::m_field_tag_error_msg = 
  "Error - PHX::MDField::fieldTag() - No tag has been set!";

template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
const std::string PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::m_field_data_error_msg = 
  "Error - PHX::MDField::operator[] - No data has been set!  Please call getFieldData(this) on all PHX::MDField objects in providers!";
#endif

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDField(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& t) :
  m_tag(name,t)
#ifdef PHX_DEBUG
  , m_tag_set(true),
  m_data_set(false)
#endif
{ }

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDField(const PHX::Tag<DataT>& v) :
  m_tag(v)
#ifdef PHX_DEBUG
  ,m_tag_set(true),
  m_data_set(false)
#endif
{ }

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::MDField() :
  m_tag("???", Teuchos::null)
#ifdef PHX_DEBUG
  ,m_tag_set(false),
  m_data_set(false)
#endif
{ }

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::~MDField()
{ }

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
const PHX::FieldTag& 
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::fieldTag() const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_tag_error_msg);
#endif
  return m_tag;
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(index_type index1, index_type index2, index_type index3, 
	   index_type index4, index_type index5, index_type index6,
	   index_type index7, index_type index8)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4,index5,index6,index7,index8);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(index_type index1, index_type index2, index_type index3, 
	   index_type index4, index_type index5, index_type index6,
	   index_type index7)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4,index5,index6,index7);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(index_type index1, index_type index2, index_type index3, 
	   index_type index4, index_type index5, index_type index6)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4,index5,index6);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(index_type index1, index_type index2, index_type index3, 
	   index_type index4, index_type index5)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4,index5);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(index_type index1, index_type index2, index_type index3, 
	   index_type index4)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(index_type index1, index_type index2, index_type index3)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(index_type index1, index_type index2)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(index_type index1)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator[](index_type index)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data[index];
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::size_type 
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::size() const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data.size();
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
void PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
setFieldTag(const PHX::Tag<DataT>& v)
{  
#ifdef PHX_DEBUG
  m_tag_set = true;
#endif
  m_tag = v;
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
void PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
setFieldData(const Teuchos::ArrayRCP<DataT>& d)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_tag_error_msg);
  m_data_set = true;
#endif

  m_array_rcp = d;

  std::vector<size_type> data_layout_dim;
  m_tag.dataLayout().dimensions(data_layout_dim);

  /* Commenting out until Carter adds new power user ctor to array.
  std::vector<size_type> array_dim;
  size_type num_cells = d.size() / m_tag.dataLayout().size();

  array_dim.push_back(num_cells);
  for (std::size_t i = 0; i < data_layout_dim.size(); ++i)
    array_dim.push_back(data_layout_dim[i]); 

  TEST_FOR_EXCEPTION(array_dim.size() != (data_layout_dim.size() + 1), 
		     std::logic_error, 
		     "This should never happen!  STL seems broken.");

  if (array_type::Rank != (data_layout_dim.size() + 1) ) {
    std::ostringstream os;
    os << "Array rank must be equal to the DataLayout rank + 1.\n"
       << "Array rank = " << array_type::Rank 
       << ", data layout rank = " << data_layout_dim.size() << " + 1 = " 
       << (data_layout_dim.size() + 1) << std::endl
       << "Offending MDField:\n"
       << *this << std::endl;
    TEST_FOR_EXCEPTION(array_type::Rank != (data_layout_dim.size() + 1), 
		       std::logic_error, os.str());
  }

  typename phdmesh::ArrayNatural<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> array(d.get(), phdmesh::ArrayStride(), &array_dim[0]);
  
  */

  // stride is multiplying dimesnsions in reverse order
  // Example: Array<cell,node,dim>(100,4,3) has a stride of
  // stride[0] = 3;
  // stride[1] = 12;
  // stride[2] = 120;
  size_type num_cells = d.size() / m_tag.dataLayout().size();
  std::vector<size_type> stride(data_layout_dim.size() + 1);

  int stride_size = 1;
  for (std::size_t i = 0; i < data_layout_dim.size(); ++i) {
    stride_size *= data_layout_dim[i];
    stride[i] = stride_size;
  }
  stride[stride.size()-1] = stride_size * num_cells;

//   std::cout << std::endl;
//   for (std::size_t i = 0; i < stride.size(); ++i)
//     std::cout << "stride[" << i << "] = " << stride[i] << std::endl;

  TEST_FOR_EXCEPTION(stride.size() != (data_layout_dim.size() + 1), 
		     std::logic_error, 
		     "This should never happen!  STL seems broken.");

  if (array_type::Rank != (data_layout_dim.size() + 1) ) {
    std::ostringstream os;
    os << "Array rank must be equal to the DataLayout rank + 1.\n"
       << "Array rank = " << array_type::Rank 
       << ", data layout rank = " << data_layout_dim.size() << " + 1 = " 
       << (data_layout_dim.size() + 1) << std::endl
       << "Offending MDField:\n"
       << *this << std::endl;
    TEST_FOR_EXCEPTION(array_type::Rank != (data_layout_dim.size() + 1), 
		       std::logic_error, os.str());
  }
  
  typename phdmesh::ArrayNatural<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> array(d.get(), phdmesh::ArrayStride(), &stride[0]);
  
  m_field_data = array;

}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
void PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
print(std::ostream& os,	bool printValues) const
{

  std::vector<const char*> dim_names;

  PHX::PrintDimension<Tag0,array_type> pd0;
  pd0.addName(dim_names);
  PHX::PrintDimension<Tag1,array_type> pd1;
  pd1.addName(dim_names);
  PHX::PrintDimension<Tag2,array_type> pd2;
  pd2.addName(dim_names);
  PHX::PrintDimension<Tag3,array_type> pd3;
  pd3.addName(dim_names);
  PHX::PrintDimension<Tag4,array_type> pd4;
  pd4.addName(dim_names);
  PHX::PrintDimension<Tag5,array_type> pd5;
  pd5.addName(dim_names);
  PHX::PrintDimension<Tag6,array_type> pd6;
  pd6.addName(dim_names);
  PHX::PrintDimension<Tag7,array_type> pd7;
  pd7.addName(dim_names);

  os << "MDField<";

  for (std::size_t i=0; i < dim_names.size(); ++i) {
    if (i > 0)
      os << ",";
    os << std::string(dim_names[i]);
  }
  os << ">(";
  for (std::size_t i=0; i < dim_names.size(); ++i) {
    if (i > 0)
      os << ",";
    os << m_field_data.dimension(i);
  }
  os << "): ";
  
  os << m_tag;

  if (printValues) {
    os << std::endl;
    for (std::size_t i = 0; i < m_field_data.size(); ++i)
      os << "value[" << i << "] = " << m_field_data[i] << std::endl;
  }

}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
std::ostream& PHX::operator<<(std::ostream& os, 
			      const PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,
			      Tag4,Tag5,Tag6,Tag7>& f)
{
  f.print(os, false);
  return os;
}

//**********************************************************************

#endif
