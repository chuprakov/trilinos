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

#ifndef PHX_MD_FIELD_H
#define PHX_MD_FIELD_H

#include <iostream>
#include <string>
#include "Teuchos_ArrayRCP.hpp"
#include "Phalanx_Array.hpp"
#include "Phalanx_FieldTag_Tag.hpp"

namespace PHX {

  template<typename DataT,
	   typename Tag0 = void, typename Tag1 = void, typename Tag2 = void, 
	   typename Tag3 = void, typename Tag4 = void, typename Tag5 = void,
	   typename Tag6 = void, typename Tag7 = void>
  class MDField {
    
  public:

    typedef DataT value_type;

    typedef typename phdmesh::ArrayNatural<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> array_type;
    
    typedef typename array_type::size_type size_type;

    typedef typename array_type::index_type index_type;

    MDField(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& t);
    
    MDField(const PHX::Tag<DataT>& v);
    
    MDField();
    
    ~MDField();
    
    const PHX::FieldTag& fieldTag() const;

    DataT& operator()(index_type index1, index_type index2, index_type index3, 
		      index_type index4, index_type index5, index_type index6,
		      index_type index7, index_type index8);

    DataT& operator()(index_type index1, index_type index2, index_type index3, 
		      index_type index4, index_type index5, index_type index6,
		      index_type index7);

    DataT& operator()(index_type index1, index_type index2, index_type index3, 
		      index_type index4, index_type index5, index_type index6);
    
    DataT& operator()(index_type index1, index_type index2, index_type index3, 
		      index_type index4, index_type index5);
    
    DataT& operator()(index_type index1, index_type index2, index_type index3, 
		      index_type index4);
    
    DataT& operator()(index_type index1, index_type index2, index_type index3);
    
    DataT& operator()(index_type index1, index_type index2);
    
    DataT& operator()(index_type index1);
    
    DataT& operator[](index_type index);

    //typename phdmesh::ArrayNatural<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::size_type size() const;

    size_type size() const;

    void setFieldTag(const PHX::Tag<DataT>& t);
    
    void setFieldData(const Teuchos::ArrayRCP<DataT>& d);
    
    void print(std::ostream& os, bool printValues = false) const;

  private:
    
    PHX::Tag<DataT> m_tag;
    
    array_type m_field_data;

    Teuchos::ArrayRCP<DataT> m_array_rcp;

#ifdef PHX_DEBUG
    bool m_tag_set;
    bool m_data_set;
    static const std::string m_field_tag_error_msg;
    static const std::string m_field_data_error_msg;
#endif

  };
  
  template<typename DataT,
	   typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	   typename Tag4, typename Tag5, typename Tag6, typename Tag7>
  std::ostream& operator<<(std::ostream& os, 
			   const PHX::MDField<DataT, Tag0, Tag1, Tag2, Tag3,
			   Tag4, Tag5, Tag6, Tag7>& h);
  
} 

#include "Phalanx_MDField_Def.hpp"

#endif 
