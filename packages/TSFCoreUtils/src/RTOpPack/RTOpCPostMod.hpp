// /////////////////////////////////////////////////////////////////////////
// RTOpCPostMod.hpp

#ifndef RTOP_C_POST_MOD_HPP
#define RTOP_C_POST_MOD_HPP

#include "RTOpCppC.hpp"

namespace RTOpPack {

class RTOpCPostMod {
public:

	///
	RTOpCPostMod( const RTOp_RTOp_vtbl_t *vtbl ) : vtbl_(vtbl)
		{
#ifdef _DEBUG
			THROW_EXCEPTION( !(vtbl && vtbl->obj_data_vtbl && vtbl->obj_data_vtbl->obj_create)
							 , std::logic_error
							 , "Error!"	);
#endif			
		}
	///
	void initialize(RTOpC *op) const
		{
			op->op().vtbl = vtbl_;
			op->op().vtbl->obj_data_vtbl->obj_create(NULL,NULL,&op->op().obj_data);
		}
	
private:
	
	const RTOp_RTOp_vtbl_t *vtbl_;
	
	RTOpCPostMod(); // Not defined and not to be called.

};

} // namespace RTOpPack

#endif // RTOP_C_POST_MOD_HPP
