// /////////////////////////////////////////////////////////////////////////////
// RTOp_obj_free_free.c
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#include <malloc.h>

#include "RTOp_obj_free_free.h"

int RTOp_obj_free_free( const struct RTOp_obj_type_vtbl_t* vtbl, const void* instance_data, void** obj )
{
  free( *obj );
  *obj = NULL;
  return 0;
}
