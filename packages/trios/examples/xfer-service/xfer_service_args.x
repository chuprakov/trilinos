/* ************************************************************************

                   Trios: Trilinos I/O Support
                 Copyright 2011 Sandia Corporation

 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 the U.S. Government retains certain rights in this software.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 1. Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.

 3. Neither the name of the Corporation nor the names of the
 contributors may be used to endorse or promote products derived from
 this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)

*************************************************************************/
/* -------------------------------------------------------------------------- */
/**
 *   @file xfer_xdr.x
 *
 *   @brief Type definitions for a simple rpc service.
 *
 *   @author Ron Oldfield (raoldfi\@cs.sandia.gov).
 *   $Revision: 342 $.
 *   $Date: 2005-05-01 23:30:57 -0600 (Sun, 01 May 2005) $.
 *
 */

/**
 * @defgroup xfer_example  Nessie Data Transfer Example
 *
 * The data-transfer example demonstrates a simple client
 * and server that transfer an array of 16-byte \ref data_t
 * data structures from a parallel application to a set of
 * servers.  We implemented two variations:
 * one that transfers the array of data structures
 * with the request, and a second method that has each server
 * pull the data using the \ref nssi_get_data() function.  Although
 * this example is fairly simple, it makes a decent benchmark code
 * to evaluate overheads of the Nessie transfer protocols and encoding
 * schemes.
 *
*/

/**
 * @defgroup xfer_types  Nessie Example Types
 * @ingroup xfer_example
 *
 * @{
 */

/* Extra stuff to put at the beginning of the header file */
#ifdef RPC_HDR
#endif

/* Extra stuff to put at the beginning of the C file. */
#ifdef RPC_XDR
%#include <xfer_service_args.h>
#endif



/**
 * @brief Opcodes for the types of transfer operations.  
 */
enum xfer_op {
    /** Opcode for the operation to transfer the array with the request. */
    XFER_PUSH = 1,
    /** Opcode for the operation that pulls the data array to the server. */
    XFER_PULL,
    /**  */
    XFER_ROUNDTRIP,
    /**  */
    XFER_GET,
    /**  */
    XFER_PUT
};

/**
 * @brief A 16-byte structure that contains an int, float, and double.
 *
 * This structure contains an int, float, and double as an example
 * of a complex structure with multiple types.  This will exercise the
 * encoding/decoding features of Nessie.
 */
struct data_t {
    /** An integer value. */
    uint32_t int_val;
    /** A floating point value. */
    float float_val;
    /** A double value. */
    double double_val;
};

/**
 * @brief Array of 16-byte structures that we can send with the request.
 *
 * Rpcgen will use this definition to define encoding functions to
 * encode and decode an array of \ref data_t structures.  We will
 * use these functions when sending the array with the request.
 */
typedef data_t data_array_t<>;

/**
 * @brief Arguments for the first transfer operation (PUSH).
 *
 * The first transfer operation includes the array of \ref data_t
 * structures as an argument of the remote operation.  This will
 * cause the array to be sent to the server as part of the request.
 */
struct xfer_push_args {
    /** The array of \ref data_t structures, including length. */
    data_array_t array;
};

/**
 * @brief Arguments for the second transfer operation (PULL).
 *
 * The second transfer operation only needs to send the length
 * of the data array.  It uses the data argument in the nssi_call_rpc()
 * to identify the raw data for the server to fetch.
 */
struct xfer_pull_args {
    /** The length of the data array. */
    int32_t len;
};

/**
 * @brief Arguments for the third transfer operation (ROUNDTRIP).
 *
 * The third transfer operation includes the array of \ref data_t
 * structures as an argument of the remote operation.  This will
 * cause the array to be sent to the server as part of the request.
 * If the size of array is large, the address of array will be sent
 * as part of the request and the server will fetch array via RDMA.
 */
struct xfer_roundtrip_args {
    /** The array of \ref data_t structures, including length. */
    data_array_t array;
};
/**
 * @brief Results for the third transfer operation (ROUNDTRIP).
 *
 * The third transfer operation includes the array of \ref data_t
 * structures as a result of the remote operation.  This will
 * cause the array to be sent to the client as part of the result.
 * If the size of array is large, the address of array will be sent
 * as part of the result and the client will fetch array via RDMA.
 */
struct xfer_roundtrip_res {
    /** The array of \ref data_t structures, including length. */
    data_array_t array;
};

/**
 * @brief Arguments for the fourth transfer operation (GET).
 *
 * The fourth transfer operation only needs to send the length
 * of the data array.  It uses the data argument in the nssi_call_rpc()
 * to identify the raw data for the server to fetch.
 */
struct xfer_get_args {
    /** The length of the data array. */
    int32_t len;
};

/**
 * @brief Results for the fourth transfer operation (GET).
 *
 * The fourth transfer operation includes the array of \ref data_t
 * structures as a result of the remote operation.  This will
 * cause the array to be sent to the client as part of the result.
 * If the size of array is large, the address of array will be sent
 * as part of the result and the client will fetch array via RDMA.
 */
struct xfer_get_res {
    /** The array of \ref data_t structures, including length. */
    data_array_t array;
};

/**
 * @brief Arguments for the fifth transfer operation (PUT).
 *
 * The fifth transfer operation includes the array of \ref data_t
 * structures as an argument of the remote operation.  This will
 * cause the array to be sent to the server as part of the request.
 */
struct xfer_put_args {
    /** The array of \ref data_t structures, including length. */
    data_array_t array;
};

/**
 * @brief Results for the fifth transfer operation (PUT).
 *
 * The fifth transfer operation only needs to send the length
 * of the data array as the result.  The client uses the data
 * argument in the nssi_call_rpc() to identify the address of
 * the raw data buffer for the server to put into.
 */
struct xfer_put_res {
    /** The length of the data array. */
    int32_t len;
};

/**
 * @}
 */
