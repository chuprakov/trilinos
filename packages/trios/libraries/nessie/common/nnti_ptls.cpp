/**
 * nnti_ptls.c
 *
 *  Created on: Jan 13, 2011
 *      Author: thkorde
 */

#include "Trios_config.h"
#include "Trios_threads.h"
#include "Trios_signal.h"

// MPI is only used to increment the PID
#ifdef HAVE_TRIOS_MPI
#include <mpi.h>
#endif

#include <assert.h>
#include <string.h>

#include "nnti_ptls.h"
#include "nnti_utils.h"



#define PTL_OP_PUT_INITIATOR  1
#define PTL_OP_GET_INITIATOR  2
#define PTL_OP_PUT_TARGET     3
#define PTL_OP_GET_TARGET     4
#define PTL_OP_SEND           5
#define PTL_OP_NEW_REQUEST    6
#define PTL_OP_RESULT         7
#define PTL_OP_RECEIVE        8

/**
 * @brief Track the state of a PtlPut (i'm the initiator).
 */
typedef struct {
    int send_start;
    int send_end;
    int ack;
} ptl_put_initiator_state_t;
/**
 * @brief Track the state of a PtlGet (i'm the initiator).
 */
typedef struct {
    int send_start;
    int send_end;
    int reply_start;
    int reply_end;
} ptl_get_initiator_state_t;
/**
 * @brief Track the state of a PtlPut (i'm the target).
 */
typedef struct {
    int put_start;
    int put_end;
} ptl_put_target_state_t;
/**
 * @brief Track the state of a PtlGet (i'm the target).
 */
typedef struct {
    int get_start;
    int get_end;
} ptl_get_target_state_t;

typedef union {
    ptl_put_initiator_state_t put_initiator;
    ptl_get_initiator_state_t get_initiator;
    ptl_put_target_state_t    put_target;
    ptl_get_target_state_t    get_target;
} ptl_op_state_t;

typedef enum {
    REQUEST_BUFFER,
    RESULT_BUFFER,
    RECEIVE_BUFFER,
    SEND_BUFFER,
    GET_SRC_BUFFER,
    GET_DST_BUFFER,
    PUT_SRC_BUFFER,
    PUT_DST_BUFFER,
    RDMA_TARGET_BUFFER,
    UNKNOWN_BUFFER
} ptl_buffer_type;

typedef struct portals_memory_handle {
    ptl_buffer_type  type;

    ptl_handle_eq_t  eq_h;
    ptl_pt_index_t   buffer_id;
    ptl_process_id_t match_id;
    ptl_match_bits_t match_bits;
    ptl_match_bits_t ignore_bits;
    ptl_handle_me_t  me_h;
    ptl_md_t         md;
    ptl_handle_md_t  md_h;

    uint8_t          last_op;
    ptl_op_state_t   op_state;
    uint8_t          is_last_op_complete;
} portals_memory_handle;


#define NUM_REQ_QUEUES 2
typedef struct portals_request_queue_handle {
    NNTI_buffer_t *reg_buf;

    /* incoming queues */
    char *req_queue[NUM_REQ_QUEUES];

    /* each message is no larger than req_size */
    int req_size;

    /* each md can recv reqs_per_queue messages */
    int reqs_per_queue;

    /* keep track of the queue index and count */
    int indices[NUM_REQ_QUEUES];
    int queue_count[NUM_REQ_QUEUES];

    /* for each queue, we need these structs */
    ptl_handle_me_t  me_h[NUM_REQ_QUEUES];
    ptl_md_t         md[NUM_REQ_QUEUES];
    ptl_handle_md_t  md_h[NUM_REQ_QUEUES];
} portals_request_queue_handle;

typedef struct portals_transport_global {

    ptl_handle_ni_t  ni_h;

    ptl_process_id_t me;

    portals_request_queue_handle req_queue;

} portals_transport_global;



static nthread_mutex_t nnti_ptl_lock;


int process_event(
        const NNTI_buffer_t  *reg_buf,
        const NNTI_buf_ops_t  remote_op,
        const ptl_event_t    *event);
int is_buf_op_complete(
        const NNTI_buffer_t *reg_buf,
        const NNTI_buf_ops_t  remote_op);
void create_peer(
        NNTI_peer_t *peer,
        ptl_nid_t nid,
        ptl_pid_t pid);
void copy_peer(
        NNTI_peer_t *src,
        NNTI_peer_t *dest);






static portals_transport_global transport_global_data;
static const int MIN_TIMEOUT = 1000;  /* in milliseconds */

/**
 * @brief Initialize NNTI to use a specific transport.
 *
 * Enable the use of a particular transport by this process.  <tt>my_url</tt>
 * allows the process to have some control (if possible) over the
 * URL assigned for the transport.  For example, a Portals URL to put
 * might be "ptl://-1,128".  This would tell Portals to use the default
 * network ID, but use PID=128.  If the transport
 * can be initialized without this info (eg. a Portals client), <tt>my_url</tt> can
 * be NULL or empty.
 */
int NNTI_ptl_init (
        const NNTI_transport_id_t  trans_id,
        const char                *my_url,
        NNTI_transport_t          *trans_hdl)
{
    int rc=NNTI_OK;

    static uint8_t initialized=FALSE;
    int max_interfaces;
    ptl_ni_limits_t actual;


    char transport[NNTI_URL_LEN];
    char address[NNTI_URL_LEN];
    char memdesc[NNTI_URL_LEN];
    char *sep;

    NNTI_nid nid;
    NNTI_pid pid;

    assert(trans_hdl);


    if (!initialized) {

        nthread_mutex_init(&nnti_ptl_lock, NTHREAD_MUTEX_NORMAL);

        if (my_url != NULL) {
            if ((rc=nnti_url_get_transport(my_url, transport, NNTI_URL_LEN)) != NNTI_OK) {
                return(rc);
            }
            if (0!=strcmp(transport, "ptl")) {
                return(NNTI_EINVAL);
            }

            if ((rc=nnti_url_get_address(my_url, address, NNTI_URL_LEN)) != NNTI_OK) {
                return(rc);
            }

//            if ((rc=nnti_url_get_memdesc(my_url, memdesc, NNTI_URL_LEN)) != NNTI_OK) {
//                return(rc);
//            }

            sep=strchr(address, ':');
//            nid=strtol(address, NULL, 0);
            pid=strtol(sep+1, NULL, 0);
        } else {
            pid=PTL_PID_ANY;
#ifdef HAVE_TRIOS_MPI
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            pid = 128 + rank;
#endif
        }


        log_debug(nnti_debug_level, "initializing portals");

//        /* The UTCP NAL requires that the application defines where the Portals
//         * API and library should send any output. We'll send the output to
//         * the same file as logger.
//         */
//        utcp_api_out = logger_get_file();
//        utcp_lib_out = logger_get_file();
//
//
//        /* register trace groups (let someone else enable) */
//        trace_counter_gid = trace_get_gid(TRACE_RPC_COUNTER_GNAME);
//        trace_interval_gid = trace_get_gid(TRACE_RPC_INTERVAL_GNAME);


        memset(&transport_global_data, 0, sizeof(portals_transport_global));
//        transport_global_data.req_queue.eq_h=PTL_EQ_NONE;

        /* initialize the portals library */
        log_debug(nnti_debug_level, "initializing portals library");
        rc = PtlInit(&max_interfaces);
        if (rc) {
            log_fatal(nnti_debug_level,"PtlInit() failed, %s", ptl_err_str[rc]);
            abort();
        }

        /* initialize the portals interface */
        log_debug(nnti_debug_level, "initializing portals network interface - pid=%d", (int)pid);
        rc = PtlNIInit(PTL_IFACE_DEFAULT, pid, NULL, &actual, &transport_global_data.ni_h);
        if ((rc != PTL_OK) && (rc != PTL_IFACE_DUP)) {
            log_fatal(nnti_debug_level, "PtlNIInit() failed, %s", ptl_err_str[rc]);
            abort();
        }

        ptl_process_id_t ptl_id;
        rc = PtlGetId(transport_global_data.ni_h, &ptl_id);
        if (rc != PTL_OK) {
            log_error(nnti_debug_level,"failed %s", ptl_err_str[rc]);
        }

        transport_global_data.me.nid = ptl_id.nid;
        transport_global_data.me.pid = ptl_id.pid;

        if (logging_info(nnti_debug_level)) {
            fprintf(logger_get_file(), "Portals Initialized: nid=%llu, pid=%llu\n",
                    (unsigned long long)transport_global_data.me.nid,
                    (unsigned long long)transport_global_data.me.pid);
        }

        log_debug(nnti_debug_level, "sizeof(trans_hdl)=%d", trans_hdl);

        create_peer(&trans_hdl->me, transport_global_data.me.nid, transport_global_data.me.pid);

        initialized = TRUE;
    }



    return(rc);
}


/**
 * @brief Return the URL field of this transport.
 *
 * Return the URL field of this transport.  After initialization, the transport will
 * have a specific location on the network where peers can contact it.  The
 * transport will convert this location to a string that other instances of the
 * transport will recognize.
 *
 * URL format: "transport://address/memory_descriptor"
 *    - transport - (required) identifies how the URL should parsed
 *    - address   - (required) uniquely identifies a location on the network
 *                - ex. "ptl://nid:pid/", "ib://ip_addr:port", "luc://endpoint_id/"
 *    - memory_descriptor - (optional) transport-specific representation of RMA params
 */
int NNTI_ptl_get_url (
        const NNTI_transport_t *trans_hdl,
        char                   *url,
        const uint64_t          maxlen)
{
    int rc=NNTI_OK;

    assert(trans_hdl);
    assert(url);
    assert(maxlen>0);

    strncpy(url, trans_hdl->me.url, maxlen);
    url[maxlen-1]='\0';

    return(rc);
}


/**
 * @brief Prepare for communication with the peer identified by <tt>url</tt>.
 *
 * Parse <tt>url</tt> in a transport specific way.  Perform any transport specific
 * actions necessary to begin communication with this peer.
 *
// * If the peer is found and responds
// * to a ping, a handle will be allocated and assigned to the pointer.  This
// * handle should be used to move data to/from the peer.
 *
 * Connectionless transport: parse and populate
 * Connected transport: parse, connection and populate
 *
 */
int NNTI_ptl_connect (
        const NNTI_transport_t *trans_hdl,
        const char             *url,
        const int               timeout,
        NNTI_peer_t            *peer_hdl)
{
    int rc=NNTI_OK;

    char transport[NNTI_URL_LEN];
    char address[NNTI_URL_LEN];
    char memdesc[NNTI_URL_LEN];
    char *sep;

    NNTI_nid nid;
    NNTI_pid pid;

    assert(trans_hdl);
    assert(peer_hdl);

    if (url != NULL) {
        if ((rc=nnti_url_get_transport(url, transport, NNTI_URL_LEN)) != NNTI_OK) {
            return(rc);
        }
        if (0!=strcmp(transport, "ptl")) {
            /* the peer described by 'url' is not a Portals peer */
            return(NNTI_EINVAL);
        }

        if ((rc=nnti_url_get_address(url, address, NNTI_URL_LEN)) != NNTI_OK) {
            return(rc);
        }

//        if ((rc=nnti_url_get_memdesc(url, memdesc, NNTI_URL_LEN)) != NNTI_OK) {
//            return(rc);
//        }

        sep=strchr(address, ':');
        nid=strtol(address, NULL, 0);
        pid=strtol(sep+1, NULL, 0);
    } else {
        /*  */
        return(NNTI_EINVAL);
    }

    create_peer(
            peer_hdl,
            nid,
            pid);
//    peer_hdl->url = strdup(url);
//
//    peer_hdl->peer.transport_id                         = NNTI_TRANSPORT_PORTALS;
//    peer_hdl->peer.NNTI_remote_process_t_u.portals.nid  = nid;
//    peer_hdl->peer.NNTI_remote_process_t_u.portals.pid  = pid;

    return(rc);
}


/**
 * @brief Terminate communication with this peer.
 *
 * Perform any transport specific actions necessary to end communication with
 * this peer.
 */
int NNTI_ptl_disconnect (
        const NNTI_transport_t *trans_hdl,
        NNTI_peer_t            *peer_hdl)
{
    int rc=NNTI_OK;

    assert(trans_hdl);
    assert(peer_hdl);

    free(peer_hdl->url);

    return(rc);
}


/**
 * @brief Prepare a block of memory for network operations.
 *
 * Wrap a user allocated block of memory in an NNTI_buffer_t.  The transport
 * may take additional actions to prepare the memory for network send/receive.
 * If the memory block doesn't meet the transport's requirements for memory
 * regions, then errors or poor performance may result.
 */
int NNTI_ptl_register_memory (
        const NNTI_transport_t *trans_hdl,
        char                   *buffer,
        const uint64_t          size,
        const NNTI_buf_ops_t    ops,
        const NNTI_peer_t      *peer,
        NNTI_buffer_t          *reg_buf)
{
    int rc=NNTI_OK;
    static uint64_t mbits=1;

    portals_memory_handle *ptls_mem_hdl=NULL;

    assert(trans_hdl);
    assert(buffer);
    assert(size>0);
    assert(ops>0);
    assert(reg_buf);

    ptls_mem_hdl=(portals_memory_handle *)malloc(sizeof(portals_memory_handle));

    assert(ptls_mem_hdl);

    reg_buf->transport_id      = trans_hdl->id;
    reg_buf->buffer_owner      = trans_hdl->me;
    reg_buf->ops               = ops;
    reg_buf->payload_size      = size;
    reg_buf->payload           = (uint64_t)buffer;
    reg_buf->transport_private = (uint64_t)ptls_mem_hdl;
    if (peer != NULL) {
        reg_buf->peer = *peer;
    } else {
        PORTALS_SET_MATCH_ANY(&reg_buf->peer);
    }

    log_debug(nnti_debug_level, "rpc_buffer->payload_size=%ld",
            reg_buf->payload_size);

    if (ops == NNTI_RECV_DST) {
        if ((size > NNTI_REQUEST_BUFFER_SIZE) && (size%NNTI_REQUEST_BUFFER_SIZE) == 0) {
            /*
             * This is a receive-only buffer.  This buffer is divisible by
             * NNTI_REQUEST_BUFFER_SIZE.  This buffer can hold more than
             * one short request.  Assume this buffer is a request queue.
             */
            ptls_mem_hdl->type=REQUEST_BUFFER;

            ptls_mem_hdl->buffer_id   = NNTI_REQ_PT_INDEX;
            ptls_mem_hdl->match_bits  = 0;
            ptls_mem_hdl->ignore_bits = 0;
        } else if (size == NNTI_RESULT_BUFFER_SIZE) {
            /*
             * This is a receive-only buffer.  This buffer can hold exactly
             * one short result.  Assume this buffer is a result queue.
             */
            ptls_mem_hdl->type=RESULT_BUFFER;

            ptls_mem_hdl->buffer_id   = NNTI_RES_PT_INDEX;
            ptls_mem_hdl->match_bits  = mbits++;
            ptls_mem_hdl->ignore_bits = 0;
        } else {
            /*
             * This is a receive-only buffer.  This buffer doesn't look
             * like a request buffer or a result buffer.  I don't know
             * what it is.  Assume it is a regular data buffer.
             */
            ptls_mem_hdl->type=RECEIVE_BUFFER;

            ptls_mem_hdl->buffer_id   = NNTI_DATA_PT_INDEX;
            ptls_mem_hdl->match_bits  = mbits++;
            ptls_mem_hdl->ignore_bits = 0;
        }
    } else {
        if (ops == NNTI_SEND_SRC) {
            ptls_mem_hdl->type=SEND_BUFFER;
        } else if (ops == NNTI_GET_DST) {
            ptls_mem_hdl->type=GET_DST_BUFFER;
        } else if (ops == NNTI_GET_SRC) {
            ptls_mem_hdl->type=GET_SRC_BUFFER;
        } else if (ops == NNTI_PUT_SRC) {
            ptls_mem_hdl->type=PUT_SRC_BUFFER;
        } else if (ops == NNTI_PUT_DST) {
            ptls_mem_hdl->type=PUT_DST_BUFFER;
        } else if (ops == (NNTI_GET_SRC|NNTI_PUT_DST)) {
            ptls_mem_hdl->type=RDMA_TARGET_BUFFER;
        }

        ptls_mem_hdl->buffer_id   = NNTI_DATA_PT_INDEX;
        ptls_mem_hdl->match_bits  = mbits++;
        ptls_mem_hdl->ignore_bits = 0;
    }

    reg_buf->buffer_addr.transport_id                            = NNTI_TRANSPORT_PORTALS;
    reg_buf->buffer_addr.NNTI_remote_addr_t_u.portals.size       = size;
    reg_buf->buffer_addr.NNTI_remote_addr_t_u.portals.buffer_id  = (NNTI_portals_indices)ptls_mem_hdl->buffer_id;
    reg_buf->buffer_addr.NNTI_remote_addr_t_u.portals.match_bits = ptls_mem_hdl->match_bits;


    ptls_mem_hdl->match_id.nid = reg_buf->peer.peer.NNTI_remote_process_t_u.portals.nid;
    ptls_mem_hdl->match_id.pid = reg_buf->peer.peer.NNTI_remote_process_t_u.portals.pid;


    if (ptls_mem_hdl->buffer_id == NNTI_REQ_PT_INDEX) {
        uint32_t index=0;
        portals_request_queue_handle *q_hdl=&transport_global_data.req_queue;

        ptls_mem_hdl->type=REQUEST_BUFFER;
        ptls_mem_hdl->last_op=PTL_OP_NEW_REQUEST;


        q_hdl->reg_buf=reg_buf;

        q_hdl->req_size=NNTI_REQUEST_BUFFER_SIZE;
        q_hdl->reqs_per_queue=(size/q_hdl->req_size)/NUM_REQ_QUEUES;

        /* create an event queue */
        nthread_lock(&nnti_ptl_lock);
        rc = PtlEQAlloc(
                transport_global_data.ni_h,
                q_hdl->reqs_per_queue,
                PTL_EQ_HANDLER_NONE,
                &ptls_mem_hdl->eq_h);
        nthread_unlock(&nnti_ptl_lock);
        if (rc != NNTI_OK) {
            log_error(nnti_debug_level, "PtlEQAlloc() failed");
            goto cleanup;
        }
        log_debug(nnti_debug_level, "allocated eq=%d", ptls_mem_hdl->eq_h);

        /* Accept requests from anyone */
        ptls_mem_hdl->match_id.nid = PTL_NID_ANY;
        ptls_mem_hdl->match_id.pid = PTL_PID_ANY;
        ptls_mem_hdl->match_bits   = 0;
        ptls_mem_hdl->ignore_bits  = 0;


        for (index=0; index<NUM_REQ_QUEUES; index++) {
            /* initialize the indices stored in the MD user pointer */
            q_hdl->indices[index] = index;
            q_hdl->queue_count[index] = 0;

            /* allocate the buffer for the incoming MD */
            q_hdl->req_queue[index] = buffer + (index*q_hdl->reqs_per_queue*q_hdl->req_size);

            /* initialize the buffer */
            memset(q_hdl->req_queue[index], 0, q_hdl->reqs_per_queue*q_hdl->req_size);

            /* initialize the MD */
            memset(&q_hdl->md[index], 0, sizeof(ptl_md_t));
            q_hdl->md[index].start = q_hdl->req_queue[index];
            q_hdl->md[index].length = q_hdl->reqs_per_queue*q_hdl->req_size;
            q_hdl->md[index].threshold = q_hdl->reqs_per_queue;
            q_hdl->md[index].max_size = q_hdl->req_size;
            q_hdl->md[index].options = PTL_MD_OP_PUT | PTL_MD_MAX_SIZE;
            q_hdl->md[index].user_ptr = &q_hdl->indices[index];  /* used to store the index */
            q_hdl->md[index].eq_handle = ptls_mem_hdl->eq_h;

            log_debug(nnti_debug_level, "attaching match entry to index=%d",
                    ptls_mem_hdl->buffer_id);

            /* Attach the match entry to the portal index */
            nthread_lock(&nnti_ptl_lock);
            rc = PtlMEAttach(
                    transport_global_data.ni_h,
                    ptls_mem_hdl->buffer_id,
                    ptls_mem_hdl->match_id,
                    ptls_mem_hdl->match_bits,
                    ptls_mem_hdl->ignore_bits,
                    PTL_RETAIN,
                    PTL_INS_AFTER,
                    &q_hdl->me_h[index]);
            nthread_unlock(&nnti_ptl_lock);
            if (rc != NNTI_OK) {
                log_error(nnti_debug_level, "could not attach ME");
                goto cleanup;
            }

            /* Attach the MD to the match entry */
            nthread_lock(&nnti_ptl_lock);
            rc = PtlMDAttach(
                    q_hdl->me_h[index],
                    q_hdl->md[index],
                    PTL_RETAIN,
                    &q_hdl->md_h[index]);
            nthread_unlock(&nnti_ptl_lock);
            if (rc != NNTI_OK) {
                log_error(nnti_debug_level, "could not alloc eq: %s",
                        ptl_err_str[rc]);
                goto cleanup;
            }
            log_debug(nnti_debug_level, "attached q_hdl->md_h[%d]: %d", index, q_hdl->md_h[index]);

            reg_buf->payload_size=q_hdl->req_size;
        }
    } else {

        /* create an event queue */
        /* TODO: should we share an event queue? */
        nthread_lock(&nnti_ptl_lock);
        rc = PtlEQAlloc(
                transport_global_data.ni_h,
                5,
                PTL_EQ_HANDLER_NONE,
                &ptls_mem_hdl->eq_h);
        nthread_unlock(&nnti_ptl_lock);
        if (rc != NNTI_OK) {
            log_error(nnti_debug_level, "failed to allocate eventq");
            goto cleanup;
        }
        log_debug(nnti_debug_level, "allocated eq=%d", ptls_mem_hdl->eq_h);

        /* create a match entry (unlink with MD) */
        nthread_lock(&nnti_ptl_lock);
        rc = PtlMEAttach(
                transport_global_data.ni_h,
                ptls_mem_hdl->buffer_id,
                ptls_mem_hdl->match_id,
                ptls_mem_hdl->match_bits,
                ptls_mem_hdl->ignore_bits,
                PTL_UNLINK,
                PTL_INS_AFTER,
                &ptls_mem_hdl->me_h);
        nthread_unlock(&nnti_ptl_lock);
        if (rc != NNTI_OK) {
            log_error(nnti_debug_level, "failed to allocate eventq");
            goto cleanup;
        }
        log_debug(nnti_debug_level, "allocated me=%d with bufid=%d, match_id(%d,%d), mbits=%d",
                ptls_mem_hdl->me_h, ptls_mem_hdl->buffer_id, ptls_mem_hdl->match_id.nid, ptls_mem_hdl->match_id.pid, ptls_mem_hdl->match_bits);

        /* initialize the md */
        memset(&ptls_mem_hdl->md, 0, sizeof(ptl_md_t));
        ptls_mem_hdl->md.start     = buffer;
        ptls_mem_hdl->md.length    = size;
        ptls_mem_hdl->md.threshold = PTL_MD_THRESH_INF;
        ptls_mem_hdl->md.options   = PTL_MD_OP_PUT|PTL_MD_OP_GET|PTL_MD_TRUNCATE;
        ptls_mem_hdl->md.user_ptr  = NULL;
        ptls_mem_hdl->md.eq_handle = ptls_mem_hdl->eq_h;

        /* attach the memory descriptor (manually unlink) */
        nthread_lock(&nnti_ptl_lock);
        rc = PtlMDAttach(
                ptls_mem_hdl->me_h,
                ptls_mem_hdl->md,
                PTL_RETAIN,
                &ptls_mem_hdl->md_h);
        nthread_unlock(&nnti_ptl_lock);
        if (rc != NNTI_OK) {
            log_error(nnti_debug_level, "failed to attach md");
            goto cleanup;
        }
        log_debug(nnti_debug_level, "attached ptls_mem_hdl->md_h: %d", ptls_mem_hdl->md_h);
    }

cleanup:
//    log_debug(LOG_ALL, "registering buf(%lu) md_h(%d) eq_h(%d)",
//        reg_buf->payload, ptls_mem_hdl->md_h, ptls_mem_hdl->eq_h);

    return(rc);
}


/**
 * @brief Cleanup after network operations are complete.
 *
 * Destroy an NNTI_buffer_t that was previously created by NNTI_regsiter_buffer().
 * It is the user's responsibility to release the the memory region.
 */
int NNTI_ptl_unregister_memory (
        NNTI_buffer_t    *reg_buf)
{
    int rc=NNTI_OK, rc2=NNTI_OK;
    portals_memory_handle *ptls_mem_hdl=NULL;
    log_level debug_level = nnti_debug_level;

    assert(reg_buf);

    ptls_mem_hdl=(portals_memory_handle *)reg_buf->transport_private;

    assert(ptls_mem_hdl);

    log_debug(debug_level, "unregistering buf(%lu) md_h(%d) eq_h(%d)",
    		reg_buf->payload, ptls_mem_hdl->md_h, ptls_mem_hdl->eq_h);

    if (ptls_mem_hdl->buffer_id == NNTI_REQ_PT_INDEX) {
        uint32_t index=0;
        portals_request_queue_handle *q_hdl=&transport_global_data.req_queue;

        for (index=0; index<NUM_REQ_QUEUES; index++) {
            /* unlink the memory descriptor */
            log_debug(debug_level, "unlinking q_hdl->md_h[%d]: %d", index, q_hdl->md_h[index]);
            nthread_lock(&nnti_ptl_lock);
            rc2 = PtlMDUnlink(q_hdl->md_h[index]);
            nthread_unlock(&nnti_ptl_lock);
            if (rc2 != NNTI_OK) {
                log_warn(debug_level, "unable to unlink memory descriptor for queue %d: %s",
                        index, ptl_err_str[rc2]);
                rc = NNTI_EBADRPC;
                goto cleanup;
            }
        }

        /* free the event queue */
        log_debug(debug_level, "freeing ptls_mem_hdl->eq_h: %d", ptls_mem_hdl->eq_h);
        nthread_lock(&nnti_ptl_lock);
        rc2 = PtlEQFree(ptls_mem_hdl->eq_h);
        nthread_unlock(&nnti_ptl_lock);
        if (rc2 != NNTI_OK) {
            log_fatal(debug_level, "unable to free EQ: %s", ptl_err_str[rc2]);
            rc = NNTI_EBADRPC;
            goto cleanup;
        }
    } else {
        log_debug(debug_level, "unlinking ptls_mem_hdl->md_h: %d", ptls_mem_hdl->md_h);
        nthread_lock(&nnti_ptl_lock);
        rc = PtlMDUnlink(ptls_mem_hdl->md_h);
        nthread_unlock(&nnti_ptl_lock);
        if (rc != PTL_OK) {
            log_error(debug_level, "failed to unlink MD: %s", ptl_err_str[rc]);
            goto cleanup;
        }

        log_debug(debug_level, "freeing ptls_mem_hdl->eq_h: %d", ptls_mem_hdl->eq_h);
        nthread_lock(&nnti_ptl_lock);
        rc = PtlEQFree(ptls_mem_hdl->eq_h);
        nthread_unlock(&nnti_ptl_lock);
        if (rc != PTL_OK) {
            log_error(debug_level, "failed to free EQ: %s", ptl_err_str[rc]);
            goto cleanup;
        }
    }


cleanup:

    //  This was caught by valgrind. Allocated in NNTI_ptl_register_memory
    if (ptls_mem_hdl) free (ptls_mem_hdl);

    reg_buf->transport_id      = NNTI_TRANSPORT_NULL;
    PORTALS_SET_MATCH_ANY(&reg_buf->buffer_owner);
    reg_buf->ops               = (NNTI_buf_ops_t)0;
    PORTALS_SET_MATCH_ANY(&reg_buf->peer);
    reg_buf->payload_size      = 0;
    reg_buf->payload           = 0;
    reg_buf->transport_private = 0;

    log_debug(debug_level, "Finished unregistering, rc=%d",rc);

    return(rc);
}


/**
 * @brief Send a message to a peer.
 *
 * Send a message (<tt>msg_hdl</tt>) to a peer (<tt>peer_hdl</tt>).  It is expected that the
 * message is small, but the exact maximum size is transport dependent.
 */
int NNTI_ptl_send (
        const NNTI_peer_t   *peer_hdl,
        const NNTI_buffer_t *msg_hdl,
        const NNTI_buffer_t *dest_hdl)
{
    int rc=NNTI_OK;

    portals_memory_handle *ptls_mem_hdl=NULL;
    ptl_process_id_t dest_id;
    ptl_pt_index_t   buffer_id;
    ptl_match_bits_t match_bits;

    assert(peer_hdl);
    assert(msg_hdl);

    ptls_mem_hdl=(portals_memory_handle *)msg_hdl->transport_private;

    if (dest_hdl == NULL) {
        dest_id.nid = peer_hdl->peer.NNTI_remote_process_t_u.portals.nid;
        dest_id.pid = peer_hdl->peer.NNTI_remote_process_t_u.portals.pid;
        buffer_id   = NNTI_REQ_PT_INDEX;
        match_bits  = 0;
    } else {
        dest_id.nid = dest_hdl->buffer_owner.peer.NNTI_remote_process_t_u.portals.nid;
        dest_id.pid = dest_hdl->buffer_owner.peer.NNTI_remote_process_t_u.portals.pid;
        buffer_id   = dest_hdl->buffer_addr.NNTI_remote_addr_t_u.portals.buffer_id;
        match_bits  = dest_hdl->buffer_addr.NNTI_remote_addr_t_u.portals.match_bits;
    }

    log_debug(nnti_debug_level, "sending to (nid=%d, pid=%d, buffer_id=%d, mbits=%d)", dest_id.nid, dest_id.pid, buffer_id, match_bits);

    rc=PtlPut(
            ptls_mem_hdl->md_h,
            PTL_ACK_REQ,
            dest_id,
            buffer_id,
            0,
            match_bits,
            0,
            0);

    ptls_mem_hdl->last_op=PTL_OP_SEND;

    return(rc);
}


/**
 * @brief Transfer data to a peer.
 *
 * Put the contents of <tt>src_buffer_hdl</tt> into <tt>dest_buffer_hdl</tt>.  It is
 * assumed that the destination is at least <tt>src_length</tt> bytes in size.
 *
 */
int NNTI_ptl_put (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset)
{
    int rc=NNTI_OK;

    portals_memory_handle *ptls_mem_hdl=NULL;
    ptl_process_id_t dest_id;

    assert(src_buffer_hdl);
    assert(dest_buffer_hdl);

    ptls_mem_hdl=(portals_memory_handle *)src_buffer_hdl->transport_private;

    dest_id.nid=dest_buffer_hdl->buffer_owner.peer.NNTI_remote_process_t_u.portals.nid;
    dest_id.pid=dest_buffer_hdl->buffer_owner.peer.NNTI_remote_process_t_u.portals.pid;

    rc=PtlPutRegion(
            ptls_mem_hdl->md_h,
            src_offset,
            src_length,
            PTL_ACK_REQ,
            dest_id,
            dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.portals.buffer_id,
            0,
            dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.portals.match_bits,
            dest_offset,
            0);

    ptls_mem_hdl->last_op=PTL_OP_PUT_INITIATOR;

    return(rc);
}


/**
 * @brief Transfer data from a peer.
 *
 * Get the contents of <tt>src_buffer_hdl</tt> into <tt>dest_buffer_hdl</tt>.  It is
 * assumed that the destination is at least <tt>src_length</tt> bytes in size.
 *
 */
int NNTI_ptl_get (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset)
{
    int rc=NNTI_OK;

    portals_memory_handle *ptls_mem_hdl=NULL;
    ptl_process_id_t src_id;

    assert(src_buffer_hdl);
    assert(dest_buffer_hdl);

    ptls_mem_hdl=(portals_memory_handle *)dest_buffer_hdl->transport_private;

    src_id.nid=src_buffer_hdl->buffer_owner.peer.NNTI_remote_process_t_u.portals.nid;
    src_id.pid=src_buffer_hdl->buffer_owner.peer.NNTI_remote_process_t_u.portals.pid;

    rc=PtlGetRegion(
            ptls_mem_hdl->md_h,
            dest_offset,
            src_length,
            src_id,
            src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.portals.buffer_id,
            0,
            src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.portals.match_bits,
            src_offset);

    ptls_mem_hdl->last_op=PTL_OP_GET_INITIATOR;

    return(rc);
}


/**
 * @brief Wait for <tt>remote_op</tt> on <tt>reg_buf</tt> to complete.
 *
 * Wait for <tt>remote_op</tt> on <tt>reg_buf</tt> to complete or timeout
 * waiting.  This is typically used to wait for a result or a bulk data
 * transfer.  The timeout is specified in milliseconds.  A timeout of <tt>-1</tt>
 * means wait forever.  A timeout of <tt>0</tt> means do not wait.
 *
 */
int NNTI_ptl_wait (
        const NNTI_buffer_t  *reg_buf,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        NNTI_status_t        *status)
{
    int nnti_rc=NNTI_OK;
    portals_memory_handle *ptls_mem_hdl=NULL;

    int rc = PTL_EQ_EMPTY;
    int elapsed_time = 0;
    int timeout_per_call;
    ptl_event_t event;
    int which;


    assert(reg_buf);
    assert(status);

    ptls_mem_hdl=(portals_memory_handle *)reg_buf->transport_private;

    assert(ptls_mem_hdl);

    if (timeout < 0)
        timeout_per_call = MIN_TIMEOUT;
    else
        timeout_per_call = (timeout < MIN_TIMEOUT)? MIN_TIMEOUT : timeout;

    memset(&ptls_mem_hdl->op_state, 0, sizeof(ptl_op_state_t));

    while (1)   {
        if (trios_exit_now()) {
            log_debug(nnti_debug_level, "caught abort signal");
            return NNTI_ECANCELED;
        }

        log_debug(nnti_debug_level, "waiting on eq_h(%d)", ptls_mem_hdl->eq_h);

        memset(&event, 0, sizeof(ptl_event_t));
        log_debug(nnti_debug_level, "lock before poll");
//        nthread_lock(&nnti_ptl_lock);
        rc = PtlEQPoll(&ptls_mem_hdl->eq_h, 1, timeout_per_call, &event, &which);
//        nthread_lock(&nnti_ptl_lock);
        log_debug(nnti_debug_level, "polling status is %s", ptl_err_str[rc]);

        log_debug(nnti_debug_level, "Poll Event= {");
        log_debug(nnti_debug_level, "\ttype         = %d", event.type);
        log_debug(nnti_debug_level, "\tinitiator    = (%llu, %llu)", (unsigned long long)event.initiator.nid, (unsigned long long)event.initiator.pid);
        log_debug(nnti_debug_level, "\tuid          = %d", event.uid);
        log_debug(nnti_debug_level, "\tjid          = %d", event.jid);
        log_debug(nnti_debug_level, "\tpt_index     = %d", event.pt_index);
        log_debug(nnti_debug_level, "\tmatch_bits   = %d", event.match_bits);
        log_debug(nnti_debug_level, "\trlength      = %llu", (unsigned long long)event.rlength);
        log_debug(nnti_debug_level, "\tmlength      = %llu", (unsigned long long)event.mlength);
        log_debug(nnti_debug_level, "\toffset       = %llu", (unsigned long long)event.offset);
        log_debug(nnti_debug_level, "\tmd_handle    = %d", event.md_handle);
        log_debug(nnti_debug_level, "\tmd.start     = %p", event.md.start);
        log_debug(nnti_debug_level, "\tmd.length    = %d", event.md.length);
        log_debug(nnti_debug_level, "\tmd.max_size  = %d", event.md.max_size);
        log_debug(nnti_debug_level, "\tmd.threshold = %d", event.md.threshold);


        /* case 1: success */
        if (rc == PTL_OK) {
            nnti_rc = NNTI_OK;
        }
        /* case 2: success, but some events were dropped */
        else if (rc == PTL_EQ_DROPPED) {
            log_warn(nnti_debug_level, "PtlEQPoll dropped some events");
            log_warn(nnti_debug_level, "PtlEQPoll succeeded, but at least one event was dropped");
            nnti_rc = NNTI_OK;
        }
        /* case 3: timed out */
        else if (rc == PTL_EQ_EMPTY) {
            elapsed_time += timeout_per_call;

            /* if the caller asked for a legitimate timeout, we need to exit */
            if (((timeout > 0) && (elapsed_time >= timeout))) {
                log_debug(nnti_debug_level, "PtlEQPoll timed out: %s",
                        ptl_err_str[rc]);
                nnti_rc = NNTI_ETIMEDOUT;
                break;
            }
            /* continue if the timeout has not expired */
            /* log_debug(nnti_debug_level, "timedout... continuing"); */



            continue;
        }
        /* case 4: failure */
        else {
            log_error(nnti_debug_level, "PtlEQPoll failed (eq_handle[%d]==%d): %s",
                    which, ptls_mem_hdl->eq_h, ptl_err_str[rc]);
            nnti_rc = NNTI_EIO;
            break;
        }

        process_event(reg_buf, remote_op, &event);
        if (is_buf_op_complete(reg_buf, remote_op) == TRUE) {
            break;
        }
    }

    if ((rc!=0) && (event.ni_fail_type != PTL_NI_OK)) {
        log_error(nnti_debug_level, "NI reported error: ni_fail_type=%s",
                PtlNIFailStr(transport_global_data.ni_h, event.ni_fail_type));
        nnti_rc = NNTI_EIO;
    }

    status->op     = remote_op;
    status->start  = (uint64_t)event.md.start;
    status->offset = event.offset;
    status->length = event.mlength;
    status->result = (NNTI_result_t)nnti_rc;
    switch (ptls_mem_hdl->last_op) {
        case PTL_OP_PUT_INITIATOR:
        case PTL_OP_GET_TARGET:
        case PTL_OP_SEND:
            create_peer(&status->src, transport_global_data.me.nid, transport_global_data.me.pid); // allocates url
            create_peer(&status->dest, event.initiator.nid, event.initiator.pid); // allocates url
            break;
        case PTL_OP_GET_INITIATOR:
        case PTL_OP_PUT_TARGET:
        case PTL_OP_NEW_REQUEST:
        case PTL_OP_RESULT:
        case PTL_OP_RECEIVE:
            create_peer(&status->src, event.initiator.nid, event.initiator.pid); // allocates url
            create_peer(&status->dest, transport_global_data.me.nid, transport_global_data.me.pid); // allocates url
            break;
    }

    if ((nnti_rc==NNTI_OK) && (ptls_mem_hdl->buffer_id == NNTI_REQ_PT_INDEX)) {
        portals_request_queue_handle *q_hdl=&transport_global_data.req_queue;

        int index = *(int *)event.md.user_ptr;
        /* get the index of the queue */
        q_hdl->queue_count[index]++;

        /* if we've processed all we can on this queue, reset */
        if (q_hdl->queue_count[index] >= q_hdl->reqs_per_queue) {

            log_debug(LOG_ALL, "Resetting MD on queue[%d]", index);

            /* Unlink the ME (also unlinks the MD) */
            nthread_lock(&nnti_ptl_lock);
            rc = PtlMEUnlink(q_hdl->me_h[index]);
            nthread_unlock(&nnti_ptl_lock);
            if (rc != PTL_OK) {
                log_error(nnti_debug_level, "Could not unlink ME: %s", ptl_err_str[rc]);
                goto cleanup;
            }

            /* Re-attach the match-list entry */
            nthread_lock(&nnti_ptl_lock);
            rc = PtlMEAttach(
                    transport_global_data.ni_h,
                    ptls_mem_hdl->buffer_id,
                    ptls_mem_hdl->match_id,
                    ptls_mem_hdl->match_bits,
                    ptls_mem_hdl->ignore_bits,
                    PTL_RETAIN,
                    PTL_INS_AFTER,
                    &q_hdl->me_h[index]);
            nthread_unlock(&nnti_ptl_lock);
            if (rc != PTL_OK) {
                log_error(nnti_debug_level, "Could not reset ME: %s", ptl_err_str[rc]);
                goto cleanup;
            }

            /* Re-attach the MD */
            nthread_lock(&nnti_ptl_lock);
            rc = PtlMDAttach(
                    q_hdl->me_h[index],
                    q_hdl->md[index],
                    PTL_RETAIN,
                    &q_hdl->md_h[index]);
            nthread_unlock(&nnti_ptl_lock);
            if (rc != PTL_OK) {
                log_error(nnti_debug_level, "Could not reset MD: %s", ptl_err_str[rc]);
                goto cleanup;
            }

            q_hdl->queue_count[index] = 0;
        }
    }

cleanup:
    return(nnti_rc);
}


/**
 * @brief Disable this transport.
 *
 * Shutdown the transport.  Any outstanding sends, gets and puts will be
 * canceled.  Any new transport requests will fail.
 *
 */
int NNTI_ptl_fini (
        const NNTI_transport_t *trans_hdl)
{
    PtlFini();

    return(NNTI_OK);
}









int process_event(
        const NNTI_buffer_t  *reg_buf,
        const NNTI_buf_ops_t  remote_op,
        const ptl_event_t    *event)
{
    int rc=NNTI_OK;
    portals_memory_handle *ptls_mem_hdl=NULL;
    log_level debug_level = nnti_debug_level;

    ptls_mem_hdl=(portals_memory_handle *)reg_buf->transport_private;

    log_debug(debug_level, "reg_buf=%p; ptls_mem_hdl->last_op=%d; remote_op=%d", reg_buf, ptls_mem_hdl->last_op, remote_op);
    switch (ptls_mem_hdl->type) {
        case SEND_BUFFER:
        case PUT_SRC_BUFFER:
            switch (event->type) {
                case PTL_EVENT_SEND_START:
                    log_debug(debug_level, "got PTL_EVENT_SEND_START - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptls_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptls_mem_hdl->op_state.put_initiator.send_start = TRUE;
                    break;
                case PTL_EVENT_SEND_END:
                    log_debug(debug_level, "got PTL_EVENT_SEND_END   - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptls_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptls_mem_hdl->op_state.put_initiator.send_end = TRUE;
                    break;
                case PTL_EVENT_ACK:
                    log_debug(debug_level, "got PTL_EVENT_ACK        - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptls_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptls_mem_hdl->op_state.put_initiator.ack = TRUE;
                    break;
                default:
                    log_error(debug_level, "unrecognized event type: %d - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            event->type, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    rc = NNTI_EINVAL;
                    goto cleanup;
            }
            break;
        case GET_DST_BUFFER:
            switch (event->type) {
                case PTL_EVENT_SEND_START:
                    log_debug(debug_level, "got PTL_EVENT_SEND_START - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptls_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptls_mem_hdl->op_state.get_initiator.send_start = TRUE;
                    break;
                case PTL_EVENT_SEND_END:
                    log_debug(debug_level, "got PTL_EVENT_SEND_END   - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptls_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptls_mem_hdl->op_state.get_initiator.send_end = TRUE;
                    break;
                case PTL_EVENT_REPLY_START:
                    log_debug(debug_level,"got PTL_EVENT_REPLY_START - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptls_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptls_mem_hdl->op_state.get_initiator.reply_start = TRUE;
                    break;
                case PTL_EVENT_REPLY_END:
                    log_debug(debug_level,"got PTL_EVENT_REPLY_END   - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptls_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptls_mem_hdl->op_state.get_initiator.reply_end = TRUE;
                    break;
                default:
                    log_error(debug_level, "unrecognized event type: %d - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            event->type, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    rc = NNTI_EINVAL;
                    goto cleanup;
            }
            break;
        case REQUEST_BUFFER:
        case RESULT_BUFFER:
        case RECEIVE_BUFFER:
            switch (event->type) {
                case PTL_EVENT_PUT_START:
                    log_debug(debug_level, "got PTL_EVENT_PUT_START  - new request - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptls_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);

                    break;
                case PTL_EVENT_PUT_END:
                    log_debug(debug_level, "got PTL_EVENT_PUT_END    - new request - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptls_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptls_mem_hdl->op_state.put_target.put_start = TRUE;
                    ptls_mem_hdl->op_state.put_target.put_end = TRUE;
                    break;
                default:
                    log_error(debug_level, "unrecognized event type: %d - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            event->type, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    rc = NNTI_EINVAL;
                    goto cleanup;
            }
            break;
        case PUT_DST_BUFFER:
            switch (event->type) {
                case PTL_EVENT_PUT_START:
                    log_debug(debug_level, "got PTL_EVENT_PUT_START  - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptls_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptls_mem_hdl->op_state.put_target.put_start = TRUE;
                    break;
                case PTL_EVENT_PUT_END:
                    log_debug(debug_level, "got PTL_EVENT_PUT_END    - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptls_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptls_mem_hdl->op_state.put_target.put_end = TRUE;
                    break;
                default:
                    log_error(debug_level, "unrecognized event type: %d - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            event->type, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    rc = NNTI_EINVAL;
                    goto cleanup;
            }
            break;
        case GET_SRC_BUFFER:
            switch (event->type) {
                case PTL_EVENT_GET_START:
                    log_debug(debug_level, "got PTL_EVENT_GET_START  - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptls_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptls_mem_hdl->op_state.get_target.get_start = TRUE;
                    break;
                case PTL_EVENT_GET_END:
                    log_debug(debug_level, "got PTL_EVENT_GET_END    - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptls_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptls_mem_hdl->op_state.get_target.get_end = TRUE;
                    break;
                default:
                    log_error(debug_level, "unrecognized event type: %d - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            event->type, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    rc = NNTI_EINVAL;
                    goto cleanup;
            }
            break;
            case RDMA_TARGET_BUFFER:
                switch (event->type) {
                    case PTL_EVENT_PUT_START:
                        log_debug(debug_level, "got PTL_EVENT_PUT_START  - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                                ptls_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                        ptls_mem_hdl->op_state.put_target.put_start = TRUE;
                        break;
                    case PTL_EVENT_PUT_END:
                        log_debug(debug_level, "got PTL_EVENT_PUT_END    - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                                ptls_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                        ptls_mem_hdl->op_state.put_target.put_end = TRUE;
                        break;
                    case PTL_EVENT_GET_START:
                        log_debug(debug_level, "got PTL_EVENT_GET_START  - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                                ptls_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                        ptls_mem_hdl->op_state.get_target.get_start = TRUE;
                        break;
                    case PTL_EVENT_GET_END:
                        log_debug(debug_level, "got PTL_EVENT_GET_END    - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                                ptls_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                        ptls_mem_hdl->op_state.get_target.get_end = TRUE;
                        break;
                    default:
                        log_error(debug_level, "unrecognized event type: %d - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                                event->type, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                        rc = NNTI_EINVAL;
                        goto cleanup;
                }
                break;
    }

    if (event->ni_fail_type != PTL_NI_OK) {
        log_error(debug_level, "failed on put end: ni_fail_type=%d\n",
                event->ni_fail_type);
        rc = event->ni_fail_type;
    }

cleanup:
    return (rc);
}

int is_buf_op_complete(
        const NNTI_buffer_t *reg_buf,
        const NNTI_buf_ops_t  remote_op)
{
    int rc=FALSE;
    portals_memory_handle *ptls_mem_hdl=NULL;
    log_level debug_level = nnti_debug_level;

    ptls_mem_hdl=(portals_memory_handle *)reg_buf->transport_private;

    switch (ptls_mem_hdl->type) {
        case SEND_BUFFER:
        case PUT_SRC_BUFFER:
            if ((ptls_mem_hdl->op_state.put_initiator.send_start==TRUE) &&
                (ptls_mem_hdl->op_state.put_initiator.send_end==TRUE)   &&
                (ptls_mem_hdl->op_state.put_initiator.ack==TRUE))       {
                ptls_mem_hdl->last_op=PTL_OP_PUT_INITIATOR;
                rc = TRUE;
            }
            break;
        case GET_DST_BUFFER:
            /* cray portals */
            if ((ptls_mem_hdl->op_state.get_initiator.send_start==TRUE)  &&
                (ptls_mem_hdl->op_state.get_initiator.send_end==TRUE)    &&
                (ptls_mem_hdl->op_state.get_initiator.reply_start==TRUE) &&
                (ptls_mem_hdl->op_state.get_initiator.reply_end==TRUE))  {
                ptls_mem_hdl->last_op=PTL_OP_GET_INITIATOR;
                rc = TRUE;
                break;
            }
            /* schutt portals */
            if ((ptls_mem_hdl->op_state.get_initiator.reply_start==TRUE) &&
                (ptls_mem_hdl->op_state.get_initiator.reply_end==TRUE))  {
                ptls_mem_hdl->last_op=PTL_OP_GET_INITIATOR;
                rc = TRUE;
                break;
            }
            break;
        case PUT_DST_BUFFER:
            if ((ptls_mem_hdl->op_state.put_target.put_start==TRUE) &&
                (ptls_mem_hdl->op_state.put_target.put_end==TRUE))  {
                ptls_mem_hdl->last_op=PTL_OP_PUT_TARGET;
                rc = TRUE;
            }
            break;
        case GET_SRC_BUFFER:
            if ((ptls_mem_hdl->op_state.get_target.get_start==TRUE) &&
                (ptls_mem_hdl->op_state.get_target.get_end==TRUE))  {
                ptls_mem_hdl->last_op=PTL_OP_GET_TARGET;
                rc = TRUE;
            }
            break;
        case REQUEST_BUFFER:
            if ((ptls_mem_hdl->op_state.put_target.put_start==TRUE) &&
                (ptls_mem_hdl->op_state.put_target.put_end==TRUE)) {
                ptls_mem_hdl->last_op=PTL_OP_NEW_REQUEST;
                rc = TRUE;
            }
            break;
        case RESULT_BUFFER:
            if ((ptls_mem_hdl->op_state.put_target.put_start==TRUE) &&
                (ptls_mem_hdl->op_state.put_target.put_end==TRUE)) {
                ptls_mem_hdl->last_op=PTL_OP_RESULT;
                rc = TRUE;
            }
            break;
        case RECEIVE_BUFFER:
            if ((ptls_mem_hdl->op_state.put_target.put_start==TRUE) &&
                (ptls_mem_hdl->op_state.put_target.put_end==TRUE)) {
                ptls_mem_hdl->last_op=PTL_OP_RECEIVE;
                rc = TRUE;
            }
            break;
        case RDMA_TARGET_BUFFER:
            if ((ptls_mem_hdl->op_state.get_target.get_start==TRUE) &&
                (ptls_mem_hdl->op_state.get_target.get_end==TRUE))  {
                ptls_mem_hdl->last_op=PTL_OP_GET_TARGET;
                rc = TRUE;
            }
            if ((ptls_mem_hdl->op_state.put_target.put_start==TRUE) &&
                (ptls_mem_hdl->op_state.put_target.put_end==TRUE)) {
                ptls_mem_hdl->last_op=PTL_OP_PUT_TARGET;
                rc = TRUE;
            }
            break;
    }
    if (rc==TRUE) {
        log_debug(nnti_debug_level, "op is complete");
    }

    return(rc);
}

void create_peer(NNTI_peer_t *peer, ptl_nid_t nid, ptl_pid_t pid)
{
    //peer->url = (char *)malloc(NNTI_URL_LEN);
    log_debug(nnti_debug_level, "sizeof(peer->url)=%d", sizeof(peer->url));
    sprintf(peer->url, "ptl://%u:%u/", nid, pid);

    peer->peer.transport_id                        = NNTI_TRANSPORT_PORTALS;
    peer->peer.NNTI_remote_process_t_u.portals.nid = nid;
    peer->peer.NNTI_remote_process_t_u.portals.pid = pid;
}

void copy_peer(NNTI_peer_t *src, NNTI_peer_t *dest)
{
    //dest->url = strdup(src->url);
    strcpy(dest->url, src->url);

    dest->peer.transport_id                        = NNTI_TRANSPORT_PORTALS;
    dest->peer.NNTI_remote_process_t_u.portals.nid = src->peer.NNTI_remote_process_t_u.portals.nid;
    dest->peer.NNTI_remote_process_t_u.portals.pid = src->peer.NNTI_remote_process_t_u.portals.nid;
}
