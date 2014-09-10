/**
//@HEADER
// ************************************************************************
//
//                   Trios: Trilinos I/O Support
//                 Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)
//
// *************************************************************************
//@HEADER
 */
/**
 * nnti_gni.c
 *
 *  Created on: Jan 13, 2011
 *      Author: thkorde
 */

#include "Trios_config.h"
#include "Trios_threads.h"
#include "Trios_timer.h"
#include "Trios_signal.h"
#include "Trios_nnti_fprint_types.h"

#include <assert.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/poll.h>
#include <sys/mman.h>
#include <ifaddrs.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netdb.h>

#include <sched.h>
#include <alps/libalpslli.h>
#include <gni_pub.h>

#include <map>
#include <deque>
#include <algorithm>

#ifdef HAVE_TRIOS_HPCTOOLKIT
#include <hpctoolkit.h>
#define SAMPLING_IS_ACTIVE() hpctoolkit_sampling_is_active()
#define SAMPLING_STOP() hpctoolkit_sampling_stop()
#define SAMPLING_START() hpctoolkit_sampling_start()
#else
#define SAMPLING_IS_ACTIVE() 0
#define SAMPLING_STOP()
#define SAMPLING_START()
#endif

#include "nnti_gni.h"
#include "nnti_utils.h"



///* if undefined, the ACK message is NOT sent to the RDMA target when
// * the RDMA op is complete.  this creates one-side semantics for RDMA
// * ops.  in this mode, the target has no idea when the RDMA op is
// * complete and what data was addressed.
// */
//#undef USE_RDMA_TARGET_ACK
///* if defined, the RDMA initiator will send an ACK message to the RDMA
// * target when the RDMA op is complete.  the target process must wait
// * on the target buffer in order to get the ACK.  this creates two-side
// * semantics for RDMA ops.   in this mode, when the wait returns the
// * the RDMA op is complete and status indicates what data was addressed.
// */
//#define USE_RDMA_TARGET_ACK



//#undef USE_ANONYMOUS_MEMORY_REGISTRATION
//#define USE_ANONYMOUS_MEMORY_REGISTRATION

//#undef USE_RDMA_EVENTS
//#define USE_RDMA_EVENTS

// NSSI
//#define USE_RDMA_TARGET_ACK
//#define USE_ANONYMOUS_MEMORY_REGISTRATION
//#define USE_RDMA_EVENTS
//#undef USE_RDMA_FENCE

// ORNL
//#undef USE_RDMA_TARGET_ACK
//#define USE_ANONYMOUS_MEMORY_REGISTRATION
//#undef USE_RDMA_EVENTS
//#undef USE_RDMA_FENCE


//#define USE_FMA
//#define USE_BTE
//#define USE_CROSSOVER
//#define USE_MIXED

//#define FMA_BTE_CROSSOVER 4096

///* mode 1 - client uses GNI params from ALPS to create a Commumnication Domain and attach to it */
//#define USE_ALPS_PTAG
///* mode 2 - client uses a mix of GNI params from ALPS and the server (cookie/pTag) to create a Commumnication Domain and attach to it */
////#undef USE_ALPS_PTAG



#define NIC_ADDR_BITS    22
#define NIC_ADDR_SHIFT   (32-NIC_ADDR_BITS)
#define NIC_ADDR_MASK    0x3FFFFF
#define CPU_NUM_BITS     7
#define CPU_NUM_SHIFT    (NIC_ADDR_SHIFT-CPU_NUM_BITS)
#define CPU_NUM_MASK     0x7F
#define THREAD_NUM_BITS  3
#define THREAD_NUM_SHIFT (CPU_NUM_SHIFT-THREAD_NUM_BITS)
#define THREAD_NUM_MASK  0x7

#define GNI_INSTID(nic_addr, cpu_num, thr_num) (((nic_addr&NIC_ADDR_MASK)<<NIC_ADDR_SHIFT)|((cpu_num&CPU_NUM_MASK)<<CPU_NUM_SHIFT)|(thr_num&THREAD_NUM_MASK))
#define GNI_NIC_ADDRESS(inst_id)               ((inst_id>>NIC_ADDR_SHIFT)&NIC_ADDR_MASK)
#define GNI_CPU_NUMBER(inst_id)                ((inst_id>>CPU_NUM_SHIFT)&CPU_NUM_MASK)
#define GNI_THREAD_NUMBER(inst_id)             (inst_id&THREAD_NUM_MASK)



enum nnti_gni_pi_ordering_t {
    PI_ORDERING_DEFAULT,
    PI_ORDERING_STRICT,
    PI_ORDERING_RELAXED
};
enum nnti_gni_rdma_mode_t {
    RDMA_FMA,
    RDMA_BTE,
    RDMA_MIXED,
    RDMA_CROSSOVER
};
typedef struct {

    bool use_alps_ptag;

    bool     use_wr_pool;
    uint32_t wr_pool_initial_size;
    uint32_t wr_pool_max_size;
    bool     wr_pool_create_if_empty;

    bool use_rdma_target_ack;
    bool use_rdma_events;
    bool use_rdma_fence;

    enum nnti_gni_pi_ordering_t pi_ordering; /* DEFAULT, STRICT, RELAXED */
    enum nnti_gni_rdma_mode_t   rdma_mode;   /* FMA, BTE, MIXED, CROSSOVER */

    uint32_t fma_bte_crossover_size;

    int32_t max_timeout_ms;

	uint32_t min_atomics_vars;

} nnti_gni_config_t;


/**
 * These states are used to signal events between the completion handler
 * and the main client or server thread.
 *
 * Once CONNECTED, they cycle through RDMA_READ_ADV, RDMA_WRITE_ADV,
 * and RDMA_WRITE_COMPLETE for each ping.
 */
typedef enum {
    IDLE = 1,
    CONNECT_REQUEST,
    ADDR_RESOLVED,
    ROUTE_RESOLVED,
    CONNECTED,
    DISCONNECTED,
    ERROR
} nnti_gni_connection_state_t;

typedef enum {
    SERVER_CONNECTION,
    CLIENT_CONNECTION
} nnti_gni_connection_type_t;

typedef enum {
    REQUEST_BUFFER,
    RECEIVE_BUFFER,
    SEND_BUFFER,
    GET_SRC_BUFFER,
    GET_DST_BUFFER,
    PUT_SRC_BUFFER,
    PUT_DST_BUFFER,
    RDMA_TARGET_BUFFER,
    UNKNOWN_BUFFER
} nnti_gni_buffer_type_t;


#define GNI_OP_PUT_INITIATOR  1
#define GNI_OP_GET_INITIATOR  2
#define GNI_OP_PUT_TARGET     3
#define GNI_OP_GET_TARGET     4
#define GNI_OP_SEND_REQUEST   5
#define GNI_OP_SEND_BUFFER    6
#define GNI_OP_NEW_REQUEST    7
#define GNI_OP_RECEIVE        8
#define GNI_OP_FETCH_ADD      9
#define GNI_OP_COMPARE_SWAP  10


typedef enum {
	NNTI_GNI_SGE_STATE_RESET,
	NNTI_GNI_SGE_STATE_STARTED,
	NNTI_GNI_SGE_STATE_RDMA_COMPLETE,
	NNTI_GNI_SGE_STATE_WC_COMPLETE,
	NNTI_GNI_SGE_STATE_CANCELING,
	NNTI_GNI_SGE_STATE_COMPLETE,
	NNTI_GNI_SGE_STATE_ERROR
} nnti_gni_sge_state_t;
typedef enum {
	NNTI_GNI_WR_STATE_RESET,         // this work request is idle
	NNTI_GNI_WR_STATE_POSTED,        // target is ready to receive data
	NNTI_GNI_WR_STATE_STARTED,       // initiator has posted the FMA/RDMA descriptor
	NNTI_GNI_WR_STATE_ATTACHED,      // target has attached this work request to an NNTI_work_request_t
	NNTI_GNI_WR_STATE_CANCELING,     // NNTI_cancel() called, but not NNTI_wait()
	NNTI_GNI_WR_STATE_RDMA_COMPLETE, // RDMA op complete
//	NNTI_GNI_WR_STATE_WC_COMPLETE,   // work completion transfer complete
	NNTI_GNI_WR_STATE_WAIT_COMPLETE, // NNTI_wait() called and completed successfully
	NNTI_GNI_WR_STATE_ERROR          // something went wrong.  check wr->result
} nnti_gni_wr_state_t;

typedef struct {
    gni_cq_handle_t cq_hdl;
    gni_ep_handle_t ep_hdl;
} nnti_gni_conn_ep_t;


typedef struct {
    uint64_t ack_received;
    uint64_t inst_id;
    uint64_t byte_len;
    uint64_t byte_offset;
    uint64_t src_index;
    uint64_t src_offset;
    uint64_t dest_index;
    uint64_t dest_offset;
    uint16_t op;
} nnti_gni_work_completion_t;

/**
 * attrs to send to the client
 */
typedef struct {
    uint64_t          req_index_addr; /* address of the request index var on the server */
    gni_mem_handle_t  req_index_mem_hdl;

    uint64_t         req_buffer_addr; /* address of the request buffer */
    uint64_t         req_size;        /* the maximum size of a request */
    uint64_t         req_count;       /* the number of requests that will fit in the queue */
    gni_mem_handle_t req_mem_hdl;

} nnti_gni_recv_queue_attrs_t;

typedef struct {
    uint64_t         req_index;      /* after AMO Fetch-Add, this buffer contains the current index of the request queue */
    uint64_t         req_index_addr; /* address of the request index var on the client. */
    gni_cq_handle_t  req_index_mem_cq_hdl;
    gni_mem_handle_t req_index_mem_hdl;
    gni_cq_handle_t  req_index_cq_hdl;
    gni_ep_handle_t  req_index_ep_hdl;

    gni_ep_handle_t  req_ep_hdl;

    uint64_t last_offset;
} nnti_gni_client_queue_t;

typedef struct {
    NNTI_peer_t       peer;
    char             *peer_name;
    NNTI_ip_addr      peer_addr;
    NNTI_tcp_port     peer_port;
    uint32_t          peer_cookie;
    uint32_t          peer_ptag;
    NNTI_instance_id  peer_instance;

    alpsAppGni_t      peer_alps_info;

    nnti_gni_client_queue_t     queue_local_attrs;
    nnti_gni_recv_queue_attrs_t queue_remote_attrs;

    uint64_t         atomics_addr;
    gni_mem_handle_t atomics_mem_hdl;

    gni_cdm_handle_t     cdm_hdl;  /* a client creates this comm domain with params from the server */
    gni_nic_handle_t     nic_hdl;
    gni_ep_handle_t      ep_hdl;

    nnti_gni_connection_state_t state;

    nnti_gni_connection_type_t connection_type;
} nnti_gni_connection_t;


// forward declaration
struct nnti_gni_work_request_t;

typedef struct {
	nnti_gni_work_request_t *gni_wr;
    gni_post_descriptor_t    post_desc;
    nnti_gni_sge_state_t     state;
} nnti_gni_sge_t;  /* scatter-gather element */

typedef struct nnti_gni_work_request_t {
    NNTI_work_request_t       *nnti_wr;

    nnti_gni_wr_state_t state;

    uint8_t                  is_initiator;

    nthread_lock_t           lock;

    const NNTI_buffer_t     *reg_buf;

    nnti_gni_sge_t           sge;
    nnti_gni_sge_t          *sge_list;
    uint32_t                 sge_count;

    nnti_gni_work_completion_t  rdma_wc;  /* RDMA buffers don't have a work completion footer.  use this work for RDMA ops. */
    nnti_gni_work_completion_t *wc;

    uint8_t                  last_op;

    NNTI_instance_id         peer_instance;
} nnti_gni_work_request_t;

typedef std::deque<nnti_gni_work_request_t *>           wr_queue_t;
typedef std::deque<nnti_gni_work_request_t *>::iterator wr_queue_iter_t;

typedef struct {
    nnti_gni_buffer_type_t   type;
    gni_mem_handle_t  mem_hdl;
    gni_mem_handle_t *mem_hdl_list;
    uint32_t          mem_hdl_count;
    wr_queue_t       *wr_queue;
    nthread_lock_t    wr_queue_lock;
    uint32_t          ref_count;
    uint64_t          extra;
} nnti_gni_memory_handle_t;

typedef struct {
    NNTI_buffer_t   *reg_buf;

    uint64_t         last_index_before_reset;
    uint64_t         req_index;       /* index of the next available slot in the request queue */
    uint64_t         req_index_addr;  /* address of the request index var on the server */
    gni_cq_handle_t  req_index_mem_cq_hdl;
    gni_mem_handle_t req_index_mem_hdl;

    char            *req_buffer;      /* pointer to the head of the request buffer */
    uint64_t         req_size;        /* the maximum size of a request */
    uint64_t         req_count;       /* the number of requests that will fit in the queue */
    uint64_t         req_buffer_size; /* the size of the request buffer in bytes (req_size*req_count) */

    int64_t          head;  /* index of the next work request to assign */
    int64_t          tail;  /* index of lowest assigned work request */
} nnti_gni_request_queue_handle_t;

typedef struct {
    uint16_t         delivery_mode;

    gni_cdm_handle_t cdm_hdl;
    gni_nic_handle_t nic_hdl;
    NNTI_instance_id instance;

    gni_cq_handle_t  ep_cq_hdl;
    gni_cq_handle_t  mem_cq_hdl;

    char             interrupt_buf;
    gni_mem_handle_t interrupt_mem_hdl;
    gni_cq_handle_t  interrupt_mem_cq_hdl;
    gni_ep_handle_t  interrupt_ep_hdl;
    gni_cq_handle_t  interrupt_ep_cq_hdl;

    uint64_t         apid;
    alpsAppGni_t     alps_info;

    int              listen_sock;
    char             listen_name[NNTI_HOSTNAME_LEN];
    uint32_t         listen_addr;  /* in NBO */
    uint16_t         listen_port;  /* in NBO */

    int64_t         *atomics;
    gni_mem_handle_t atomics_mem_hdl;
    nthread_lock_t   atomics_lock;

    nnti_gni_request_queue_handle_t req_queue;
} nnti_gni_transport_global_t;




static nthread_lock_t nnti_gni_lock;
static nthread_lock_t nnti_mem_lock;
static nthread_lock_t nnti_progress_lock;
static nthread_cond_t nnti_progress_cond;


static void *aligned_malloc(
        size_t size);
static NNTI_result_t setup_atomics(void);
static gni_mem_handle_t register_memory_segment(
		nnti_gni_buffer_type_t type,
		void *buf,
		uint64_t len,
		uint64_t extra);
static NNTI_result_t register_memory(
        const NNTI_transport_t *trans_hdl,
        char                   *buffer,
        const uint64_t          element_size,
        const uint64_t          extra,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf);
static NNTI_result_t unregister_memory(
        gni_mem_handle_t mem_hdl);
static nnti_gni_sge_t *decode_sge(
        gni_cq_entry_t      *ev_data);
static nnti_gni_work_request_t *decode_work_request(
        gni_cq_entry_t      *ev_data);
static int cancel_wr(
        nnti_gni_work_request_t *gni_wr);
static int process_event(
        const NNTI_buffer_t   *reg_buf,
        nnti_gni_sge_t        *gni_sge,
        gni_cq_handle_t        cq_hdl,
        gni_cq_entry_t        *ev_data,
        gni_post_descriptor_t *post_desc_ptr);
static NNTI_result_t post_recv_queue_work_requests(
        NNTI_buffer_t *reg_buf);
static NNTI_result_t post_recv_work_request(
        NNTI_buffer_t *reg_buf);
static NNTI_result_t repost_recv_work_request(
        NNTI_buffer_t    *reg_buf,
        nnti_gni_work_request_t *wr);
static int8_t is_wr_canceling(
        nnti_gni_work_request_t *gni_wr);
static int8_t is_wr_complete(
        nnti_gni_work_request_t *wr);
static int8_t is_any_wr_complete(
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count,
        uint32_t             *which_wr);
static int8_t is_all_wr_complete(
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count);
static void create_status(
        NNTI_work_request_t  *wr,
        nnti_gni_work_request_t     *gni_wr,
        NNTI_result_t         nnti_rc,
        gni_cq_entry_t       *ev_data,
        NNTI_status_t        *status);
static void create_peer(NNTI_peer_t *peer,
        char *name,
        NNTI_ip_addr addr,
        NNTI_tcp_port port,
        uint32_t ptag,
        uint32_t cookie,
        NNTI_instance_id instance);
static int init_server_listen_socket(void);
static int check_listen_socket_for_new_connections(void);
static uint32_t get_cpunum(void);
static void get_alps_info(alpsAppGni_t *alps_info);
static int tcp_read(int sock, void *incoming, size_t len);
static int tcp_write(int sock, const void *outgoing, size_t len);
static int tcp_exchange(int sock, int is_server, void *incoming, void *outgoing, size_t len);
static void transition_connection_to_ready(
        int sock,
        nnti_gni_connection_t *conn);
static NNTI_result_t get_ipaddr(
        char *ipaddr,
        int maxlen);
static NNTI_result_t init_connection(
        nnti_gni_connection_t **conn,
        const int sock,
        const int is_server);
static void close_connection(nnti_gni_connection_t *c);
static void print_wc(
        const nnti_gni_work_completion_t *wc);
static void print_cq_event(
        const gni_cq_entry_t *event,
        const bool            force);
static void print_post_desc(
        const gni_post_descriptor_t *post_desc_ptr);
static int need_mem_cq(nnti_gni_buffer_type_t type);
static int need_wc_mem_cq(nnti_gni_memory_handle_t *gni_mem_hdl);
static void print_gni_conn(nnti_gni_connection_t *c);
static NNTI_result_t insert_conn_peer(const NNTI_peer_t *peer, nnti_gni_connection_t *conn);
static NNTI_result_t insert_conn_instance(const NNTI_instance_id instance, nnti_gni_connection_t *conn);
static nnti_gni_connection_t *get_conn_peer(const NNTI_peer_t *peer);
static nnti_gni_connection_t *get_conn_instance(const NNTI_instance_id instance);
static NNTI_peer_t *get_peer_by_url(const char *url);
static nnti_gni_connection_t *del_conn_peer(const NNTI_peer_t *peer);
static nnti_gni_connection_t *del_conn_instance(const NNTI_instance_id instance);
static void print_peer_map(void);
static void print_instance_map(void);
static void close_all_conn(void);
//static void print_put_buf(void *buf, uint32_t size);
static void print_raw_buf(void *buf, uint32_t size);

static NNTI_result_t insert_buf_bufhash(NNTI_buffer_t *buf);
static NNTI_buffer_t *get_buf_bufhash(const uint32_t bufhash);
static NNTI_buffer_t *del_buf_bufhash(NNTI_buffer_t *victim);
static void print_bufhash_map(void);

static NNTI_result_t insert_sge_sgehash(nnti_gni_sge_t *);
static nnti_gni_sge_t *get_sge_sgehash(const uint32_t bufhash);
static nnti_gni_sge_t *del_sge_sgehash(nnti_gni_sge_t *victim);
static void print_sgehash_map(void);

static NNTI_result_t insert_wr_wrhash(nnti_gni_work_request_t *);
static nnti_gni_work_request_t *get_wr_wrhash(const uint32_t bufhash);
static nnti_gni_work_request_t *del_wr_wrhash(nnti_gni_work_request_t *victim);
static void print_wrhash_map(void);

static NNTI_result_t wr_pool_init(void);
static nnti_gni_work_request_t *wr_pool_target_pop(void);
static nnti_gni_work_request_t *wr_pool_initiator_pop(void);
static void wr_pool_target_push(nnti_gni_work_request_t *wr);
static void wr_pool_initiator_push(nnti_gni_work_request_t *wr);
static NNTI_result_t wr_pool_fini(void);

static uint16_t get_dlvr_mode_from_env();
static void set_dlvr_mode(
        gni_post_descriptor_t *pd);
static void set_rdma_mode(
        gni_post_descriptor_t *pd);
static void set_post_desc(
        gni_post_descriptor_t *pd,
        nnti_gni_buffer_type_t target_buf_type,
        uint32_t buf_length);
static void set_wc_post_desc(
        gni_post_descriptor_t *pd,
        uint32_t buf_length);

static void server_req_queue_init(
        nnti_gni_request_queue_handle_t *q,
        char                     *buffer,
        uint64_t                  req_size,
        uint64_t                  extra,
        uint64_t                  req_count);
static void server_req_queue_destroy(
        nnti_gni_request_queue_handle_t *q);

static void client_req_queue_init(
        nnti_gni_connection_t *c);
static void client_req_queue_destroy(
        nnti_gni_connection_t *c);

static int reset_req_index(
        nnti_gni_request_queue_handle_t  *req_queue_attrs);

static void send_rdma_wc (
        nnti_gni_work_request_t    *wr,
        const NNTI_buffer_t *local_buf,
        const NNTI_buffer_t *remote_buf);

static int fetch_add_buffer_offset(
        nnti_gni_client_queue_t     *local_req_queue_attrs,
        nnti_gni_recv_queue_attrs_t *remote_req_queue_attrs,
        uint64_t                   addend,
        uint64_t                  *prev_offset);
static int send_req(
        const NNTI_peer_t         *peer_hdl,
        nnti_gni_client_queue_t     *local_req_queue_attrs,
        nnti_gni_recv_queue_attrs_t *remote_req_queue_attrs,
        uint64_t                   offset,
        const NNTI_buffer_t       *reg_buf,
        nnti_gni_work_request_t          *gni_wr);
static int request_send(
        const NNTI_peer_t         *peer_hdl,
        nnti_gni_client_queue_t     *client_q,
        nnti_gni_recv_queue_attrs_t *server_q,
        const NNTI_buffer_t       *reg_buf,
        int                        req_num,
        nnti_gni_work_request_t          *gni_wr);
static int send_buffer(
        const NNTI_peer_t   *peer_hdl,
        const NNTI_buffer_t *src_hdl,
        const NNTI_buffer_t *dest_hdl,
        nnti_gni_work_request_t    *gni_wr);
static int buffer_send(
        const NNTI_peer_t   *peer_hdl,
        const NNTI_buffer_t *src_hdl,
        const NNTI_buffer_t *dest_hdl,
        nnti_gni_work_request_t    *gni_wr);

static NNTI_result_t progress(int timeout);

static void config_init(
        nnti_gni_config_t *c);
static void config_get_from_env(
        nnti_gni_config_t *c);

inline int DEQUEUE_POST_DESCRIPTOR(
        gni_cq_handle_t           cq_hdl,
        gni_cq_entry_t           *ev_data,
        gni_post_descriptor_t   **post_desc_ptr);
inline int IS_RDMA_INITIATOR_BUFFER(
        const nnti_gni_work_request_t *wr);
inline int IS_RDMA_TARGET_BUFFER(
        const nnti_gni_work_request_t *wr);


#define GNI_MEM_HDL(_b)           ((nnti_gni_memory_handle_t *)((_b)->transport_private))
#define GNI_WORK_REQUEST(_gwr)    ((nnti_gni_work_request_t *)((_gwr)->transport_private))

#define GNI_ELEMENT_ADDRESS(_b,_i)  (_b->payload+((_b->payload_size+sizeof(nnti_gni_work_completion_t))*_i))
#define GNI_WC_ADDRESS(_b,_i)       ((nnti_gni_work_completion_t *)(_b->payload + \
		                            ((_b->payload_size+sizeof(nnti_gni_work_completion_t))*_i) + \
		                            _b->payload_size))

#define GNI_ATTACH_WR(_nwr,_gwr)  {             \
        _nwr->transport_private=(uint64_t)_gwr; \
        _gwr->nnti_wr          =_nwr;           \
        }
#define GNI_DETACH_WR(_nwr,_gwr)  {          \
	    log_debug(nnti_debug_level, "_nwr=%p ; _gwr=%p", _nwr, _gwr); \
        _gwr->nnti_wr          =NULL;        \
        _nwr->transport_private=0;           \
        }

static bool gni_initialized=false;

static nnti_gni_transport_global_t transport_global_data;
static const int MIN_TIMEOUT = 100;  /* in milliseconds.  must be >0 for GNI_CqVectorWaitEvent(). */

static log_level nnti_cq_debug_level;
static log_level nnti_event_debug_level;
static log_level nnti_ee_debug_level;



/**
 * This custom key is used to look up existing connections.
 */
struct addrport_key {
    NNTI_ip_addr    addr;       /* part1 of a compound key */
    NNTI_tcp_port   port;       /* part2 of a compound key */

    // Need this operators for the hash map
    bool operator<(const addrport_key &key1) const {
        if (addr < key1.addr) {
            return true;
        } else if (addr == key1.addr) {
            if (port < key1.port) {
                return true;
            }
        }
        return false;
    }
    bool operator>(const addrport_key &key1) const {
        if (addr > key1.addr) {
            return true;
        } else if (addr == key1.addr) {
            if (port > key1.port) {
                return true;
            }
        }
        return false;
    }
};


/* Thomas Wang's 64 bit to 32 bit Hash Function (http://www.concentric.net/~ttwang/tech/inthash.htm) */
static uint32_t hash6432shift(uint64_t key)
{
  key = (~key) + (key << 18); // key = (key << 18) - key - 1;
  key = key ^ (key >> 31);
  key = key * 21;             // key = (key + (key << 2)) + (key << 4);
  key = key ^ (key >> 11);
  key = key + (key << 6);
  key = key ^ (key >> 22);
  return (uint32_t)key;
}

/*
 * We need a couple of maps to keep track of connections.  Servers need to find
 * connections by QP number when requests arrive.  Clients need to find connections
 * by peer address and port.  Setup those maps here.
 */
static std::map<addrport_key, nnti_gni_connection_t *> connections_by_peer;
typedef std::map<addrport_key, nnti_gni_connection_t *>::iterator conn_by_peer_iter_t;
typedef std::pair<addrport_key, nnti_gni_connection_t *> conn_by_peer_t;
static nthread_lock_t nnti_conn_peer_lock;

static std::map<NNTI_instance_id, nnti_gni_connection_t *> connections_by_instance;
typedef std::map<NNTI_instance_id, nnti_gni_connection_t *>::iterator conn_by_inst_iter_t;
typedef std::pair<NNTI_instance_id, nnti_gni_connection_t *> conn_by_inst_t;
static nthread_lock_t nnti_conn_instance_lock;

static std::map<uint32_t, NNTI_buffer_t *> buffers_by_bufhash;
typedef std::map<uint32_t, NNTI_buffer_t *>::iterator buf_by_bufhash_iter_t;
typedef std::pair<uint32_t, NNTI_buffer_t *> buf_by_bufhash_t;
static nthread_lock_t nnti_buf_bufhash_lock;

static std::map<uint32_t, nnti_gni_sge_t *> sge_by_sgehash;
typedef std::map<uint32_t, nnti_gni_sge_t *>::iterator sge_by_sgehash_iter_t;
typedef std::pair<uint32_t, nnti_gni_sge_t *> sge_by_sgehash_t;
static nthread_lock_t nnti_sge_sgehash_lock;

static std::map<uint32_t, nnti_gni_work_request_t *> wr_by_wrhash;
typedef std::map<uint32_t, nnti_gni_work_request_t *>::iterator wr_by_wrhash_iter_t;
typedef std::pair<uint32_t, nnti_gni_work_request_t *> wr_by_wrhash_t;
static nthread_lock_t nnti_wr_wrhash_lock;

typedef std::deque<nnti_gni_work_request_t *>           wr_pool_t;
typedef std::deque<nnti_gni_work_request_t *>::iterator wr_pool_iter_t;
static nthread_lock_t nnti_wr_pool_lock;

static wr_pool_t target_wr_pool;
static wr_pool_t initiator_wr_pool;


static nnti_gni_config_t config;


/* ---------------- Wrappers to protect HPCToolkit issues ---------- */
static
gni_return_t GNI_CqWaitEvent_wrapper(
        gni_cq_handle_t cq_hdl,
        uint64_t        timeout,
        gni_cq_entry_t *ev_data)
{
    gni_return_t  rc;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    rc = GNI_CqWaitEvent(cq_hdl, timeout, ev_data);
    if (sampling) SAMPLING_START();

    return rc;
}

static
gni_return_t GNI_CqVectorWaitEvent_wrapper(
        gni_cq_handle_t *cq_list,
        uint32_t         num_cqs,
        uint64_t         timeout,
        gni_cq_entry_t  *event_data,
        uint32_t        *which)
{
    gni_return_t  rc;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    rc = GNI_CqVectorWaitEvent(cq_list, num_cqs, timeout, event_data, which);
    if (sampling) SAMPLING_START();

    return rc;
}

static
gni_return_t GNI_CdmAttach_wrapper(
        gni_cdm_handle_t cdm_hndl,
        uint32_t device_id,
        uint32_t *local_addr,
        gni_nic_handle_t *nic_hndl)
{
    gni_return_t rc;
    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    rc=GNI_CdmAttach (cdm_hndl, device_id, local_addr, nic_hndl);
    if (sampling) SAMPLING_START();
    return rc;
}

static
gni_return_t GNI_GetCompleted_wrapper(
        gni_cq_handle_t cq_hndl,
        gni_cq_entry_t event_data,
        gni_post_descriptor_t **post_descr)
{
    gni_return_t rc;
    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    rc=GNI_GetCompleted (cq_hndl, event_data, post_descr);
    if (sampling) SAMPLING_START();
    return rc;
}

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
NNTI_result_t NNTI_gni_init (
        const NNTI_transport_id_t  trans_id,
        const char                *my_url,
        NNTI_transport_t          *trans_hdl)
{
    int rc=NNTI_OK;

    trios_declare_timer(call_time);

    char transport[NNTI_URL_LEN];
    char address[NNTI_URL_LEN];
    char *sep;

    char hostname[NNTI_HOSTNAME_LEN];

    uint32_t nic_addr  =0;
    uint32_t cpu_num   =0;
    uint32_t thread_num=0;
    uint32_t gni_cpu_id=0;


    assert(trans_hdl);

    log_debug(nnti_ee_debug_level, "enter");

    log_debug(nnti_debug_level, "my_url=%s", my_url);
    log_debug(nnti_debug_level, "initialized=%d, FALSE==%d", (int)gni_initialized, (int)FALSE);

    if (!gni_initialized) {

        memset(&transport_global_data, 0, sizeof(nnti_gni_transport_global_t));

        nnti_cq_debug_level=nnti_debug_level;
        nnti_event_debug_level=nnti_debug_level;
        nnti_ee_debug_level=nnti_debug_level;

        nthread_lock_init(&nnti_gni_lock);
        nthread_lock_init(&nnti_mem_lock);

        nthread_lock_init(&nnti_progress_lock);
        nthread_cond_init(&nnti_progress_cond);

        // initialize the mutexes for the connection maps
        nthread_lock_init(&nnti_conn_peer_lock);
        nthread_lock_init(&nnti_conn_instance_lock);
        nthread_lock_init(&nnti_wr_wrhash_lock);
        nthread_lock_init(&nnti_buf_bufhash_lock);

        nthread_lock_init(&nnti_wr_pool_lock);

        nthread_lock_init(&transport_global_data.atomics_lock);

        config_init(&config);
        config_get_from_env(&config);

        log_debug(nnti_debug_level, "my_url=%s", my_url);

        hostname[0]='\0';
        if (my_url != NULL) {
            if ((rc=nnti_url_get_transport(my_url, transport, NNTI_URL_LEN)) != NNTI_OK) {
                goto cleanup;
            }
            if (0!=strcmp(transport, "gni")) {
                rc=NNTI_EINVAL;
                goto cleanup;
            }

            if ((rc=nnti_url_get_address(my_url, address, NNTI_URL_LEN)) != NNTI_OK) {
                goto cleanup;
            }

            sep=strchr(address, ':');
            if (sep == address) {
                /* no hostname given; try gethostname */
                gethostname(hostname, NNTI_HOSTNAME_LEN);
            } else {
                strncpy(hostname, address, sep-address);
                hostname[sep-address]='\0';
            }
            sep++;
        } else {
            rc=get_ipaddr(hostname, NNTI_HOSTNAME_LEN);
            if (rc != NNTI_OK) {
                log_error(nnti_debug_level, "could not find IP address to listen on");
                goto cleanup;
            }
        }
        strcpy(transport_global_data.listen_name, hostname);


        transport_global_data.delivery_mode = get_dlvr_mode_from_env();


        log_debug(nnti_debug_level, "initializing Gemini");

//        /* register trace groups (let someone else enable) */
//        trace_counter_gid = trace_get_gid(TRACE_RPC_COUNTER_GNAME);
//        trace_interval_gid = trace_get_gid(TRACE_RPC_INTERVAL_GNAME);


        trios_start_timer(call_time);
        get_alps_info(&transport_global_data.alps_info);
        trios_stop_timer("get_alps_info", call_time);

        trios_start_timer(call_time);
        rc=GNI_CdmGetNicAddress (transport_global_data.alps_info.device_id, &nic_addr, &gni_cpu_id);
        trios_stop_timer("CdmGetNicAddress", call_time);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "CdmGetNicAddress() failed: %d", rc);
            if (rc==GNI_RC_NO_MATCH)
                rc=NNTI_EEXIST;
            else
                rc=NNTI_EINVAL;

            goto cleanup;
        }

        trios_start_timer(call_time);
        cpu_num = get_cpunum();
        trios_stop_timer("get_cpunum", call_time);

        transport_global_data.instance=GNI_INSTID(nic_addr, cpu_num, thread_num);
        log_debug(nnti_debug_level, "nic_addr(%llu), cpu_num(%llu), thread_num(%llu), inst_id(%llu), "
                "derived.nic_addr(%llu), derived.cpu_num(%llu), derived.thread_num(%llu)",
                (uint64_t)nic_addr, (uint64_t)cpu_num, (uint64_t)thread_num,
                (uint64_t)transport_global_data.instance,
                (uint64_t)GNI_NIC_ADDRESS(transport_global_data.instance),
                (uint64_t)GNI_CPU_NUMBER(transport_global_data.instance),
                (uint64_t)GNI_THREAD_NUMBER(transport_global_data.instance));

        log_debug(nnti_debug_level, "global_nic_hdl - host(%s) device_id(%llu) local_addr(%lld) cookie(%llu) ptag(%llu) "
                    "apid(%llu) inst_id(%llu) gni_nic_addr(%llu) gni_cpu_id(%llu) linux_cpu_num(%llu) omp_thread_num(%llu)",
                    transport_global_data.listen_name,
                    (unsigned long long)transport_global_data.alps_info.device_id,
                    (long long)transport_global_data.alps_info.local_addr,
                    (unsigned long long)transport_global_data.alps_info.cookie,
                    (unsigned long long)transport_global_data.alps_info.ptag,
                    (unsigned long long)transport_global_data.apid,
                    (unsigned long long)transport_global_data.instance,
                    (unsigned long long)nic_addr,
                    (unsigned long long)gni_cpu_id,
                    (unsigned long long)cpu_num,
                    (unsigned long long)thread_num);

        trios_start_timer(call_time);
        rc=GNI_CdmCreate(transport_global_data.instance,
                transport_global_data.alps_info.ptag,
                transport_global_data.alps_info.cookie,
                GNI_CDM_MODE_ERR_NO_KILL,
//                GNI_CDM_MODE_ERR_NO_KILL|GNI_CDM_MODE_BTE_SINGLE_CHANNEL,
                &transport_global_data.cdm_hdl);
        trios_stop_timer("CdmCreate", call_time);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "CdmCreate() failed: %d", rc);
            rc=NNTI_EINVAL;
            goto cleanup;
        }

        trios_start_timer(call_time);

        rc=GNI_CdmAttach_wrapper (transport_global_data.cdm_hdl,
                transport_global_data.alps_info.device_id,
                (uint32_t*)&transport_global_data.alps_info.local_addr, /* ALPS and GNI disagree about the type of local_addr.  cast here. */
                &transport_global_data.nic_hdl);

        trios_stop_timer("CdmAttach", call_time);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "CdmAttach() failed: %d", rc);
            if (rc==GNI_RC_PERMISSION_ERROR)
                rc=NNTI_EPERM;
            else
                rc=NNTI_EINVAL;

            goto cleanup;
        }

        rc=GNI_CqCreate (transport_global_data.nic_hdl, 10000, 0, GNI_CQ_BLOCKING, NULL, NULL, &transport_global_data.ep_cq_hdl);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "CqCreate(transport_global_data.ep_cq_hdl) failed: %d", rc);
            goto cleanup;
        }
        rc=GNI_CqCreate (transport_global_data.nic_hdl, 10000, 0, GNI_CQ_BLOCKING, NULL, NULL, &transport_global_data.mem_cq_hdl);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "CqCreate(transport_global_data.mem_cq_hdl) failed: %d", rc);
            goto cleanup;
        }


        rc=GNI_CqCreate (transport_global_data.nic_hdl, 10, 0, GNI_CQ_BLOCKING, NULL, NULL, &transport_global_data.interrupt_mem_cq_hdl);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "CqCreate(transport_global_data.interrupt_cq_hdl) failed: %d", rc);
            goto cleanup;
        }
        rc=GNI_MemRegister (transport_global_data.nic_hdl,
        				    (uint64_t)&transport_global_data.interrupt_buf,
        				    sizeof(uint32_t),
        				    transport_global_data.interrupt_mem_cq_hdl,
        				    GNI_MEM_READWRITE,
        				    (uint32_t)-1,
        				    &transport_global_data.interrupt_mem_hdl);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "MemRegister(transport_global_data.interrupt_mem_hdl) failed: rc=%d, %s", rc, strerror(errno));
            goto cleanup;
        }

        rc=GNI_CqCreate (transport_global_data.nic_hdl, 1, 0, GNI_CQ_BLOCKING, NULL, NULL, &transport_global_data.interrupt_ep_cq_hdl);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "CqCreate(transport_global_data.interrupt_cq_hdl) failed: %d", rc);
            goto cleanup;
        }
        rc=GNI_EpCreate (transport_global_data.nic_hdl, transport_global_data.interrupt_ep_cq_hdl, &transport_global_data.interrupt_ep_hdl);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "EpCreate(transport_global_data.interrupt_ep_hdl) failed: %d", rc);
            goto cleanup;
        }
        rc=GNI_EpBind (transport_global_data.interrupt_ep_hdl, transport_global_data.alps_info.local_addr, transport_global_data.instance);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "EpBind(transport_global_data.interrupt_ep_hdl) failed: %d", rc);
            goto cleanup;
        }


        if (config.use_wr_pool) {
            rc=wr_pool_init();
            if (rc!=NNTI_OK) {
                log_error(nnti_debug_level, "wr_pool_init(): %d", rc);
                goto cleanup;
            }
        }

        setup_atomics();

        trios_start_timer(call_time);
        init_server_listen_socket();
        trios_stop_timer("init_server_listen_socket", call_time);

        if (logging_info(nnti_debug_level)) {
            fprintf(logger_get_file(), "Gemini Initialized: host(%s) port(%u)\n",
                    transport_global_data.listen_name,
                    ntohs(transport_global_data.listen_port));
        }

        trios_start_timer(call_time);
        create_peer(
                &trans_hdl->me,
                transport_global_data.listen_name,
                transport_global_data.listen_addr,
                transport_global_data.listen_port,
                transport_global_data.alps_info.ptag,
                transport_global_data.alps_info.cookie,
                transport_global_data.instance);
        trios_stop_timer("create_peer", call_time);

        gni_initialized = true;
    }

cleanup:
    log_debug(nnti_ee_debug_level, "exit");

    return((NNTI_result_t)rc);
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
 *                - ex. "ptl://nid:pid/", "gni://ip_addr:port"
 *    - memory_descriptor - (optional) transport-specific representation of RMA params
 */
NNTI_result_t NNTI_gni_get_url (
        const NNTI_transport_t *trans_hdl,
        char                   *url,
        const uint64_t          maxlen)
{
    NNTI_result_t rc=NNTI_OK;

    assert(trans_hdl);
    assert(url);
    assert(maxlen>0);

    log_debug(nnti_ee_debug_level, "enter");

    strncpy(url, trans_hdl->me.url, maxlen);
    url[maxlen-1]='\0';

    log_debug(nnti_ee_debug_level, "exit");

    return(rc);
}


/**
 * @brief Prepare for communication with the peer identified by <tt>url</tt>.
 *
 * Parse <tt>url</tt> in a transport specific way.  Perform any transport specific
 * actions necessary to begin communication with this peer.
 *
 *
 * Connectionless transport: parse and populate
 * Connected transport: parse, connection and populate
 *
 */
NNTI_result_t NNTI_gni_connect (
        const NNTI_transport_t *trans_hdl,
        const char             *url,
        const int               timeout,
        NNTI_peer_t            *peer_hdl)
{
    int rc=NNTI_OK;

    trios_declare_timer(call_time);

    char transport[NNTI_URL_LEN];
    char address[NNTI_URL_LEN];
    char params[NNTI_URL_LEN];
    char *sep;

    NNTI_peer_t *existing_peer=NULL;

    char          hostname[NNTI_HOSTNAME_LEN];
    char          port_str[NNTI_HOSTNAME_LEN];
    NNTI_tcp_port port;

    char     *cookie_str;
    char     *ptag_str;

    int s;
    struct hostent *host_entry;
    struct sockaddr_in skin;
    socklen_t skin_size=sizeof(struct sockaddr_in);

    nnti_gni_connection_t *conn=NULL;

    long start_time;
    long elapsed_time = 0;
    long timeout_per_call;


    assert(trans_hdl);
    assert(peer_hdl);

    log_debug(nnti_ee_debug_level, "enter (url=%s)", url);

    existing_peer=get_peer_by_url(url);
    if (existing_peer!=NULL) {
        *peer_hdl=*existing_peer;
        return(NNTI_OK);
    }

    conn = (nnti_gni_connection_t *)calloc(1, sizeof(nnti_gni_connection_t));
    log_debug(nnti_debug_level, "calloc returned conn=%p.", conn);
    if (conn == NULL) {
        log_error(nnti_debug_level, "calloc returned NULL.  out of memory?: %s", strerror(errno));
        rc=NNTI_ENOMEM;
        goto cleanup;
    }

    if (url != NULL) {
        if ((rc=nnti_url_get_transport(url, transport, NNTI_URL_LEN)) != NNTI_OK) {
            goto cleanup;
        }
        if (0!=strcmp(transport, "gni")) {
            /* the peer described by 'url' is not an Gemini peer */
            rc=NNTI_EINVAL;
            goto cleanup;
        }

        if ((rc=nnti_url_get_address(url, address, NNTI_URL_LEN)) != NNTI_OK) {
            /* the peer described by 'url' is not an Gemini peer */
            rc=NNTI_EINVAL;
            goto cleanup;
        }

        if ((rc=nnti_url_get_params(url, params, NNTI_URL_LEN)) != NNTI_OK) {
            /* the peer described by 'url' is not an Gemini peer */
            rc=NNTI_EINVAL;
            goto cleanup;
        }

        sep=strchr(address, ':');
        strncpy(hostname, address, sep-address);
        hostname[sep-address]='\0';
        strcpy(port_str, sep+1);
        port=strtol(port_str, NULL, 0);

        log_debug(nnti_ee_debug_level, "params=%s", params);

        ptag_str=strstr(params, "ptag=");
        sep=strchr(ptag_str, '&');
        *sep='\0';
        log_debug(nnti_ee_debug_level, "ptag_str=%s", ptag_str+5);
        conn->peer_ptag=strtol(ptag_str+5, NULL, 10);
        *sep='&';

        cookie_str=strstr(params, "cookie=");
        log_debug(nnti_ee_debug_level, "cookie_str=%s", cookie_str+7);
        conn->peer_cookie=strtoull(cookie_str+7, NULL, 10);

        log_debug(nnti_ee_debug_level, "url=%s", url);

    } else {
        /*  */
        rc=NNTI_EINVAL;
        goto cleanup;
    }

    if (config.use_alps_ptag) {
        conn->cdm_hdl = transport_global_data.cdm_hdl;
        conn->nic_hdl = transport_global_data.nic_hdl;
    } else {
        log_debug(nnti_ee_debug_level, "conn_nic_hdl - host(%s) device_id(%llu) local_addr(%lld) cookie(%llu) ptag(%llu) "
                "apid(%llu) inst_id(%llu) gni_nic_addr(%llu) linux_cpu_num(%llu) omp_thread_num(%llu)",
                transport_global_data.listen_name,
                (unsigned long long)transport_global_data.alps_info.device_id,
                (long long)transport_global_data.alps_info.local_addr,
                (unsigned long long)conn->peer_cookie,
                (unsigned long long)conn->peer_ptag,
                (unsigned long long)transport_global_data.apid,
                (unsigned long long)transport_global_data.instance,
                (uint64_t)GNI_NIC_ADDRESS(transport_global_data.instance),
                (uint64_t)GNI_CPU_NUMBER(transport_global_data.instance),
                (uint64_t)GNI_THREAD_NUMBER(transport_global_data.instance));

        trios_start_timer(call_time);
        rc=GNI_CdmCreate(transport_global_data.instance,
                conn->peer_ptag,
                conn->peer_cookie,
                GNI_CDM_MODE_ERR_ALL_KILL,
                &conn->cdm_hdl);
        trios_stop_timer("CdmCreate", call_time);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "CdmCreate() failed: %d", rc);
            rc=NNTI_EINVAL;
            goto cleanup;
        }

        trios_start_timer(call_time);
        rc=GNI_CdmAttach_wrapper (conn->cdm_hdl,
                transport_global_data.alps_info.device_id,
                (uint32_t*)&transport_global_data.alps_info.local_addr, /* ALPS and GNI disagree about the type of local_addr.  cast here. */
                &conn->nic_hdl);
        trios_stop_timer("CdmAttach", call_time);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "CdmAttach() failed: %d", rc);
            if (rc==GNI_RC_PERMISSION_ERROR)
                rc=NNTI_EPERM;
            else
                rc=NNTI_EINVAL;
            goto cleanup;
        }

        rc=GNI_CqCreate (conn->nic_hdl, 10000, 0, GNI_CQ_BLOCKING, NULL, NULL, &transport_global_data.ep_cq_hdl);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "CqCreate(transport_global_data.ep_cq_hdl) failed: %d", rc);
            goto cleanup;
        }
        rc=GNI_CqCreate (conn->nic_hdl, 10000, 0, GNI_CQ_BLOCKING, NULL, NULL, &transport_global_data.mem_cq_hdl);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "CqCreate(transport_global_data.mem_cq_hdl) failed: %d", rc);
            goto cleanup;
        }
    }

    host_entry = gethostbyname(hostname);
    if (!host_entry) {
        log_warn(nnti_debug_level, "failed to resolve server name (%s): %s", hostname, strerror(errno));
        rc=NNTI_ENOENT;
        goto cleanup;
    }
    memset(&skin, 0, sizeof(skin));
    skin.sin_family = host_entry->h_addrtype;
    memcpy(&skin.sin_addr, host_entry->h_addr_list[0], (size_t) host_entry->h_length);
    skin.sin_port = htons(port);

    elapsed_time=0;
    timeout_per_call = MIN_TIMEOUT;

    s = socket(AF_INET, SOCK_STREAM, 0);
    if (s < 0) {
        log_warn(nnti_debug_level, "failed to create tcp socket: errno=%d (%s)", errno, strerror(errno));
        rc=NNTI_EIO;
        goto cleanup;
    }
    trios_start_timer(call_time);
    start_time=trios_get_time_ms();
    while((timeout==-1) || (elapsed_time < timeout)) {
        log_debug(nnti_debug_level, "calling connect");
        if (connect(s, (struct sockaddr *)&skin, skin_size) == 0) {
            log_debug(nnti_debug_level, "connected");
            break;
        }
        elapsed_time=trios_get_time_ms()-start_time;
        log_warn(nnti_debug_level, "failed to connect to server (%s:%u): errno=%d (%s)", hostname, port, errno, strerror(errno));
        if ((timeout>0) && (elapsed_time >= timeout)) {
            rc=NNTI_EIO;
            goto cleanup;
        }
        nnti_sleep(timeout_per_call);
    }
    trios_stop_timer("socket connect", call_time);

    trios_start_timer(call_time);
    rc = init_connection(&conn, s, 0);
    trios_stop_timer("gni init connection", call_time);
    if (conn==NULL) {
        rc=NNTI_EIO;
        goto cleanup;
    }

    create_peer(
            peer_hdl,
            conn->peer_name,
            conn->peer_addr,
            conn->peer_port,
            conn->peer_ptag,
            conn->peer_cookie,
            conn->peer_instance);

    conn->peer=*peer_hdl;

    insert_conn_peer(peer_hdl, conn);
    insert_conn_instance(conn->peer_instance, conn);

    transition_connection_to_ready(s, conn);

    if (close(s) < 0) {
        log_warn(nnti_debug_level, "failed to close tcp socket: errno=%d (%s)", errno, strerror(errno));
        rc=NNTI_EIO;
        goto cleanup;
    }

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_peer(logger_get_file(), "peer_hdl",
                "end of NNTI_gni_connect", peer_hdl);
    }

cleanup:
    if (rc != NNTI_OK) {
        if (conn!=NULL) free(conn);
    }
    log_debug(nnti_ee_debug_level, "exit");

    return((NNTI_result_t)rc);
}


/**
 * @brief Terminate communication with this peer.
 *
 * Perform any transport specific actions necessary to end communication with
 * this peer.
 */
NNTI_result_t NNTI_gni_disconnect (
        const NNTI_transport_t *trans_hdl,
        NNTI_peer_t            *peer_hdl)
{
    NNTI_result_t rc=NNTI_OK;

    assert(trans_hdl);
    assert(peer_hdl);

    log_debug(nnti_ee_debug_level, "enter");

    nnti_gni_connection_t *conn=get_conn_peer(peer_hdl);

    close_connection(conn);

    del_conn_peer(peer_hdl);
    del_conn_instance(conn->peer_instance);

    log_debug(nnti_ee_debug_level, "exit");

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
NNTI_result_t NNTI_gni_alloc (
        const NNTI_transport_t *trans_hdl,
        const uint64_t          element_size,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf)
{
    NNTI_result_t nnti_rc=NNTI_OK;

    log_debug(nnti_debug_level, "enter");

    /* each element needs some extra space to store the work completion data */
    uint32_t transport_extra=sizeof(nnti_gni_work_completion_t);

    assert(trans_hdl);
    assert(element_size>0);
    assert(num_elements>0);
    assert(ops>0);
    assert(reg_buf);

    char *buf=(char *)aligned_malloc((element_size+transport_extra)*num_elements);
    assert(buf);

    nnti_rc=register_memory(
            trans_hdl,
            buf,
            element_size,
            transport_extra,
            num_elements,
            ops,
            reg_buf);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "reg_buf",
                "end of NNTI_gni_alloc", reg_buf);
    }

    log_debug(nnti_debug_level, "exit");

    return(nnti_rc);
}


/**
 * @brief Cleanup after network operations are complete.
 *
 * Destroy an NNTI_buffer_t that was previously created by NNTI_regsiter_buffer().
 * It is the user's responsibility to release the the memory region.
 */
NNTI_result_t NNTI_gni_free (
        NNTI_buffer_t *reg_buf)
{
    NNTI_result_t nnti_rc=NNTI_OK;

    log_debug(nnti_debug_level, "enter");

    assert(reg_buf);

    char *buf=NNTI_BUFFER_C_POINTER(reg_buf);
    assert(buf);

    nnti_rc=NNTI_gni_unregister_memory(reg_buf);

    free(buf);

    log_debug(nnti_debug_level, "exit");

    return(nnti_rc);
}


/**
 * @brief Prepare a block of memory for network operations.
 *
 * Wrap a user allocated block of memory in an NNTI_buffer_t.  The transport
 * may take additional actions to prepare the memory for network send/receive.
 * If the memory block doesn't meet the transport's requirements for memory
 * regions, then errors or poor performance may result.
 */
NNTI_result_t NNTI_gni_register_memory (
        const NNTI_transport_t *trans_hdl,
        char                   *buffer,
        const uint64_t          element_size,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf)
{
    NNTI_result_t rc=NNTI_OK;

    rc=register_memory(
    		trans_hdl,
            buffer,
            element_size,
            0, //extra
            num_elements,
            ops,
            reg_buf);

    return(rc);
}


/**
 * @brief Prepare a list of memory segments for network operations.
 *
 * Wrap a list of user allocated memory segments in an NNTI_buffer_t.  The
 * transport may take additional actions to prepare the memory segments for
 * network send/receive.  If the memory segments don't meet the transport's
 * requirements for memory regions, then errors or poor performance may
 * result.
 *
 */
NNTI_result_t NNTI_gni_register_segments (
        const NNTI_transport_t *trans_hdl,
        char                  **segments,
        const uint64_t         *segment_lengths,
        const uint64_t          num_segments,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf)
{
	NNTI_result_t nnti_rc=NNTI_OK;
	gni_return_t  gni_rc =GNI_RC_SUCCESS;

	uint32_t flags=GNI_MEM_READWRITE;

	trios_declare_timer(call_time);

	gni_cq_handle_t mem_cq_hdl=NULL;

    NNTI_buffer_t     *old_buf=NULL;
    nnti_gni_memory_handle_t *gni_mem_hdl=NULL;

    assert(trans_hdl);
    assert(segments);
    assert(segment_lengths);
    assert(num_segments>0);
    assert(ops>0);
    assert(reg_buf);

    log_debug(nnti_ee_debug_level, "enter (reg_buf=%p)", reg_buf);

    if ((ops == NNTI_SEND_SRC) || (ops == NNTI_RECV_DST) || (ops == NNTI_RECV_QUEUE)) {
        log_debug(nnti_debug_level, "NNTI_SEND_SRC, NNTI_RECV_DST and NNTI_RECV_QUEUE types cannot be segmented.");
        return(NNTI_EINVAL);
    }

nthread_lock(&nnti_mem_lock);

	old_buf=get_buf_bufhash(hash6432shift((uint64_t)segments[0]));
    if (old_buf==NULL) {
        gni_mem_hdl=(nnti_gni_memory_handle_t*)calloc(1, sizeof(nnti_gni_memory_handle_t));
        assert(gni_mem_hdl);
        gni_mem_hdl->wr_queue =new wr_queue_t;
        gni_mem_hdl->ref_count=1;
        nthread_lock_init(&gni_mem_hdl->wr_queue_lock);
    } else {


		// check that the number of old_buf segments equals num_segments.
		if (old_buf->buffer_segments.NNTI_remote_addr_array_t_len != num_segments) {
			log_debug(nnti_debug_level, "Segment count mismatch (old=%d new=%d).  Aborting...",
					old_buf->buffer_segments.NNTI_remote_addr_array_t_len, num_segments);
			return(NNTI_EINVAL);
		}
		// check that all the old_buf segments match the current list of segments
		for (int i=0;i<old_buf->buffer_segments.NNTI_remote_addr_array_t_len;i++) {
			if (old_buf->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.buf != (uint64_t)segments[i]) {
				log_debug(nnti_debug_level, "Segment address mismatch (old[%d]=%p new[%d]=%p).",
						i, old_buf->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.buf, i, segments[i]);
				return(NNTI_EINVAL);
			}
			if (old_buf->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size != segment_lengths[i]) {
				log_debug(nnti_debug_level, "Segment length mismatch (old[%d]=%d new[%d]=%d).",
						i, old_buf->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size, i, segment_lengths[i]);
				return(NNTI_EINVAL);
			}
		}

        gni_mem_hdl=GNI_MEM_HDL(old_buf);
        gni_mem_hdl->ref_count++;
    }

    log_debug(nnti_ee_debug_level, "gni_mem_hdl->ref_count==%lu", gni_mem_hdl->ref_count);

	memset(reg_buf, 0, sizeof(NNTI_buffer_t));

    reg_buf->transport_id      = trans_hdl->id;
    reg_buf->buffer_owner      = trans_hdl->me;
    reg_buf->ops               = ops;
    reg_buf->payload_size=0;
    for (int i=0;i<num_segments;i++) {
        reg_buf->payload_size += segment_lengths[i];
    }
    reg_buf->payload           = (uint64_t)segments[0];
    reg_buf->transport_private = (uint64_t)gni_mem_hdl;

    log_debug(nnti_debug_level, "rpc_buffer->payload_size=%ld",
            reg_buf->payload_size);

    if (gni_mem_hdl->ref_count==1) {

    	gni_mem_hdl->extra=0;

        reg_buf->buffer_segments.NNTI_remote_addr_array_t_val=(NNTI_remote_addr_t *)calloc(num_segments, sizeof(NNTI_remote_addr_t));
        reg_buf->buffer_segments.NNTI_remote_addr_array_t_len=num_segments;

        gni_mem_hdl->mem_hdl_list =(gni_mem_handle_t *)calloc(num_segments, sizeof(gni_mem_handle_t));
        gni_mem_hdl->mem_hdl_count=num_segments;

    	if (ops == NNTI_GET_DST) {
    		gni_mem_hdl->type    =GET_DST_BUFFER;
    	} else if (ops == NNTI_GET_SRC) {
    		gni_mem_hdl->type    =GET_SRC_BUFFER;
    	} else if (ops == NNTI_PUT_SRC) {
    		gni_mem_hdl->type    =PUT_SRC_BUFFER;
    	} else if (ops == NNTI_PUT_DST) {
    		gni_mem_hdl->type    =PUT_DST_BUFFER;
    	} else if (ops == (NNTI_GET_SRC|NNTI_PUT_DST)) {
    		gni_mem_hdl->type    =RDMA_TARGET_BUFFER;
    	} else {
    		gni_mem_hdl->type=UNKNOWN_BUFFER;
    	}

        for (int i=0;i<num_segments;i++) {
        	gni_mem_hdl->mem_hdl_list[i]=register_memory_segment(
        			gni_mem_hdl->type,
        			segments[i],
        			segment_lengths[i],
        			gni_mem_hdl->extra);

        	log_debug(nnti_debug_level, "hdl->mem_hdl_list[%d]=(%llu,%llu)", i, (uint64_t)gni_mem_hdl->mem_hdl_list[i].qword1, (uint64_t)gni_mem_hdl->mem_hdl_list[i].qword2);
        }
    }

    if (config.use_rdma_target_ack) {
        if ((gni_mem_hdl->type == RDMA_TARGET_BUFFER) ||
                (gni_mem_hdl->type == GET_SRC_BUFFER) ||
                (gni_mem_hdl->type == PUT_DST_BUFFER)) {
            post_recv_work_request(reg_buf);
        }
    }

    if (nnti_rc==NNTI_OK) {
        for (int i=0;i<num_segments;i++) {
        	reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[i].transport_id                            = NNTI_TRANSPORT_GEMINI;
        	reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.mem_hdl.qword1 = gni_mem_hdl->mem_hdl_list[i].qword1;
        	reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.mem_hdl.qword2 = gni_mem_hdl->mem_hdl_list[i].qword2;

        	reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size = segment_lengths[i];
        	reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.buf  = (uint64_t)segments[i];
        	reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.type = NNTI_GNI_SEND_SRC;
        }
    }

    if (gni_mem_hdl->ref_count==1) {
        insert_buf_bufhash(reg_buf);
        log_debug(nnti_debug_level, "gni_mem_hdl->type==%llu",
                (uint64_t)gni_mem_hdl->type);
        log_debug(nnti_debug_level, "reg_buf.buf.hash==%llu",
                (uint64_t)hash6432shift(reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.buf));
        log_debug(nnti_debug_level, "ref_count==1 called insert_buf_bufhash() (reg_buf=%p, reg_buf.hash6432=%llu)",
                reg_buf, (uint64_t)hash6432shift(reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.buf));
        NNTI_buffer_t *tmp_buf=get_buf_bufhash((uint64_t)hash6432shift(reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.buf));
        log_debug(nnti_debug_level, "immediate get_buf_bufhash() says tmp_buf=%p", tmp_buf);
    }

nthread_unlock(&nnti_mem_lock);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "reg_buf",
                "end of NNTI_gni_register_memory", reg_buf);
    }

    log_debug(nnti_ee_debug_level, "exit (reg_buf=%p, reg_buf.hash6432=%llu)", reg_buf, (uint64_t)hash6432shift(reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.buf));

cleanup:
	switch(gni_rc) {
		case GNI_RC_SUCCESS:
			nnti_rc=NNTI_OK;
		default:
			nnti_rc=NNTI_EIO;
	}

	return (nnti_rc);
}


/**
 * @brief Cleanup after network operations are complete.
 *
 * Destroy an NNTI_buffer_t that was previously created by NNTI_regsiter_buffer().
 * It is the user's responsibility to release the the memory region.
 */
NNTI_result_t NNTI_gni_unregister_memory (
        NNTI_buffer_t    *reg_buf)
{
    NNTI_result_t rc=NNTI_OK;
    nnti_gni_memory_handle_t *gni_mem_hdl=NULL;

    assert(reg_buf);

    log_debug(nnti_ee_debug_level, "enter (reg_buf=%p)", reg_buf);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "reg_buf",
                "start of NNTI_gni_unregister_memory", reg_buf);
    }

nthread_lock(&nnti_mem_lock);

    gni_mem_hdl=(nnti_gni_memory_handle_t *)reg_buf->transport_private;
    assert(gni_mem_hdl);
    gni_mem_hdl->ref_count--;

    log_debug(nnti_ee_debug_level, "gni_mem_hdl->ref_count==%lu", gni_mem_hdl->ref_count);

    if (gni_mem_hdl->ref_count==0) {
        log_debug(nnti_ee_debug_level, "gni_mem_hdl->ref_count is 0.  release all resources.");

        if (gni_mem_hdl->type==REQUEST_BUFFER) {
            server_req_queue_destroy(
                    &transport_global_data.req_queue);

        } else {
        	for (int i=0;i<gni_mem_hdl->mem_hdl_count;i++) {
                unregister_memory(gni_mem_hdl->mem_hdl_list[i]);
        	}
        	if (gni_mem_hdl->mem_hdl_count > 1) {
        		free(gni_mem_hdl->mem_hdl_list);
        	}
        	gni_mem_hdl->mem_hdl_list =NULL;
        	gni_mem_hdl->mem_hdl_count=0;
        }

        del_buf_bufhash(reg_buf);

        nthread_lock(&gni_mem_hdl->wr_queue_lock);
        while (!gni_mem_hdl->wr_queue->empty()) {
            nnti_gni_work_request_t *wr=gni_mem_hdl->wr_queue->front();
            log_debug(nnti_debug_level, "removing pending wr=%p", wr);
            gni_mem_hdl->wr_queue->pop_front();
            del_wr_wrhash(wr);
            if (config.use_wr_pool) {
                if (wr->is_initiator==TRUE) {
                    wr_pool_initiator_push(wr);
                } else {
                    wr_pool_target_push(wr);
                }
            } else {
                log_debug(nnti_debug_level, "free(wr) (reg_buf=%p, wr=%p)", reg_buf, wr);
                free(wr);
            }
        }
        nthread_unlock(&gni_mem_hdl->wr_queue_lock);

        nthread_lock_fini(&gni_mem_hdl->wr_queue_lock);

        reg_buf->transport_id      = NNTI_TRANSPORT_NULL;
        GNI_SET_MATCH_ANY(&reg_buf->buffer_owner);
        reg_buf->ops               = (NNTI_buf_ops_t)0;
        reg_buf->payload_size      = 0;
        reg_buf->payload           = 0;
        reg_buf->transport_private = 0;

        if (gni_mem_hdl) {
            delete gni_mem_hdl->wr_queue;
            free(gni_mem_hdl);
        }
        if (reg_buf->buffer_segments.NNTI_remote_addr_array_t_val) {
            free(reg_buf->buffer_segments.NNTI_remote_addr_array_t_val);
            reg_buf->buffer_segments.NNTI_remote_addr_array_t_val=NULL;
            reg_buf->buffer_segments.NNTI_remote_addr_array_t_len=0;
        }

        memset(reg_buf, 0, sizeof(NNTI_buffer_t));
    }

nthread_unlock(&nnti_mem_lock);

    log_debug(nnti_ee_debug_level, "exit (reg_buf=%p)", reg_buf);

    return(rc);
}


/**
 * @brief Send a message to a peer.
 *
 * Send a message (<tt>msg_hdl</tt>) to a peer (<tt>peer_hdl</tt>).  It is expected that the
 * message is small, but the exact maximum size is transport dependent.
 */
NNTI_result_t NNTI_gni_send (
        const NNTI_peer_t   *peer_hdl,
        const NNTI_buffer_t *msg_hdl,
        const NNTI_buffer_t *dest_hdl,
        NNTI_work_request_t *wr)
{
    NNTI_result_t nnti_rc=NNTI_OK;

    trios_declare_timer(call_time);

    nnti_gni_memory_handle_t *gni_mem_hdl=NULL;
    nnti_gni_work_request_t  *gni_wr=NULL;

    log_debug(nnti_ee_debug_level, "enter");

    assert(peer_hdl);
    assert(msg_hdl);

    gni_mem_hdl=(nnti_gni_memory_handle_t *)msg_hdl->transport_private;
    assert(gni_mem_hdl);

    if (config.use_wr_pool) {
        gni_wr=wr_pool_initiator_pop();
    } else {
        gni_wr=(nnti_gni_work_request_t *)malloc(sizeof(nnti_gni_work_request_t));
        memset(gni_wr, 0, sizeof(nnti_gni_work_request_t));
        nthread_lock_init(&gni_wr->lock);
    }
    assert(gni_wr);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "msg_hdl",
                "NNTI_send", msg_hdl);
    }
    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "dest_hdl",
                "NNTI_send", dest_hdl);
    }

    if ((dest_hdl == NULL) || (dest_hdl->ops == NNTI_RECV_QUEUE)) {
        nnti_gni_connection_t *conn=get_conn_peer(peer_hdl);
        assert(conn);

        trios_start_timer(call_time);
        int rc=request_send(peer_hdl, &conn->queue_local_attrs, &conn->queue_remote_attrs, msg_hdl, 0, gni_wr);
        trios_stop_timer("send to request queue", call_time);

        if (rc == -1) {
            wr->result=NNTI_EDROPPED;
        } else {
            wr->result=NNTI_OK;
        }

    } else if (dest_hdl->ops == NNTI_RECV_DST) {
        nnti_gni_connection_t *conn=get_conn_peer(peer_hdl);
        assert(conn);

        trios_start_timer(call_time);
        buffer_send(peer_hdl, msg_hdl, dest_hdl, gni_wr);
        trios_stop_timer("send to receive dest", call_time);

        wr->result=NNTI_OK;
    }

    gni_wr->nnti_wr=wr;

    wr->transport_id     =msg_hdl->transport_id;
    wr->reg_buf          =(NNTI_buffer_t*)msg_hdl;
    wr->ops              =NNTI_SEND_SRC;
    wr->transport_private=(uint64_t)gni_wr;

    log_debug(nnti_debug_level, "sending to (%s, instance=%llu)", peer_hdl->url, (uint64_t)peer_hdl->peer.NNTI_remote_process_t_u.gni.inst_id);

    log_debug(nnti_ee_debug_level, "exit");

    return(nnti_rc);
}


/**
 * @brief Transfer data to a peer.
 *
 * Put the contents of <tt>src_buffer_hdl</tt> into <tt>dest_buffer_hdl</tt>.  It is
 * assumed that the destination is at least <tt>src_length</tt> bytes in size.
 *
 */
NNTI_result_t NNTI_gni_put (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset,
        NNTI_work_request_t *wr)
{
    int rc=NNTI_OK;
    trios_declare_timer(call_time);

    gni_ep_handle_t        ep_hdl;

    nnti_gni_memory_handle_t *gni_mem_hdl=NULL;
    nnti_gni_work_request_t  *gni_wr=NULL;

    uint64_t remaining=0;

    bool use_fma=true;

    struct {
    	uint32_t segment_index;     // index of the current segment
    	uint64_t segment_offset;    // offset into the current segment
    	uint32_t segment_remaining; // bytes remaining in the current segment
    	uint32_t last_segment;      // index of the last segment in this transfer
    } src_params, dst_params, src_copy, dst_copy;

    uint32_t sge_index=0;

    log_debug(nnti_ee_debug_level, "enter");

    assert(src_buffer_hdl);
    assert(dest_buffer_hdl);

    if (src_offset > src_buffer_hdl->payload_size) {
    }
    if (src_buffer_hdl->payload_size < src_length-src_offset) {
    }
    if (dest_offset > dest_buffer_hdl->payload_size) {
    }
    if (dest_buffer_hdl->payload_size < src_length-dest_offset) {
    }

    gni_mem_hdl=(nnti_gni_memory_handle_t *)src_buffer_hdl->transport_private;
    assert(gni_mem_hdl);

	if (config.use_wr_pool) {
		gni_wr=wr_pool_initiator_pop();
	} else {
		gni_wr=(nnti_gni_work_request_t *)malloc(sizeof(nnti_gni_work_request_t));
		memset(gni_wr, 0, sizeof(nnti_gni_work_request_t));
		nthread_lock_init(&gni_wr->lock);
	}
	assert(gni_wr);

	gni_wr->reg_buf = src_buffer_hdl;

    nnti_gni_connection_t *conn=get_conn_peer(&dest_buffer_hdl->buffer_owner);
    assert(conn);
    ep_hdl=conn->ep_hdl;

	nthread_lock(&gni_mem_hdl->wr_queue_lock);
	gni_mem_hdl->wr_queue->push_back(gni_wr);
	nthread_unlock(&gni_mem_hdl->wr_queue_lock);

	// the src_offset could exclude some local segments from this operation.  find the first segment.
	src_params.segment_offset   =src_offset;
	src_params.segment_remaining=0;
	src_params.segment_index    =0;
	for (uint32_t i=0;i<src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_len;i++) {
		log_debug(nnti_debug_level, "src_offset=%llu, src_params.segment_offset=%llu, src_params.segment_remaining=%llu, src_buffer.segment[%lu].size==%lu",
				src_offset, src_params.segment_offset, src_params.segment_remaining,
				i, src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size);
		if (src_params.segment_offset > src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size) {
			src_params.segment_offset -= src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size;

			log_debug(nnti_debug_level, "src_params.segment_index=%lu, src_params.segment_offset=%llu, src_params.segment_remaining=%llu",
					src_params.segment_index, src_params.segment_offset, src_params.segment_remaining);

			continue;
		} else {
			src_params.segment_index     = i;
			src_params.segment_remaining = src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size - src_params.segment_offset;

			log_debug(nnti_debug_level, "src_params.segment_index=%lu, src_params.segment_offset=%llu, src_params.segment_remaining=%llu",
					src_params.segment_index, src_params.segment_offset, src_params.segment_remaining);

			break;
		}
	}
	remaining=src_length;
	for (uint32_t i=src_params.segment_index;i<src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_len;i++) {
		if (remaining > src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size) {
			remaining -= src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size;
		} else {
			src_params.last_segment = i;
			break;
		}
	}

	log_debug(nnti_debug_level, "src_params.segment_index=%lu, src_params.segment_offset=%llu, src_params.segment_remaining=%llu, src_params.last_segment=%lu",
			src_params.segment_index, src_params.segment_offset, src_params.segment_remaining, src_params.last_segment);

	// the dest_offset could exclude some local segments from this operation.  find the first segment.
	dst_params.segment_offset   =dest_offset;
	dst_params.segment_remaining=0;
	dst_params.segment_index    =0;
	for (uint32_t i=0;i<dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_len;i++) {
		log_debug(nnti_debug_level, "dest_offset=%llu, dst_params.segment_offset=%llu, dst_params.segment_remaining=%llu, dest_buffer.segment[%lu].size==%lu",
				dest_offset, dst_params.segment_offset, dst_params.segment_remaining,
				i, dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size);
		if (dst_params.segment_offset > dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size) {
			dst_params.segment_offset -= dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size;

			log_debug(nnti_debug_level, "dst_params.segment_index=%lu, dst_params.segment_offset=%llu, dst_params.segment_remaining=%llu",
					dst_params.segment_index, dst_params.segment_offset, dst_params.segment_remaining);

			continue;
		} else {
			dst_params.segment_index     = i;
			dst_params.segment_remaining = dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size - dst_params.segment_offset;

			log_debug(nnti_debug_level, "dst_params.segment_index=%lu, dst_params.segment_offset=%llu, dst_params.segment_remaining=%llu",
					dst_params.segment_index, dst_params.segment_offset, dst_params.segment_remaining);

			break;
		}
	}
	remaining=src_length;
	for (uint32_t i=dst_params.segment_index;i<dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_len;i++) {
		if (remaining > dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size) {
			remaining -= dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size;
		} else {
			dst_params.last_segment = i;
			break;
		}
	}

	log_debug(nnti_debug_level, "dst_params.segment_index=%lu, dst_params.segment_offset=%llu, dst_params.segment_remaining=%llu, dst_params.last_segment=%lu",
			dst_params.segment_index, dst_params.segment_offset, dst_params.segment_remaining, dst_params.last_segment);

	/* START calculate the number of SGEs required */
	{
	src_copy=src_params;
	dst_copy=dst_params;

	NNTI_remote_addr_t *src_remote_addr=&src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[src_copy.segment_index];
	NNTI_remote_addr_t *dst_remote_addr=&dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[dst_copy.segment_index];

	gni_wr->sge_count=0;

	while (src_copy.segment_index <= src_copy.last_segment) {

		gni_wr->sge_count++;

		if (src_copy.segment_remaining < dst_copy.segment_remaining) {
			log_debug(nnti_debug_level, "the remaining source segment fits in the remaining destination segment with extra space in the destination");
			// the remaining source segment fits in the remaining destination segment with extra space
			dst_copy.segment_offset    += src_copy.segment_remaining;
			dst_copy.segment_remaining -= src_copy.segment_remaining;

			src_remote_addr++;
			src_copy.segment_index++;
			src_copy.segment_offset=0;
			src_copy.segment_remaining=src_remote_addr->NNTI_remote_addr_t_u.gni.size;

		} else if (src_copy.segment_remaining == dst_copy.segment_remaining) {
			log_debug(nnti_debug_level, "the remaining source segment fits in the remaining destination segment EXACTLY");
			// the remaining source segment fits EXACTLY in the remaining destination segment
			src_remote_addr++;
			src_copy.segment_index++;
			src_copy.segment_offset=0;
			src_copy.segment_remaining=src_remote_addr->NNTI_remote_addr_t_u.gni.size;

			dst_remote_addr++;
			dst_copy.segment_index++;
			dst_copy.segment_offset=0;
			dst_copy.segment_remaining=dst_remote_addr->NNTI_remote_addr_t_u.gni.size;

		} else {
			log_debug(nnti_debug_level, "the remaining source segment DOES NOT fit in the remaining destination segment");
			// the remaining source segment DOES NOT fit in the remaining destination segment
			src_copy.segment_offset    += dst_copy.segment_remaining;
			src_copy.segment_remaining -= dst_copy.segment_remaining;

			dst_remote_addr++;
			dst_copy.segment_index++;
			dst_copy.segment_offset=0;
			dst_copy.segment_remaining=dst_remote_addr->NNTI_remote_addr_t_u.gni.size;

		}

		log_debug(nnti_debug_level, "src_copy.segment_index=%lu, src_copy.segment_offset=%llu, src_copy.segment_remaining=%lu, src_copy.last_segment=%lu",
				src_copy.segment_index, src_copy.segment_offset, src_copy.segment_remaining, src_copy.last_segment);
		log_debug(nnti_debug_level, "dst_copy.segment_index=%lu, dst_copy.segment_offset=%llu, dst_copy.segment_remaining=%lu, dst_copy.last_segment=%lu",
				dst_copy.segment_index, dst_copy.segment_offset, dst_copy.segment_remaining, dst_copy.last_segment);

	}
	}
	/* END   calculate the number of SGEs required */

	NNTI_remote_addr_t *src_remote_addr=&src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[src_params.segment_index];
	NNTI_remote_addr_t *dst_remote_addr=&dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[dst_params.segment_index];

	gni_wr->sge_list  = (nnti_gni_sge_t *)malloc(gni_wr->sge_count*sizeof(nnti_gni_sge_t));

	log_debug(nnti_debug_level, "this get requires %d SGEs", gni_wr->sge_count);

	sge_index=0;

	while (src_params.segment_index <= src_params.last_segment) {

		memset(&gni_wr->sge_list[sge_index].post_desc, 0, sizeof(gni_post_descriptor_t));

		gni_wr->sge_list[sge_index].post_desc.src_cq_hndl=transport_global_data.ep_cq_hdl;

		gni_wr->sge_list[sge_index].post_desc.local_addr            =src_remote_addr->NNTI_remote_addr_t_u.gni.buf+src_params.segment_offset;
		gni_wr->sge_list[sge_index].post_desc.local_mem_hndl.qword1 =src_remote_addr->NNTI_remote_addr_t_u.gni.mem_hdl.qword1;
		gni_wr->sge_list[sge_index].post_desc.local_mem_hndl.qword2 =src_remote_addr->NNTI_remote_addr_t_u.gni.mem_hdl.qword2;

		gni_wr->sge_list[sge_index].post_desc.remote_addr           =dst_remote_addr->NNTI_remote_addr_t_u.gni.buf+dst_params.segment_offset;
		gni_wr->sge_list[sge_index].post_desc.remote_mem_hndl.qword1=dst_remote_addr->NNTI_remote_addr_t_u.gni.mem_hdl.qword1;
		gni_wr->sge_list[sge_index].post_desc.remote_mem_hndl.qword2=dst_remote_addr->NNTI_remote_addr_t_u.gni.mem_hdl.qword2;

		if (src_params.segment_remaining <= dst_params.segment_remaining) {
			// the remaining source segment fits in the remaining destination segment with extra space
			gni_wr->sge_list[sge_index].post_desc.length=src_params.segment_remaining;
		} else {
			// the remaining source segment DOES NOT fit in the remaining destination segment
			gni_wr->sge_list[sge_index].post_desc.length=dst_params.segment_remaining;
		}

		set_post_desc(&gni_wr->sge_list[sge_index].post_desc, PUT_DST_BUFFER, gni_wr->sge_list[sge_index].post_desc.length);

		gni_wr->sge_list[sge_index].state=NNTI_GNI_SGE_STATE_STARTED;

		gni_wr->sge_list[sge_index].gni_wr=gni_wr;

		gni_wr->wc=&gni_wr->rdma_wc;
		gni_wr->wc->op         =GNI_OP_PUT_TARGET;
		gni_wr->wc->byte_len   =gni_wr->sge_list[sge_index].post_desc.length;
		gni_wr->wc->src_index  =src_params.segment_index;
		gni_wr->wc->src_offset =src_params.segment_offset;
		gni_wr->wc->dest_index =dst_params.segment_index;
		gni_wr->wc->dest_offset=dst_params.segment_offset;
		gni_wr->wc->inst_id    =transport_global_data.instance;

		print_wc(gni_wr->wc);

		if (config.rdma_mode==RDMA_CROSSOVER) {
			if (gni_wr->sge_list[sge_index].post_desc.length > config.fma_bte_crossover_size) {
				use_fma=false;
			} else {
				use_fma=true;
			}
		} else if ((config.rdma_mode==RDMA_BTE) || (config.rdma_mode==RDMA_MIXED)) {
			use_fma=false;
		} else if (config.rdma_mode==RDMA_FMA) {
			use_fma=true;
		}

		insert_sge_sgehash(&gni_wr->sge_list[sge_index]);

	    nthread_lock(&nnti_gni_lock);
	    GNI_EpSetEventData(
	            ep_hdl,
	            hash6432shift((uint64_t)&gni_wr->sge_list[sge_index]),
	            hash6432shift((uint64_t)dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.buf));

	    if (use_fma) {
	        log_debug(nnti_event_debug_level, "calling PostFma(fma put ; ep_hdl(%llu) transport_global_data.ep_cq_hdl(%llu) local_mem_hdl(%llu, %llu) remote_mem_hdl(%llu, %llu))",
	                ep_hdl, transport_global_data.ep_cq_hdl,
	                gni_wr->sge_list[sge_index].post_desc.local_mem_hndl.qword1, gni_wr->sge_list[sge_index].post_desc.local_mem_hndl.qword2,
	                gni_wr->sge_list[sge_index].post_desc.remote_mem_hndl.qword1, gni_wr->sge_list[sge_index].post_desc.remote_mem_hndl.qword2);
	        trios_start_timer(call_time);
	        rc=GNI_PostFma(ep_hdl, &gni_wr->sge_list[sge_index].post_desc);
	        trios_stop_timer("PostFma put", call_time);
	        if (rc!=GNI_RC_SUCCESS) {
	            log_error(nnti_debug_level, "failed to post FMA (rc=%d): %s", rc, strerror(errno));
	            rc=NNTI_EIO;
	        }
	    } else {
	        log_debug(nnti_event_debug_level, "calling PostRdma(rdma put ; ep_hdl(%llu) transport_global_data.ep_cq_hdl(%llu) local_mem_hdl(%llu, %llu) remote_mem_hdl(%llu, %llu))",
	                ep_hdl, transport_global_data.ep_cq_hdl,
	                gni_wr->sge_list[sge_index].post_desc.local_mem_hndl.qword1, gni_wr->sge_list[sge_index].post_desc.local_mem_hndl.qword2,
	                gni_wr->sge_list[sge_index].post_desc.remote_mem_hndl.qword1, gni_wr->sge_list[sge_index].post_desc.remote_mem_hndl.qword2);
	        trios_start_timer(call_time);
	        rc=GNI_PostRdma(ep_hdl, &gni_wr->sge_list[sge_index].post_desc);
	        trios_stop_timer("PostRdma put", call_time);
	        if (rc!=GNI_RC_SUCCESS) {
	            log_error(nnti_debug_level, "failed to post RDMA (rc=%d): %s", rc, strerror(errno));
	            rc=NNTI_EIO;
	        }
	    }
	    nthread_unlock(&nnti_gni_lock);

		if (src_params.segment_remaining < dst_params.segment_remaining) {
			// the remaining source segment fits in the remaining destination segment with extra space
			dst_params.segment_offset    += src_params.segment_remaining;
			dst_params.segment_remaining -= src_params.segment_remaining;

			src_remote_addr++;
			src_params.segment_index++;
			src_params.segment_offset=0;
			src_params.segment_remaining=src_remote_addr->NNTI_remote_addr_t_u.gni.size;

		} else if (src_params.segment_remaining == dst_params.segment_remaining) {
			// the remaining source segment fits EXACTLY in the remaining destination segment
			src_remote_addr++;
			src_params.segment_index++;
			src_params.segment_offset=0;
			src_params.segment_remaining=src_remote_addr->NNTI_remote_addr_t_u.gni.size;

			dst_remote_addr++;
			dst_params.segment_index++;
			dst_params.segment_offset=0;
			dst_params.segment_remaining=dst_remote_addr->NNTI_remote_addr_t_u.gni.size;

		} else {
			// the remaining source segment DOES NOT fit in the remaining destination segment
			src_params.segment_offset    += dst_params.segment_remaining;
			src_params.segment_remaining -= dst_params.segment_remaining;

			dst_remote_addr++;
			dst_params.segment_index++;
			dst_params.segment_offset=0;
			dst_params.segment_remaining=dst_remote_addr->NNTI_remote_addr_t_u.gni.size;

		}

		sge_index++;
	}

    gni_wr->last_op      =GNI_OP_PUT_INITIATOR;
    gni_wr->peer_instance=dest_buffer_hdl->buffer_owner.peer.NNTI_remote_process_t_u.gni.inst_id;
    log_debug(nnti_event_debug_level, "gni_wr->peer_instance==%lu", gni_wr->peer_instance);

    gni_wr->nnti_wr=wr;

    wr->transport_id     =src_buffer_hdl->transport_id;
    wr->reg_buf          =(NNTI_buffer_t*)src_buffer_hdl;
    wr->ops              =NNTI_PUT_SRC;
    wr->result           =NNTI_OK;
    wr->transport_private=(uint64_t)gni_wr;

    log_debug(nnti_ee_debug_level, "exit");

    return((NNTI_result_t)rc);
}


/**
 * @brief Transfer data from a peer.
 *
 * Get the contents of <tt>src_buffer_hdl</tt> into <tt>dest_buffer_hdl</tt>.  It is
 * assumed that the destination is at least <tt>src_length</tt> bytes in size.
 *
 */
NNTI_result_t NNTI_gni_get (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset,
        NNTI_work_request_t *wr)
{
    int rc=NNTI_OK;
    trios_declare_timer(call_time);

    gni_ep_handle_t        ep_hdl;

    nnti_gni_memory_handle_t *gni_mem_hdl=NULL;
    nnti_gni_work_request_t  *gni_wr=NULL;

    uint64_t remaining=0;

    bool use_fma=true;

    struct {
    	uint32_t segment_index;     // index of the current segment
    	uint64_t segment_offset;    // offset into the current segment
    	uint64_t segment_remaining; // bytes remaining in the current segment
    	uint32_t last_segment;      // index of the last segment in this transfer
    } src_params, dst_params, src_copy, dst_copy;

    uint32_t sge_index=0;

    log_debug(nnti_ee_debug_level, "enter");

    assert(src_buffer_hdl);
    assert(dest_buffer_hdl);

    if (src_offset > src_buffer_hdl->payload_size) {
    }
    if (src_buffer_hdl->payload_size < src_length-src_offset) {
    }
    if (dest_offset > dest_buffer_hdl->payload_size) {
    }
    if (dest_buffer_hdl->payload_size < src_length-dest_offset) {
    }

    gni_mem_hdl=(nnti_gni_memory_handle_t *)dest_buffer_hdl->transport_private;
    assert(gni_mem_hdl);

    if (config.use_wr_pool) {
        gni_wr=wr_pool_initiator_pop();
    } else {
        gni_wr=(nnti_gni_work_request_t *)malloc(sizeof(nnti_gni_work_request_t));
        memset(gni_wr, 0, sizeof(nnti_gni_work_request_t));
        nthread_lock_init(&gni_wr->lock);
    }
    assert(gni_wr);

    gni_wr->reg_buf = dest_buffer_hdl;

    nnti_gni_connection_t *conn=get_conn_peer(&src_buffer_hdl->buffer_owner);
    assert(conn);
    ep_hdl=conn->ep_hdl;

	nthread_lock(&gni_mem_hdl->wr_queue_lock);
	gni_mem_hdl->wr_queue->push_back(gni_wr);
	nthread_unlock(&gni_mem_hdl->wr_queue_lock);

	// the dest_offset could exclude some local segments from this operation.  find the first segment.
	dst_params.segment_offset   =dest_offset;
	dst_params.segment_remaining=0;
	dst_params.segment_index    =0;
	for (uint32_t i=0;i<dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_len;i++) {
		log_debug(nnti_debug_level, "dest_offset=%llu, dst_params.segment_offset=%llu, dst_params.segment_remaining=%llu, dest_buffer.segment[%llu].size==%llu",
				dest_offset, dst_params.segment_offset, dst_params.segment_remaining,
				(uint64_t)i, (uint64_t)dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size);
		if (dst_params.segment_offset > dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size) {
			dst_params.segment_offset -= dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size;

			log_debug(nnti_debug_level, "dst_params.segment_index=%lu, dst_params.segment_offset=%llu, dst_params.segment_remaining=%llu",
					dst_params.segment_index, dst_params.segment_offset, dst_params.segment_remaining);

			continue;
		} else {
			dst_params.segment_index     = i;
			dst_params.segment_remaining = dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size - dst_params.segment_offset;

			log_debug(nnti_debug_level, "dst_params.segment_index=%llu, dst_params.segment_offset=%llu, dst_params.segment_remaining=%llu",
					(uint64_t)dst_params.segment_index, dst_params.segment_offset, dst_params.segment_remaining);

			break;
		}
	}
	remaining=src_length;
	for (uint32_t i=dst_params.segment_index;i<dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_len;i++) {
		if (remaining > dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size) {
			remaining -= dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size;
		} else {
			dst_params.last_segment = i;
			break;
		}
	}

	log_debug(nnti_debug_level, "dst_params.segment_index=%llu, dst_params.segment_offset=%llu, dst_params.segment_remaining=%llu, dst_params.last_segment=%llu",
			(uint64_t)dst_params.segment_index, dst_params.segment_offset, dst_params.segment_remaining, (uint64_t)dst_params.last_segment);

	// the src_offset could exclude some local segments from this operation.  find the first segment.
	src_params.segment_offset   =src_offset;
	src_params.segment_remaining=0;
	src_params.segment_index    =0;
	for (uint32_t i=0;i<src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_len;i++) {
		log_debug(nnti_debug_level, "src_offset=%llu, src_params.segment_offset=%llu, src_params.segment_remaining=%llu, src_buffer.segment[%llu].size==%llu",
				src_offset, src_params.segment_offset, src_params.segment_remaining,
				(uint64_t)i, (uint64_t)src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size);
		if (src_params.segment_offset > src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size) {
			src_params.segment_offset -= src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size;

			log_debug(nnti_debug_level, "src_params.segment_index=%llu, src_params.segment_offset=%llu, src_params.segment_remaining=%llu",
					(uint64_t)src_params.segment_index, src_params.segment_offset, src_params.segment_remaining);

			continue;
		} else {
			src_params.segment_index     = i;
			src_params.segment_remaining = src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size - src_params.segment_offset;

			log_debug(nnti_debug_level, "src_params.segment_index=%llu, src_params.segment_offset=%llu, src_params.segment_remaining=%llu",
					(uint64_t)src_params.segment_index, src_params.segment_offset, src_params.segment_remaining);

			break;
		}
	}
	remaining=src_length;
	for (uint32_t i=src_params.segment_index;i<src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_len;i++) {
		if (remaining > src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size) {
			remaining -= src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.gni.size;
		} else {
			src_params.last_segment = i;
			break;
		}
	}

	log_debug(nnti_debug_level, "src_params.segment_index=%llu, src_params.segment_offset=%llu, src_params.segment_remaining=%llu, src_params.last_segment=%llu",
			(uint64_t)src_params.segment_index, src_params.segment_offset, src_params.segment_remaining, (uint64_t)src_params.last_segment);

	/* START calculate the number of SGEs required */
	{
	dst_copy=dst_params;
	src_copy=src_params;

	NNTI_remote_addr_t *dst_remote_addr=&dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[dst_copy.segment_index];
	NNTI_remote_addr_t *src_remote_addr=&src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[src_copy.segment_index];

	gni_wr->sge_count=0;

	while (src_copy.segment_index <= src_copy.last_segment) {

		gni_wr->sge_count++;

		if (src_copy.segment_remaining < dst_copy.segment_remaining) {
			log_debug(nnti_debug_level, "the remaining source segment fits in the remaining destination segment with extra space in the destination");
			// the remaining source segment fits in the remaining destination segment with extra space
			dst_copy.segment_offset    += src_copy.segment_remaining;
			dst_copy.segment_remaining -= src_copy.segment_remaining;

			src_remote_addr++;
			src_copy.segment_index++;
			src_copy.segment_offset=0;
			src_copy.segment_remaining=src_remote_addr->NNTI_remote_addr_t_u.gni.size;

		} else if (src_copy.segment_remaining == dst_copy.segment_remaining) {
			log_debug(nnti_debug_level, "the remaining source segment fits in the remaining destination segment EXACTLY");
			// the remaining source segment fits EXACTLY in the remaining destination segment
			src_remote_addr++;
			src_copy.segment_index++;
			src_copy.segment_offset=0;
			src_copy.segment_remaining=src_remote_addr->NNTI_remote_addr_t_u.gni.size;

			dst_remote_addr++;
			dst_copy.segment_index++;
			dst_copy.segment_offset=0;
			dst_copy.segment_remaining=dst_remote_addr->NNTI_remote_addr_t_u.gni.size;

		} else {
			log_debug(nnti_debug_level, "the remaining source segment DOES NOT fit in the remaining destination segment");
			// the remaining source segment DOES NOT fit in the remaining destination segment
			src_copy.segment_offset    += dst_copy.segment_remaining;
			src_copy.segment_remaining -= dst_copy.segment_remaining;

			dst_remote_addr++;
			dst_copy.segment_index++;
			dst_copy.segment_offset=0;
			dst_copy.segment_remaining=dst_remote_addr->NNTI_remote_addr_t_u.gni.size;

		}

		log_debug(nnti_debug_level, "src_copy.segment_index=%llu, src_copy.segment_offset=%llu, src_copy.segment_remaining=%llu, src_copy.last_segment=%llu",
				(uint64_t)src_copy.segment_index, src_copy.segment_offset, (uint64_t)src_copy.segment_remaining, (uint64_t)src_copy.last_segment);
		log_debug(nnti_debug_level, "dst_copy.segment_index=%llu, dst_copy.segment_offset=%llu, dst_copy.segment_remaining=%llu, dst_copy.last_segment=%llu",
				(uint64_t)dst_copy.segment_index, dst_copy.segment_offset, (uint64_t)dst_copy.segment_remaining, (uint64_t)dst_copy.last_segment);
	}
	}
	/* END   calculate the number of SGEs required */

	NNTI_remote_addr_t *dst_remote_addr=&dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[dst_params.segment_index];
	NNTI_remote_addr_t *src_remote_addr=&src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[src_params.segment_index];

	gni_wr->sge_list  = (nnti_gni_sge_t *)malloc(gni_wr->sge_count*sizeof(nnti_gni_sge_t));

	log_debug(nnti_debug_level, "this get requires %d SGEs", gni_wr->sge_count);

	sge_index=0;

	while (src_params.segment_index <= src_params.last_segment) {

		memset(&gni_wr->sge_list[sge_index].post_desc, 0, sizeof(gni_post_descriptor_t));

		gni_wr->sge_list[sge_index].post_desc.src_cq_hndl=transport_global_data.ep_cq_hdl;

		gni_wr->sge_list[sge_index].post_desc.local_addr            =dst_remote_addr->NNTI_remote_addr_t_u.gni.buf+dst_params.segment_offset;
		gni_wr->sge_list[sge_index].post_desc.local_mem_hndl.qword1 =dst_remote_addr->NNTI_remote_addr_t_u.gni.mem_hdl.qword1;
		gni_wr->sge_list[sge_index].post_desc.local_mem_hndl.qword2 =dst_remote_addr->NNTI_remote_addr_t_u.gni.mem_hdl.qword2;

		gni_wr->sge_list[sge_index].post_desc.remote_addr           =src_remote_addr->NNTI_remote_addr_t_u.gni.buf+src_params.segment_offset;
		gni_wr->sge_list[sge_index].post_desc.remote_mem_hndl.qword1=src_remote_addr->NNTI_remote_addr_t_u.gni.mem_hdl.qword1;
		gni_wr->sge_list[sge_index].post_desc.remote_mem_hndl.qword2=src_remote_addr->NNTI_remote_addr_t_u.gni.mem_hdl.qword2;

		if (src_params.segment_remaining <= dst_params.segment_remaining) {
			// the remaining source segment fits in the remaining destination segment with extra space
			gni_wr->sge_list[sge_index].post_desc.length=src_params.segment_remaining;
		} else {
			// the remaining source segment DOES NOT fit in the remaining destination segment
			gni_wr->sge_list[sge_index].post_desc.length=dst_params.segment_remaining;
		}

		set_post_desc(&gni_wr->sge_list[sge_index].post_desc, GET_SRC_BUFFER, gni_wr->sge_list[sge_index].post_desc.length);

		gni_wr->sge_list[sge_index].state=NNTI_GNI_SGE_STATE_STARTED;

		gni_wr->sge_list[sge_index].gni_wr=gni_wr;

		gni_wr->wc=&gni_wr->rdma_wc;
		gni_wr->wc->op         =GNI_OP_GET_TARGET;
		gni_wr->wc->byte_len   =gni_wr->sge_list[sge_index].post_desc.length;
		gni_wr->wc->src_index  =src_params.segment_index;
		gni_wr->wc->src_offset =src_params.segment_offset;
		gni_wr->wc->dest_index =dst_params.segment_index;
		gni_wr->wc->dest_offset=dst_params.segment_offset;
		gni_wr->wc->inst_id    =transport_global_data.instance;

		if (config.rdma_mode==RDMA_CROSSOVER) {
			if (gni_wr->sge_list[sge_index].post_desc.length > config.fma_bte_crossover_size) {
				use_fma=false;
			} else {
				use_fma=true;
			}
		} else if ((config.rdma_mode==RDMA_BTE) || (config.rdma_mode==RDMA_MIXED)) {
			use_fma=false;
		} else if (config.rdma_mode==RDMA_FMA) {
			use_fma=true;
		}

		insert_sge_sgehash(&gni_wr->sge_list[sge_index]);

	    nthread_lock(&nnti_gni_lock);
	    GNI_EpSetEventData(
	            ep_hdl,
	            hash6432shift((uint64_t)&gni_wr->sge_list[sge_index]),
	            hash6432shift((uint64_t)src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.buf));

	    if (use_fma) {
	        log_debug(nnti_event_debug_level, "calling PostFma(fma get ; ep_hdl(%llu) transport_global_data.ep_cq_hdl(%llu) local_mem_hdl(%llu, %llu) remote_mem_hdl(%llu, %llu))",
	                ep_hdl, transport_global_data.ep_cq_hdl,
	                gni_wr->sge_list[sge_index].post_desc.local_mem_hndl.qword1, gni_wr->sge_list[sge_index].post_desc.local_mem_hndl.qword2,
	                gni_wr->sge_list[sge_index].post_desc.remote_mem_hndl.qword1, gni_wr->sge_list[sge_index].post_desc.remote_mem_hndl.qword2);
	        trios_start_timer(call_time);
	        rc=GNI_PostFma(ep_hdl, &gni_wr->sge_list[sge_index].post_desc);
	        trios_stop_timer("PostFma get", call_time);
	        if (rc!=GNI_RC_SUCCESS) {
	            log_error(nnti_debug_level, "failed to post FMA (rc=%d): %s", rc, strerror(errno));
	            rc=NNTI_EIO;
	        }
	    } else {
	        log_debug(nnti_event_debug_level, "calling PostRdma(rdma get ; ep_hdl(%llu) transport_global_data.ep_cq_hdl(%llu) local_mem_hdl(%llu, %llu) remote_mem_hdl(%llu, %llu))",
	                ep_hdl, transport_global_data.ep_cq_hdl,
	                gni_wr->sge_list[sge_index].post_desc.local_mem_hndl.qword1, gni_wr->sge_list[sge_index].post_desc.local_mem_hndl.qword2,
	                gni_wr->sge_list[sge_index].post_desc.remote_mem_hndl.qword1, gni_wr->sge_list[sge_index].post_desc.remote_mem_hndl.qword2);
	        trios_start_timer(call_time);
	        rc=GNI_PostRdma(ep_hdl, &gni_wr->sge_list[sge_index].post_desc);
	        trios_stop_timer("PostRdma get", call_time);
	        if (rc!=GNI_RC_SUCCESS) {
	            log_error(nnti_debug_level, "failed to post RDMA (rc=%d): %s", rc, strerror(errno));
	            rc=NNTI_EIO;
	        }
	    }
	    nthread_unlock(&nnti_gni_lock);

		if (src_params.segment_remaining < dst_params.segment_remaining) {
			// the remaining source segment fits in the remaining destination segment with extra space
			dst_params.segment_offset    += src_params.segment_remaining;
			dst_params.segment_remaining -= src_params.segment_remaining;

			src_remote_addr++;
			src_params.segment_index++;
			src_params.segment_offset=0;
			src_params.segment_remaining=src_remote_addr->NNTI_remote_addr_t_u.gni.size;

		} else if (src_params.segment_remaining == dst_params.segment_remaining) {
			// the remaining source segment fits EXACTLY in the remaining destination segment
			src_remote_addr++;
			src_params.segment_index++;
			src_params.segment_offset=0;
			src_params.segment_remaining=src_remote_addr->NNTI_remote_addr_t_u.gni.size;

			dst_remote_addr++;
			dst_params.segment_index++;
			dst_params.segment_offset=0;
			dst_params.segment_remaining=dst_remote_addr->NNTI_remote_addr_t_u.gni.size;

		} else {
			// the remaining source segment DOES NOT fit in the remaining destination segment
			src_params.segment_offset    += dst_params.segment_remaining;
			src_params.segment_remaining -= dst_params.segment_remaining;

			dst_remote_addr++;
			dst_params.segment_index++;
			dst_params.segment_offset=0;
			dst_params.segment_remaining=dst_remote_addr->NNTI_remote_addr_t_u.gni.size;

		}

		sge_index++;
	}

    gni_wr->last_op  =GNI_OP_GET_INITIATOR;
    gni_wr->peer_instance=src_buffer_hdl->buffer_owner.peer.NNTI_remote_process_t_u.gni.inst_id;
    log_debug(nnti_event_debug_level, "gni_wr->peer_instance==%lu", gni_wr->peer_instance);

    gni_wr->nnti_wr=wr;

    wr->transport_id     =dest_buffer_hdl->transport_id;
    wr->reg_buf          =(NNTI_buffer_t*)dest_buffer_hdl;
    wr->ops              =NNTI_GET_SRC;
    wr->result           =NNTI_OK;
    wr->transport_private=(uint64_t)gni_wr;

    log_debug(nnti_ee_debug_level, "exit");

    return((NNTI_result_t)rc);
}


/**
 * @brief Transfer data to a peer.
 *
 * \param[in] src_buffer_hdl    A buffer containing the data to put.
 * \param[in] src_length        The number of bytes to put.
 * \param[in] dest_buffer_list  A list of buffers to put the data into.
 * \param[in] dest_count        The number of destination buffers.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_gni_scatter (
        const NNTI_buffer_t  *src_buffer_hdl,
        const uint64_t        src_length,
        const NNTI_buffer_t **dest_buffer_list,
        const uint64_t        dest_count,
        NNTI_work_request_t  *wr)
{
    log_debug(nnti_debug_level, "enter");

    log_debug(nnti_debug_level, "exit");

    return NNTI_ENOTSUP;
}


/**
 * @brief Transfer data from a peer.
 *
 * \param[in] src_buffer_list  A list of buffers containing the data to get.
 * \param[in] src_length       The number of bytes to get.
 * \param[in] src_count        The number of source buffers.
 * \param[in] dest_buffer_hdl  A buffer to get the data into.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_gni_gather (
        const NNTI_buffer_t **src_buffer_list,
        const uint64_t        src_length,
        const uint64_t        src_count,
        const NNTI_buffer_t  *dest_buffer_hdl,
        NNTI_work_request_t  *wr)
{
    log_debug(nnti_debug_level, "enter");

    log_debug(nnti_debug_level, "exit");

    return NNTI_ENOTSUP;
}


NNTI_result_t NNTI_gni_atomic_set_callback (
		const NNTI_transport_t *trans_hdl,
		const uint64_t          local_atomic,
		NNTI_callback_fn_t      cbfunc,
		void                   *context)
{
    return NNTI_ENOTSUP;
}


NNTI_result_t NNTI_gni_atomic_read (
		const NNTI_transport_t *trans_hdl,
		const uint64_t          local_atomic,
		int64_t                *value)
{
    nthread_lock(&transport_global_data.atomics_lock);
    *value = transport_global_data.atomics[local_atomic];
    nthread_unlock(&transport_global_data.atomics_lock);

    return NNTI_OK;
}


NNTI_result_t NNTI_gni_atomic_fop (
		const NNTI_transport_t *trans_hdl,
		const NNTI_peer_t      *peer_hdl,
		const uint64_t          target_atomic,
		const uint64_t          result_atomic,
		const int64_t           operand,
		const NNTI_atomic_op_t  op,
		NNTI_work_request_t    *wr)
{
    int rc=0;
    trios_declare_timer(call_time);

    nnti_gni_work_request_t *gni_wr=NULL;

    log_level debug_level=nnti_debug_level;

    log_debug(nnti_ee_debug_level, "enter");

    if (config.use_wr_pool) {
        gni_wr=wr_pool_initiator_pop();
    } else {
        gni_wr=(nnti_gni_work_request_t *)malloc(sizeof(nnti_gni_work_request_t));
        memset(gni_wr, 0, sizeof(nnti_gni_work_request_t));
        nthread_lock_init(&gni_wr->lock);
    }
    assert(gni_wr);

    gni_wr->reg_buf = NULL;

    nnti_gni_connection_t *conn=get_conn_peer(peer_hdl);
    assert(conn);

    gni_wr->last_op  =GNI_OP_FETCH_ADD;
    gni_wr->peer_instance=peer_hdl->peer.NNTI_remote_process_t_u.gni.inst_id;
    log_debug(nnti_event_debug_level, "gni_wr->peer_instance==%lu", gni_wr->peer_instance);

    gni_wr->sge_count=1;
    gni_wr->sge_list =&gni_wr->sge;

    memset(&gni_wr->sge_list[0].post_desc, 0, sizeof(gni_post_descriptor_t));

    gni_wr->sge_list[0].post_desc.src_cq_hndl=transport_global_data.ep_cq_hdl;

    gni_wr->sge_list[0].post_desc.type   =GNI_POST_AMO;
    gni_wr->sge_list[0].post_desc.cq_mode=GNI_CQMODE_GLOBAL_EVENT;

    set_dlvr_mode(&gni_wr->sge_list[0].post_desc);
    set_rdma_mode(&gni_wr->sge_list[0].post_desc);

    gni_wr->sge_list[0].post_desc.local_addr     =(uint64_t)transport_global_data.atomics+(result_atomic*sizeof(int64_t));
    gni_wr->sge_list[0].post_desc.local_mem_hndl =transport_global_data.atomics_mem_hdl;
    gni_wr->sge_list[0].post_desc.remote_addr    =conn->atomics_addr+(target_atomic*sizeof(int64_t));
    gni_wr->sge_list[0].post_desc.remote_mem_hndl=conn->atomics_mem_hdl;
    gni_wr->sge_list[0].post_desc.length         =sizeof(int64_t);
    gni_wr->sge_list[0].post_desc.amo_cmd        =GNI_FMA_ATOMIC_FADD;
    gni_wr->sge_list[0].post_desc.first_operand  =operand;
    gni_wr->sge_list[0].post_desc.second_operand =0;

    gni_wr->sge_list[0].state=NNTI_GNI_SGE_STATE_STARTED;
    gni_wr->sge_list[0].gni_wr=gni_wr;

    gni_wr->nnti_wr = wr;

    wr->transport_id     =trans_hdl->id;
    wr->reg_buf          =(NNTI_buffer_t*)NULL;
    wr->ops              =NNTI_ATOMICS;
    wr->result           =NNTI_OK;
    wr->transport_private=(uint64_t)gni_wr;

    insert_sge_sgehash(&gni_wr->sge_list[0]);

    log_debug(debug_level, "calling PostFma(atomics fop - ep_hdl(%llu), local_addr=%llu, remote_addr=%llu)",
            conn->ep_hdl, gni_wr->sge_list[0].post_desc.local_addr, (uint64_t)gni_wr->sge_list[0].post_desc.remote_addr);

    nthread_lock(&nnti_gni_lock);
    GNI_EpSetEventData(
            conn->ep_hdl,
            hash6432shift((uint64_t)&gni_wr->sge_list[0]),
            0);
    rc=GNI_PostFma(conn->ep_hdl, &gni_wr->sge_list[0].post_desc);
    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "PostFma(atomics fop) failed: %d", rc);
    nthread_unlock(&nnti_gni_lock);

    return(NNTI_OK);
}


NNTI_result_t NNTI_gni_atomic_cswap (
		const NNTI_transport_t *trans_hdl,
		const NNTI_peer_t      *peer_hdl,
		const uint64_t          target_atomic,
		const uint64_t          result_atomic,
		const int64_t           compare_operand,
		const int64_t           swap_operand,
		NNTI_work_request_t    *wr)
{
    int rc=0;
    trios_declare_timer(call_time);

    nnti_gni_work_request_t *gni_wr=NULL;

    log_level debug_level=nnti_debug_level;

    log_debug(nnti_ee_debug_level, "enter");

    if (config.use_wr_pool) {
        gni_wr=wr_pool_initiator_pop();
    } else {
        gni_wr=(nnti_gni_work_request_t *)malloc(sizeof(nnti_gni_work_request_t));
        memset(gni_wr, 0, sizeof(nnti_gni_work_request_t));
        nthread_lock_init(&gni_wr->lock);
    }
    assert(gni_wr);

    gni_wr->reg_buf = NULL;

    nnti_gni_connection_t *conn=get_conn_peer(peer_hdl);
    assert(conn);

    gni_wr->last_op  =GNI_OP_FETCH_ADD;
    gni_wr->peer_instance=peer_hdl->peer.NNTI_remote_process_t_u.gni.inst_id;
    log_debug(nnti_event_debug_level, "gni_wr->peer_instance==%lu", gni_wr->peer_instance);

    gni_wr->sge_count=1;
    gni_wr->sge_list =&gni_wr->sge;

    memset(&gni_wr->sge_list[0].post_desc, 0, sizeof(gni_post_descriptor_t));

    gni_wr->sge_list[0].post_desc.src_cq_hndl=transport_global_data.ep_cq_hdl;

    gni_wr->sge_list[0].post_desc.type   =GNI_POST_AMO;
    gni_wr->sge_list[0].post_desc.cq_mode=GNI_CQMODE_GLOBAL_EVENT;

    set_dlvr_mode(&gni_wr->sge_list[0].post_desc);
    set_rdma_mode(&gni_wr->sge_list[0].post_desc);

    gni_wr->sge_list[0].post_desc.local_addr     =(uint64_t)transport_global_data.atomics+(result_atomic*sizeof(int64_t));
    gni_wr->sge_list[0].post_desc.local_mem_hndl =transport_global_data.atomics_mem_hdl;
    gni_wr->sge_list[0].post_desc.remote_addr    =conn->atomics_addr+(target_atomic*sizeof(int64_t));
    gni_wr->sge_list[0].post_desc.remote_mem_hndl=conn->atomics_mem_hdl;
    gni_wr->sge_list[0].post_desc.length         =sizeof(int64_t);
    gni_wr->sge_list[0].post_desc.amo_cmd        =GNI_FMA_ATOMIC_CSWAP;
    gni_wr->sge_list[0].post_desc.first_operand  =compare_operand;
    gni_wr->sge_list[0].post_desc.second_operand =swap_operand;

    gni_wr->sge_list[0].state=NNTI_GNI_SGE_STATE_STARTED;
    gni_wr->sge_list[0].gni_wr=gni_wr;

    gni_wr->nnti_wr = wr;

    wr->transport_id     =trans_hdl->id;
    wr->reg_buf          =(NNTI_buffer_t*)NULL;
    wr->ops              =NNTI_ATOMICS;
    wr->result           =NNTI_OK;
    wr->transport_private=(uint64_t)gni_wr;

    insert_sge_sgehash(&gni_wr->sge_list[0]);

    log_debug(debug_level, "calling PostFma(atomics cswap - ep_hdl(%llu), local_addr=%llu, remote_addr=%llu)",
            conn->ep_hdl, gni_wr->sge_list[0].post_desc.local_addr, (uint64_t)gni_wr->sge_list[0].post_desc.remote_addr);

    nthread_lock(&nnti_gni_lock);
    GNI_EpSetEventData(
            conn->ep_hdl,
            hash6432shift((uint64_t)&gni_wr->sge_list[0]),
            0);
    rc=GNI_PostFma(conn->ep_hdl, &gni_wr->sge_list[0].post_desc);
    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "PostFma(atomics cswap) failed: %d", rc);
    nthread_unlock(&nnti_gni_lock);

    return(NNTI_OK);
}


/**
 * @brief Create a receive work request that can be used to wait for buffer
 * operations to complete.
 *
 */
NNTI_result_t NNTI_gni_create_work_request (
        NNTI_buffer_t        *reg_buf,
        NNTI_work_request_t  *wr)
{
	nnti_gni_memory_handle_t *gni_mem_hdl=NULL;
	nnti_gni_work_request_t  *gni_wr=NULL;

	bool found=false;

	gni_mem_hdl=GNI_MEM_HDL(reg_buf);
	assert(gni_mem_hdl);

    log_debug(nnti_debug_level, "enter (reg_buf=%p ; wr=%p)", reg_buf, wr);

    if (reg_buf->ops == NNTI_RECV_QUEUE) {
    	nnti_gni_request_queue_handle_t *q_hdl=&transport_global_data.req_queue;
    	assert(q_hdl);

    	// head and tail are NOT adjacent.  search between them for a work request in the NNTI_RESET state.
    	if (q_hdl->tail <= q_hdl->head) {
    		// head has not wrapped around yet.  simple forward search.
    		for (int64_t index=q_hdl->tail; index<=q_hdl->head; index++) {
    			gni_wr=gni_mem_hdl->wr_queue->at(index);

    			log_debug(nnti_debug_level, "wr_queue[%lld]==%p", index, gni_wr);

    			log_debug(nnti_debug_level, "gni_wr(%p) gni_wr->state=%d", gni_wr, gni_wr->state);
    			if ((gni_wr->state==NNTI_GNI_WR_STATE_POSTED) || (gni_wr->state==NNTI_GNI_WR_STATE_RDMA_COMPLETE)) {
    				GNI_ATTACH_WR(wr,gni_wr);

    				if (index == q_hdl->head) {
    					q_hdl->head = (q_hdl->head+1) % q_hdl->req_count;
    					log_debug(nnti_debug_level, "new head==%lld", q_hdl->head);
    				}

    				found=true;

    				break;
    			}
    		}
    	} else {
    		// head has wrapped around.  search from tail to end, then front to head.
    		for (int64_t index=q_hdl->tail; index<q_hdl->req_count; index++) {
    			gni_wr=gni_mem_hdl->wr_queue->at(index);

    			log_debug(nnti_debug_level, "wr_queue[%lld]==%p", index, gni_wr);

    			log_debug(nnti_debug_level, "gni_wr(%p) gni_wr->state=%d", gni_wr, gni_wr->state);
    			if ((gni_wr->state==NNTI_GNI_WR_STATE_POSTED) || (gni_wr->state==NNTI_GNI_WR_STATE_RDMA_COMPLETE)) {
    				GNI_ATTACH_WR(wr,gni_wr);

    				found=true;

    				break;
    			}
    		}
    		if (!found) {
    			// didn't find anything at the end of the queue.  look at the front.
    			for (int64_t index=0; index<=q_hdl->head; index++) {
    				gni_wr=gni_mem_hdl->wr_queue->at(index);

    				log_debug(nnti_debug_level, "wr_queue[%lld]==%p", index, gni_wr);

        			log_debug(nnti_debug_level, "gni_wr(%p) gni_wr->state=%d", gni_wr, gni_wr->state);
        			if ((gni_wr->state==NNTI_GNI_WR_STATE_POSTED) || (gni_wr->state==NNTI_GNI_WR_STATE_RDMA_COMPLETE)) {
    					GNI_ATTACH_WR(wr,gni_wr);

    					if (index == q_hdl->head) {
    						q_hdl->head = (q_hdl->head+1) % q_hdl->req_count;
    						log_debug(nnti_debug_level, "new head==%lld", q_hdl->head);
    					}

    					found=true;

    					break;
    				}
    			}
    		}
    	}
    	log_debug(nnti_debug_level, "q->head=%lld ; q->tail==%lld", q_hdl->head, q_hdl->tail);
    } else if (reg_buf->ops == NNTI_RECV_DST) {
    	gni_wr=gni_mem_hdl->wr_queue->front();
		GNI_ATTACH_WR(wr,gni_wr);
    } else {

    }

    wr->transport_id     =reg_buf->transport_id;
    wr->reg_buf          =reg_buf;
    wr->ops              =reg_buf->ops;
    wr->result           =NNTI_OK;

    log_debug(nnti_debug_level, "assigned nnti_gni_work_request_t(%p) to NNTI_work_request(%p)", gni_wr, wr);

    log_debug(nnti_debug_level, "exit (reg_buf=%p ; wr=%p)", reg_buf, wr);

    return(NNTI_OK);
}


/**
 * @brief Disassociates a receive work request from a previous receive
 * and prepares it for reuse.
 *
 */
NNTI_result_t NNTI_gni_clear_work_request (
        NNTI_work_request_t  *wr)
{
    log_debug(nnti_debug_level, "enter (wr=%p)", wr);

    wr->result           =NNTI_OK;
    wr->transport_private=0;

    log_debug(nnti_debug_level, "exit (wr=%p)", wr);

    return(NNTI_OK);
}


/**
 * @brief Disassociates a receive work request from reg_buf.
 *
 */
NNTI_result_t NNTI_gni_destroy_work_request (
        NNTI_work_request_t  *wr)
{
	nnti_gni_memory_handle_t *gni_mem_hdl=NULL;
	nnti_gni_work_request_t  *gni_wr     =NULL;
	nnti_gni_work_request_t  *tail_wr    =NULL;

    wr_queue_iter_t q_victim;

	log_debug(nnti_debug_level, "enter (wr=%p)", wr);

	gni_mem_hdl=GNI_MEM_HDL(wr->reg_buf);
	assert(gni_mem_hdl);
	if ((gni_mem_hdl->type == SEND_BUFFER) || (gni_mem_hdl->type == GET_DST_BUFFER) || (gni_mem_hdl->type == PUT_SRC_BUFFER)) {
		// work requrests from initiator buffers are cleaned up in NNTI_gni_wait*().  nothing to do here.
		return(NNTI_OK);
	}

	gni_wr=GNI_WORK_REQUEST(wr);
	assert(gni_wr);

	if (wr->ops == NNTI_RECV_QUEUE) {
    	nnti_gni_request_queue_handle_t *q_hdl=&transport_global_data.req_queue;
    	assert(q_hdl);

		if (gni_wr->state == NNTI_GNI_WR_STATE_ATTACHED) {
			gni_wr->state=NNTI_GNI_WR_STATE_POSTED;
		} else {
			tail_wr=gni_mem_hdl->wr_queue->at(q_hdl->tail);
			assert(tail_wr);

			if (gni_wr == tail_wr) {
				if (q_hdl->tail <= q_hdl->head) {
					// head has not wrapped around yet.  move tail forward to the next work request NOT in the NNTI_COMPLETE state.
					for (int64_t index=q_hdl->tail; index<q_hdl->head; index++) {
						gni_wr=gni_mem_hdl->wr_queue->at(index);
						if (gni_wr->state == NNTI_GNI_WR_STATE_WAIT_COMPLETE) {
							gni_wr->state=NNTI_GNI_WR_STATE_POSTED;
						    gni_wr->last_op       =0;

							q_hdl->tail = (q_hdl->tail+1) % q_hdl->req_count;
						} else {
							break;
						}
					}
				} else {
					// head has wrapped around.  search from tail to end, then front to head.
					for (int64_t index=q_hdl->tail; index<q_hdl->req_count; index++) {
						gni_wr=gni_mem_hdl->wr_queue->at(index);
						if (gni_wr->state == NNTI_GNI_WR_STATE_WAIT_COMPLETE) {
							gni_wr->state=NNTI_GNI_WR_STATE_POSTED;
						    gni_wr->last_op       =0;

							q_hdl->tail = (q_hdl->tail+1) % q_hdl->req_count;
						} else {
							break;
						}
					}
					if (q_hdl->tail == 0) {
						// tail wrapped around in the previous loop.  keep going at the front.
						for (int64_t index=0; index<q_hdl->head; index++) {
							gni_wr=gni_mem_hdl->wr_queue->at(index);
							if (gni_wr->state == NNTI_GNI_WR_STATE_WAIT_COMPLETE) {
								gni_wr->state=NNTI_GNI_WR_STATE_POSTED;
							    gni_wr->last_op       =0;

								q_hdl->tail = (q_hdl->tail+1) % q_hdl->req_count;
							} else {
								break;
							}
						}
					}
				}
			}
		}
        if ((q_hdl->req_index >= q_hdl->req_count) &&
        	(q_hdl->head == 0)                 &&
        	(q_hdl->tail == 0)) {
        	reset_req_index(q_hdl);
        	log_debug(nnti_event_debug_level, "resetting req_index(%llu)", q_hdl->req_index);
        }

    	log_debug(nnti_debug_level, "q->head=%lld ; q->tail==%lld", q_hdl->head, q_hdl->tail);
    } else if (wr->ops == NNTI_RECV_DST) {
        nthread_lock(&gni_mem_hdl->wr_queue_lock);
		gni_wr->state=NNTI_GNI_WR_STATE_POSTED;
        q_victim=find(gni_mem_hdl->wr_queue->begin(), gni_mem_hdl->wr_queue->end(), gni_wr);
        if (q_victim != gni_mem_hdl->wr_queue->end()) {
            log_debug(nnti_debug_level, "erasing gni_wr=%p from the wr_queue", gni_wr);
            gni_mem_hdl->wr_queue->erase(q_victim);
        }
        repost_recv_work_request(wr->reg_buf, gni_wr);
        nthread_unlock(&gni_mem_hdl->wr_queue_lock);
    }

	GNI_DETACH_WR(wr, gni_wr);

    wr->transport_id=NNTI_TRANSPORT_NULL;
    wr->reg_buf     =NULL;
    wr->ops         =(NNTI_buf_ops_t)0;
    wr->result      =NNTI_OK;

    log_debug(nnti_debug_level, "exit (wr=%p)", wr);

    return(NNTI_OK);
}


/**
 * @brief Attempts to cancel an NNTI opertion.
 *
 */
NNTI_result_t NNTI_gni_cancel (
        NNTI_work_request_t *wr)
{
    log_debug(nnti_debug_level, "enter");

    log_debug(nnti_debug_level, "exit");

    return NNTI_ENOTSUP;
}


/**
 * @brief Attempts to cancel a list of NNTI opertions.
 *
 */
NNTI_result_t NNTI_gni_cancelall (
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count)
{
    log_debug(nnti_debug_level, "enter");

    log_debug(nnti_debug_level, "exit");

    return NNTI_ENOTSUP;
}


/**
 * @brief Interrupts NNTI_wait*()
 *
 */
NNTI_result_t NNTI_gni_interrupt (
        const NNTI_transport_t *trans_hdl)
{
	NNTI_result_t rc=NNTI_OK;
	int gni_rc=GNI_RC_SUCCESS;

    uint32_t  dummy=0xAAAAAAAA;
    uint32_t *ptr32=NULL;
	gni_post_descriptor_t  post_desc;
	gni_post_descriptor_t *post_desc_ptr;
    gni_cq_entry_t   ev_data;

    log_debug(nnti_debug_level, "enter");

    memset(&post_desc, 0, sizeof(gni_post_descriptor_t));
    post_desc.type           =GNI_POST_CQWRITE;
    post_desc.cq_mode        =GNI_CQMODE_GLOBAL_EVENT|GNI_CQMODE_REMOTE_EVENT;
    post_desc.dlvr_mode      =GNI_DLVMODE_IN_ORDER;
    post_desc.remote_mem_hndl=transport_global_data.interrupt_mem_hdl;

    ptr32=(uint32_t*)&post_desc.cqwrite_value;
    *ptr32=dummy;

    log_debug(nnti_debug_level, "calling PostCqWrite(cqwrite_value=%X)\n", *(uint32_t*)&post_desc.cqwrite_value);
    gni_rc=GNI_PostCqWrite(transport_global_data.interrupt_ep_hdl, &post_desc);
    if (gni_rc!=GNI_RC_SUCCESS) {
        log_error(nnti_debug_level, "PostCqWrite(post_desc) failed: %d\n", rc);
        rc=NNTI_EIO;
    }

    gni_rc=GNI_CqWaitEvent_wrapper(transport_global_data.interrupt_ep_cq_hdl, -1, &ev_data);
    DEQUEUE_POST_DESCRIPTOR(transport_global_data.interrupt_ep_cq_hdl, &ev_data, &post_desc_ptr);

    log_debug(nnti_debug_level, "exit");

    return rc;
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
NNTI_result_t NNTI_gni_wait (
        NNTI_work_request_t  *wr,
        const int             timeout,
        NNTI_status_t        *status)
{
    NNTI_result_t nnti_rc=NNTI_OK;

    nnti_gni_memory_handle_t *gni_mem_hdl=NULL;
    nnti_gni_work_request_t  *gni_wr=NULL;
    const NNTI_buffer_t      *wait_buf=NULL;

    wr_queue_iter_t q_victim;

    gni_cq_handle_t cq_hdl=0;
    gni_cq_entry_t  ev_data;

    gni_post_descriptor_t   *post_desc_ptr=NULL;

    NNTI_result_t rc=NNTI_OK;
    long elapsed_time = 0;
    long timeout_per_call = timeout;

    long entry_time=trios_get_time_ms();

    int8_t complete=FALSE;
    bool first_pass=true;

    trios_declare_timer(call_time);

    log_debug(nnti_ee_debug_level, "enter (wr=%p)", wr);

    assert(wr);
    assert(status);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "wr->reg_buf",
                "start of NNTI_gni_wait", wr->reg_buf);
    }

    gni_wr=GNI_WORK_REQUEST(wr);
    assert(gni_wr);

    check_listen_socket_for_new_connections();

    if ((config.max_timeout_ms > 0) && (config.max_timeout_ms < timeout)) {
        timeout_per_call=config.max_timeout_ms;
    }

    if (wr->result == NNTI_EDROPPED) {
    	log_debug(nnti_debug_level, "send was dropped (wr=%p ; gni_wr=%p)", wr, GNI_WORK_REQUEST(wr));
    	nnti_rc = NNTI_EDROPPED;
    } else if (is_wr_canceling(gni_wr)) {
    	log_debug(nnti_debug_level, "wr canceling (wr=%p ; gni_wr=%p)", wr, GNI_WORK_REQUEST(wr));
    	cancel_wr(gni_wr);
    } else if (is_wr_complete(gni_wr) == TRUE) {
    	log_debug(nnti_debug_level, "wr already complete (wr=%p ; gni_wr=%p)", wr, GNI_WORK_REQUEST(wr));
    	nnti_rc = NNTI_OK;
    } else {
        log_debug(nnti_debug_level, "wr NOT complete (wr=%p ; gni_wr=%p)", wr, GNI_WORK_REQUEST(wr));

        while (1) {
            rc=progress(timeout-elapsed_time);

            elapsed_time = (trios_get_time_ms() - entry_time);

            /* case 1: success */
            if (rc==NNTI_OK) {
//                logger_set_default_level(old_log_level);
                log_debug(nnti_debug_level, "progress() successful...");
                nnti_rc = rc;
            }
            /* case 2: timed out */
            else if (rc==NNTI_ETIMEDOUT) {
                log_debug(nnti_debug_level, "progress() timed out...");
//                logger_set_default_level(LOG_OFF);
            }
            /* case 3: interrupted */
            else if (rc==NNTI_EINTR) {
                log_debug(nnti_debug_level, "progress() interrupted...");
//                logger_set_default_level(LOG_OFF);
                nnti_rc = rc;
                break;
            }
            /* case 4: failure */
            else {
//                logger_set_default_level(old_log_level);
                log_debug(nnti_debug_level, "progress() failed: %s", strerror(errno));
                nnti_rc = rc;
                break;
            }

            if (is_wr_complete(gni_wr) == TRUE) {
                log_debug(nnti_debug_level, "wr completed (wr=%p ; gni_wr=%p)", wr, GNI_WORK_REQUEST(wr));
                nnti_rc = NNTI_OK;
                break;
            }

            /* if the caller asked for a legitimate timeout, we need to exit */
            if ((timeout >= 0) && (elapsed_time >= timeout)) {
                nnti_rc = NNTI_ETIMEDOUT;
                break;
            }

            if (trios_exit_now()) {
                log_debug(nnti_debug_level, "caught abort signal");
                nnti_rc = NNTI_ECANCELED;
                break;
            }
        }
    }

    create_status(wr, gni_wr, nnti_rc, &ev_data, status);

//    print_wc(&wc);
//    if (nnti_rc == NNTI_OK) {
//        print_raw_buf((char *)reg_buf->payload+gni_wr->wc.byte_offset, gni_mem_hdl->last_wc.byte_len);
//    }

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_status(logger_get_file(), "status",
                "end of NNTI_wait", status);
    }

    if (is_wr_complete(gni_wr)) {

        gni_wr->state=NNTI_GNI_WR_STATE_WAIT_COMPLETE;

        if (gni_wr->nnti_wr->ops == NNTI_ATOMICS) {
            del_sge_sgehash(&gni_wr->sge_list[0]);

            if (config.use_wr_pool) {
                wr_pool_initiator_push(gni_wr);
            } else {
                log_debug(nnti_debug_level, "free(gni_wr) (wr=%p, gni_wr=%p)", wr, gni_wr);
                free(gni_wr);
            }
            /*
             * This work request (wr) has reached a completed (final) state.  wr is reset here.
             */
            wr->transport_private=(uint64_t)NULL;

        } else {
            gni_mem_hdl=GNI_MEM_HDL(wr->reg_buf);
            assert(gni_mem_hdl);

            switch (gni_mem_hdl->type) {
                case REQUEST_BUFFER:
                case RECEIVE_BUFFER:
                    // defer cleanup to NNTI_gni_destroy_work_request()
                    break;
                case GET_SRC_BUFFER:
                case PUT_DST_BUFFER:
                case RDMA_TARGET_BUFFER:
                    if (config.use_rdma_target_ack) {
                        nthread_lock(&gni_mem_hdl->wr_queue_lock);
                        q_victim=find(gni_mem_hdl->wr_queue->begin(), gni_mem_hdl->wr_queue->end(), gni_wr);
                        if (q_victim != gni_mem_hdl->wr_queue->end()) {
                            log_debug(nnti_debug_level, "erasing gni_wr=%p from the wr_queue", gni_wr);
                            gni_mem_hdl->wr_queue->erase(q_victim);
                        }
                        repost_recv_work_request(wr->reg_buf, gni_wr);
                        nthread_unlock(&gni_mem_hdl->wr_queue_lock);
                    }
                    break;
                case SEND_BUFFER:
                case GET_DST_BUFFER:
                case PUT_SRC_BUFFER:
                    nthread_lock(&gni_mem_hdl->wr_queue_lock);
                    q_victim=find(gni_mem_hdl->wr_queue->begin(), gni_mem_hdl->wr_queue->end(), gni_wr);
                    if (q_victim != gni_mem_hdl->wr_queue->end()) {
                        log_debug(nnti_debug_level, "erasing gni_wr=%p from the wr_queue", gni_wr);
                        gni_mem_hdl->wr_queue->erase(q_victim);
                    }
                    nthread_unlock(&gni_mem_hdl->wr_queue_lock);

                    for (int j=0;j<gni_wr->sge_count;j++) {
                        del_sge_sgehash(&gni_wr->sge_list[j]);
                    }
                    if (gni_wr->sge_count > 1) {
                        free(gni_wr->sge_list);
                    }
                    del_wr_wrhash(gni_wr);
                    if (config.use_wr_pool) {
                        wr_pool_initiator_push(gni_wr);
                    } else {
                        log_debug(nnti_debug_level, "free(gni_wr) (wr=%p, gni_wr=%p)", wr, gni_wr);
                        free(gni_wr);
                    }
                    /*
                     * This work request (wr) has reached a completed (final) state.  wr is reset here.
                     */
                    wr->transport_private=(uint64_t)NULL;
                    break;
                case UNKNOWN_BUFFER:
                default:
                    log_error(nnti_debug_level, "unknown buffer type(%llu).", gni_mem_hdl->type);
                    break;
            }
        }
    }

    log_debug(nnti_ee_debug_level, "exit (wr=%p)", wr);
    return(nnti_rc);
}

/**
 * @brief Wait for <tt>remote_op</tt> on any buffer in <tt>buf_list</tt> to complete.
 *
 * Wait for <tt>remote_op</tt> on any buffer in <tt>buf_list</tt> to complete or timeout
 * waiting.  This is typically used to wait for a result or a bulk data
 * transfer.  The timeout is specified in milliseconds.  A timeout of <tt>-1</tt>
 * means wait forever.  A timeout of <tt>0</tt> means do not wait.
 *
 * Caveats:
 *   1) All buffers in buf_list must be registered with the same transport.
 *   2) You can't wait on the request queue and RDMA buffers in the same call.  Will probably be fixed in the future.
 */
NNTI_result_t NNTI_gni_waitany (
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count,
        const int             timeout,
        uint32_t             *which,
        NNTI_status_t        *status)
{
    NNTI_result_t nnti_rc=NNTI_OK;

    nnti_gni_memory_handle_t *gni_mem_hdl=NULL;
    const NNTI_buffer_t      *wait_buf=NULL;

    wr_queue_iter_t q_victim;

    nnti_gni_work_request_t *gni_wr=NULL;
//    int32_t                  which_wr=0;

    uint32_t which_cq=0;

    gni_cq_entry_t  ev_data;

    gni_post_descriptor_t   *post_desc_ptr=NULL;

    NNTI_result_t rc=NNTI_OK;
    long elapsed_time = 0;
    long timeout_per_call=timeout;

    int8_t complete=FALSE;
    bool first_pass=true;

    long entry_time=trios_get_time_ms();

    trios_declare_timer(call_time);

    log_debug(nnti_ee_debug_level, "enter");

    check_listen_socket_for_new_connections();

    assert(wr_list);
    assert(wr_count > 0);
    assert(status);

    if (wr_count == 1) {
        nnti_rc=NNTI_gni_wait(wr_list[0], timeout, status);
        *which=0;
        goto cleanup;
    }

    if ((config.max_timeout_ms > 0) && (config.max_timeout_ms < timeout)) {
        timeout_per_call=config.max_timeout_ms;
    }

    if (wr_count > 1) {
        for (uint32_t i=0;i<wr_count;i++) {
            if (wr_list[i] != NULL) {
                gni_wr=GNI_WORK_REQUEST(wr_list[i]);
                assert(gni_wr);
                if (is_wr_canceling(gni_wr)) {
                    cancel_wr(gni_wr);
                }
            }
        }
    }

    *which=0;

    if (is_any_wr_complete(wr_list, wr_count, which) == TRUE) {
        log_debug(nnti_debug_level, "any wr already complete (wr_list[%d]=%p ; gni_wr=%p)", *which, wr_list[*which], GNI_WORK_REQUEST(wr_list[*which]));
        nnti_rc = NNTI_OK;
    } else {
        log_debug(nnti_debug_level, "any wr NOT complete");

        while (1) {
            rc=progress(timeout-elapsed_time);

            elapsed_time = (trios_get_time_ms() - entry_time);

            /* case 1: success */
            if (rc==NNTI_OK) {
//                logger_set_default_level(old_log_level);
                log_debug(nnti_debug_level, "progress() successful...");
                nnti_rc = rc;
            }
            /* case 2: timed out */
            else if (rc==NNTI_ETIMEDOUT) {
                log_debug(nnti_debug_level, "progress() timed out...");
//                logger_set_default_level(LOG_OFF);
            }
            /* case 3: interrupted */
            else if (rc==NNTI_EINTR) {
                log_debug(nnti_debug_level, "progress() interrupted...");
//                logger_set_default_level(LOG_OFF);
                nnti_rc = rc;
                break;
            }
            /* case 4: failure */
            else {
//                logger_set_default_level(old_log_level);
                log_debug(nnti_debug_level, "progress() failed: %s", strerror(errno));
                nnti_rc = rc;
                break;
            }

            if (is_any_wr_complete(wr_list, wr_count, which) == TRUE) {
                log_debug(nnti_debug_level, "wr complete (wr_list[%d]=%p ; gni_wr=%p)", *which, wr_list[*which], GNI_WORK_REQUEST(wr_list[*which]));
                nnti_rc = NNTI_OK;
                break;
            }

            /* if the caller asked for a legitimate timeout, we need to exit */
            if ((timeout >= 0) && (elapsed_time >= timeout)) {
                nnti_rc = NNTI_ETIMEDOUT;
                break;
            }

            if (trios_exit_now()) {
                log_debug(nnti_debug_level, "caught abort signal");
                nnti_rc = NNTI_ECANCELED;
                break;
            }
        }
    }

    create_status(wr_list[*which], GNI_WORK_REQUEST(wr_list[*which]), nnti_rc, &ev_data, status);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_status(logger_get_file(), "status",
                "end of NNTI_waitany", status);
    }

    if (is_wr_complete(GNI_WORK_REQUEST(wr_list[*which]))) {

        GNI_WORK_REQUEST(wr_list[*which])->state=NNTI_GNI_WR_STATE_WAIT_COMPLETE;

        gni_mem_hdl=GNI_MEM_HDL(wr_list[*which]->reg_buf);
        assert(gni_mem_hdl);
        switch (gni_mem_hdl->type) {
            case REQUEST_BUFFER:
            case RECEIVE_BUFFER:
            	// defer cleanup to NNTI_gni_destroy_work_request()
                break;
            case GET_SRC_BUFFER:
            case PUT_DST_BUFFER:
            case RDMA_TARGET_BUFFER:
                if (config.use_rdma_target_ack) {
                    nthread_lock(&gni_mem_hdl->wr_queue_lock);
                    q_victim=find(gni_mem_hdl->wr_queue->begin(), gni_mem_hdl->wr_queue->end(), GNI_WORK_REQUEST(wr_list[*which]));
                    if (q_victim != gni_mem_hdl->wr_queue->end()) {
                        log_debug(nnti_debug_level, "erasing gni_wr=%p from the wr_queue", GNI_WORK_REQUEST(wr_list[*which]));
                        gni_mem_hdl->wr_queue->erase(q_victim);
                    }
                    repost_recv_work_request(wr_list[*which]->reg_buf, GNI_WORK_REQUEST(wr_list[*which]));
                    nthread_unlock(&gni_mem_hdl->wr_queue_lock);
                }
                break;
            case SEND_BUFFER:
            case GET_DST_BUFFER:
            case PUT_SRC_BUFFER:
                nthread_lock(&gni_mem_hdl->wr_queue_lock);
                q_victim=find(gni_mem_hdl->wr_queue->begin(), gni_mem_hdl->wr_queue->end(), GNI_WORK_REQUEST(wr_list[*which]));
                if (q_victim != gni_mem_hdl->wr_queue->end()) {
                    log_debug(nnti_debug_level, "erasing gni_wr=%p from the wr_queue", GNI_WORK_REQUEST(wr_list[*which]));
                    gni_mem_hdl->wr_queue->erase(q_victim);
                }
                nthread_unlock(&gni_mem_hdl->wr_queue_lock);

                for (int j=0;j<GNI_WORK_REQUEST(wr_list[*which])->sge_count;j++) {
                	del_sge_sgehash(&GNI_WORK_REQUEST(wr_list[*which])->sge_list[j]);
                }
                if (GNI_WORK_REQUEST(wr_list[*which])->sge_count > 1) {
                	free(GNI_WORK_REQUEST(wr_list[*which])->sge_list);
                }
                del_wr_wrhash(GNI_WORK_REQUEST(wr_list[*which]));
                if (config.use_wr_pool) {
                    wr_pool_initiator_push(GNI_WORK_REQUEST(wr_list[*which]));
                } else {
                    free(GNI_WORK_REQUEST(wr_list[*which]));
                }
                /*
                 * This work request (wr) has reached a completed (final) state.  wr is reset here.
                 */
                wr_list[*which]->transport_private=(uint64_t)NULL;
                break;
            case UNKNOWN_BUFFER:
            default:
                log_error(nnti_debug_level, "unknown buffer type(%llu).", gni_mem_hdl->type);
                break;
        }
    }

cleanup:
    log_debug(nnti_ee_debug_level, "exit");

    return(nnti_rc);
}

/**
 * @brief Wait for <tt>remote_op</tt> on all buffers in <tt>buf_list</tt> to complete.
 *
 * Wait for <tt>remote_op</tt> on all buffers in <tt>buf_list</tt> to complete or timeout
 * waiting.  This is typically used to wait for a result or a bulk data
 * transfer.  The timeout is specified in milliseconds.  A timeout of <tt>-1</tt>
 * means wait forever.  A timeout of <tt>0</tt> means do not wait.
 *
 * Caveats:
 *   1) All buffers in buf_list must be registered with the same transport.
 *   2) You can't wait on the receive queue and RDMA buffers in the same call.  Will probably be fixed in the future.
 */
NNTI_result_t NNTI_gni_waitall (
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count,
        const int             timeout,
        NNTI_status_t       **status)
{
    NNTI_result_t nnti_rc=NNTI_OK;

    nnti_gni_memory_handle_t *gni_mem_hdl=NULL;
    const NNTI_buffer_t      *wait_buf=NULL;

    wr_queue_iter_t q_victim;

    nnti_gni_work_request_t *gni_wr=NULL;
//    int32_t                  which_wr=0;

    uint32_t which_cq=0;

    gni_cq_entry_t  ev_data;

    gni_post_descriptor_t   *post_desc_ptr=NULL;

    NNTI_result_t rc=NNTI_OK;
    long elapsed_time = 0;
    long timeout_per_call=timeout;

    int8_t complete=FALSE;
    bool first_pass=true;

    long entry_time=trios_get_time_ms();

    trios_declare_timer(call_time);

    log_debug(nnti_ee_debug_level, "enter");

    check_listen_socket_for_new_connections();

    assert(wr_list);
    assert(wr_count > 0);
    assert(status);

    if (wr_count == 1) {
        nnti_rc=NNTI_gni_wait(wr_list[0], timeout, status[0]);
        goto cleanup;
    }

    if ((config.max_timeout_ms > 0) && (config.max_timeout_ms < timeout)) {
        timeout_per_call=config.max_timeout_ms;
    }

    if (wr_count > 1) {
        for (uint32_t i=0;i<wr_count;i++) {
            if (wr_list[i] != NULL) {
                gni_wr=GNI_WORK_REQUEST(wr_list[i]);
                assert(gni_wr);
                if (is_wr_canceling(gni_wr)) {
                    cancel_wr(gni_wr);
                }
            }
        }
    }

    if (is_all_wr_complete(wr_list, wr_count) == TRUE) {
        log_debug(nnti_debug_level, "all work requests already complete");
        nnti_rc = NNTI_OK;
    } else {
        log_debug(nnti_debug_level, "all work requests NOT complete");

        while (1) {
            rc=progress(timeout-elapsed_time);

            elapsed_time = (trios_get_time_ms() - entry_time);

            /* case 1: success */
            if (rc==NNTI_OK) {
//                logger_set_default_level(old_log_level);
                log_debug(nnti_debug_level, "progress() successful...");
                nnti_rc = rc;
            }
            /* case 2: timed out */
            else if (rc==NNTI_ETIMEDOUT) {
                log_debug(nnti_debug_level, "progress() timed out...");
//                logger_set_default_level(LOG_OFF);
            }
            /* case 3: interrupted */
            else if (rc==NNTI_EINTR) {
                log_debug(nnti_debug_level, "progress() interrupted...");
//                logger_set_default_level(LOG_OFF);
                nnti_rc = rc;
                break;
            }
            /* case 4: failure */
            else {
//                logger_set_default_level(old_log_level);
                log_debug(nnti_debug_level, "progress() failed: %s", strerror(errno));
                nnti_rc = rc;
                break;
            }

            if (is_all_wr_complete(wr_list, wr_count) == TRUE) {
                log_debug(nnti_debug_level, "all work requests complete");
                nnti_rc = NNTI_OK;
                break;
            }

            /* if the caller asked for a legitimate timeout, we need to exit */
            if ((timeout >= 0) && (elapsed_time >= timeout)) {
                nnti_rc = NNTI_ETIMEDOUT;
                break;
            }

            if (trios_exit_now()) {
                log_debug(nnti_debug_level, "caught abort signal");
                nnti_rc = NNTI_ECANCELED;
                break;
            }
        }
    }

    for (uint32_t i=0;i<wr_count;i++) {

    	if (wr_list[i] == NULL) {
    		continue;
    	}

        create_status(wr_list[i], GNI_WORK_REQUEST(wr_list[i]), nnti_rc, NULL, status[i]);

        if (logging_debug(nnti_debug_level)) {
            fprint_NNTI_status(logger_get_file(), "status",
                    "end of NNTI_waitall", status[i]);
        }

        if (is_wr_complete(GNI_WORK_REQUEST(wr_list[i]))) {

        	GNI_WORK_REQUEST(wr_list[i])->state=NNTI_GNI_WR_STATE_WAIT_COMPLETE;

            gni_mem_hdl=GNI_MEM_HDL(wr_list[i]->reg_buf);
            assert(gni_mem_hdl);
            switch (gni_mem_hdl->type) {
                case REQUEST_BUFFER:
                case RECEIVE_BUFFER:
                	// defer cleanup to NNTI_gni_destroy_work_request()
                    break;
                case GET_SRC_BUFFER:
                case PUT_DST_BUFFER:
                case RDMA_TARGET_BUFFER:
                    if (config.use_rdma_target_ack) {
                        nthread_lock(&gni_mem_hdl->wr_queue_lock);
                        q_victim=find(gni_mem_hdl->wr_queue->begin(), gni_mem_hdl->wr_queue->end(), GNI_WORK_REQUEST(wr_list[i]));
                        if (q_victim != gni_mem_hdl->wr_queue->end()) {
                            log_debug(nnti_debug_level, "erasing gni_wr=%p from the wr_queue", GNI_WORK_REQUEST(wr_list[i]));
                            gni_mem_hdl->wr_queue->erase(q_victim);
                        }
                        repost_recv_work_request(wr_list[i]->reg_buf, GNI_WORK_REQUEST(wr_list[i]));
                        nthread_unlock(&gni_mem_hdl->wr_queue_lock);
                    }
                    break;
                case SEND_BUFFER:
                case GET_DST_BUFFER:
                case PUT_SRC_BUFFER:
                    nthread_lock(&gni_mem_hdl->wr_queue_lock);
                    q_victim=find(gni_mem_hdl->wr_queue->begin(), gni_mem_hdl->wr_queue->end(), GNI_WORK_REQUEST(wr_list[i]));
                    if (q_victim != gni_mem_hdl->wr_queue->end()) {
                        log_debug(nnti_debug_level, "erasing gni_wr=%p from the wr_queue", GNI_WORK_REQUEST(wr_list[i]));
                        gni_mem_hdl->wr_queue->erase(q_victim);
                    }
                    nthread_unlock(&gni_mem_hdl->wr_queue_lock);

                    for (int j=0;j<GNI_WORK_REQUEST(wr_list[i])->sge_count;j++) {
                    	del_sge_sgehash(&GNI_WORK_REQUEST(wr_list[i])->sge_list[j]);
                    }
                    if (GNI_WORK_REQUEST(wr_list[i])->sge_count > 1) {
                    	free(GNI_WORK_REQUEST(wr_list[i])->sge_list);
                    }
                    del_wr_wrhash(GNI_WORK_REQUEST(wr_list[i]));
                    if (config.use_wr_pool) {
                        wr_pool_initiator_push(GNI_WORK_REQUEST(wr_list[i]));
                    } else {
                        free(GNI_WORK_REQUEST(wr_list[i]));
                    }
                    /*
                     * This work request (wr) has reached a completed (final) state.  wr is reset here.
                     */
                    wr_list[i]->transport_private=(uint64_t)NULL;
                    break;
                case UNKNOWN_BUFFER:
                default:
                    log_error(nnti_debug_level, "unknown buffer type(%llu).", gni_mem_hdl->type);
                    break;
            }
        }
    }


cleanup:
    log_debug(nnti_ee_debug_level, "exit");

    return(nnti_rc);
}

/**
 * @brief Disable this transport.
 *
 * Shutdown the transport.  Any outstanding sends, gets and puts will be
 * canceled.  Any new transport requests will fail.
 *
 */
NNTI_result_t NNTI_gni_fini (
        const NNTI_transport_t *trans_hdl)
{
    int rc=GNI_RC_SUCCESS; /* return code */

    log_debug(nnti_ee_debug_level, "enter");

    close_all_conn();

    if (config.use_wr_pool) {
        rc=wr_pool_fini();
        if (rc!=NNTI_OK) {
            log_error(nnti_debug_level, "wr_pool_fini() failed: %d", rc);
            rc=NNTI_EINVAL;
        }
    }

    close(transport_global_data.listen_sock);
    transport_global_data.listen_name[0]='\0';
    transport_global_data.listen_addr=0;
    transport_global_data.listen_port=0;

    rc=GNI_CdmDestroy(transport_global_data.cdm_hdl);
    if (rc!=GNI_RC_SUCCESS) {
        log_error(nnti_debug_level, "CdmCreate() failed: %d", rc);
        rc=NNTI_EINVAL;
    }

    nthread_lock_fini(&nnti_gni_lock);
    nthread_lock_fini(&nnti_mem_lock);

    nthread_lock_fini(&nnti_progress_lock);
    nthread_cond_fini(&nnti_progress_cond);

    nthread_lock_fini(&nnti_conn_peer_lock);
    nthread_lock_fini(&nnti_conn_instance_lock);
    nthread_lock_fini(&nnti_wr_wrhash_lock);
    nthread_lock_fini(&nnti_buf_bufhash_lock);
    nthread_lock_fini(&nnti_wr_pool_lock);

    nthread_lock_fini(&transport_global_data.atomics_lock);

    gni_initialized = false;

    log_debug(nnti_ee_debug_level, "exit");

    return(NNTI_OK);
}

static NNTI_result_t setup_atomics(void)
{
    NNTI_result_t rc=NNTI_OK; /* return code */
    int gni_rc=GNI_RC_SUCCESS; /* return code */

    trios_declare_timer(call_time);

    uint32_t atomics_bytes;

    log_debug(nnti_debug_level, "enter");

    atomics_bytes=config.min_atomics_vars * sizeof(int64_t);
    trios_start_timer(call_time);
    transport_global_data.atomics=(int64_t*)aligned_malloc(atomics_bytes);
    if (transport_global_data.atomics == NULL) {
        return(NNTI_ENOMEM);
    }
    memset(transport_global_data.atomics, 0, atomics_bytes);
    trios_stop_timer("malloc and memset", call_time);

    trios_start_timer(call_time);
    gni_rc=GNI_MemRegister (transport_global_data.nic_hdl,
            (uint64_t)transport_global_data.atomics,
            atomics_bytes,
            NULL,
            GNI_MEM_READWRITE,
            (uint32_t)-1,
            &transport_global_data.atomics_mem_hdl);
    if (gni_rc!=GNI_RC_SUCCESS) {
        log_error(nnti_debug_level, "MemRegister(atomics_mem_hdl) failed: gni_rc=%d, %s", gni_rc, strerror(errno));
        rc=NNTI_EIO;
    }
    trios_stop_timer("atomics register", call_time);

    log_debug(nnti_debug_level, "exit (atomics==%p)...", transport_global_data.atomics);

    return(rc);
}

static void *aligned_malloc(
        size_t size)
{
    int rc = 0;
    void * ptr = NULL;

#ifdef NO_MEMALIGN
      ptr = malloc( size );
#else
      rc = posix_memalign( &ptr , 128 , size );
#endif

    if (rc!=0) {
        log_error(nnti_debug_level, "posix_memalign() failed: rc=%d (%s)", rc, strerror(errno));
        ptr=NULL;
    }

    return ptr;
}

static gni_mem_handle_t register_memory_segment(nnti_gni_buffer_type_t type, void *buf, uint64_t len, uint64_t extra)
{
	NNTI_result_t nnti_rc=NNTI_OK;
    int gni_rc=GNI_RC_SUCCESS; /* return code */

    uint32_t flags=GNI_MEM_READWRITE;

    trios_declare_timer(call_time);

    gni_mem_handle_t mem_hdl;
    gni_cq_handle_t mem_cq_hdl;

    mem_cq_hdl=NULL;

    log_debug(nnti_debug_level, "enter type(%d) buffer(%p) len(%llu) extra(%llu)", type, buf, len, extra);

    if (config.pi_ordering==PI_ORDERING_STRICT) {
        log_debug(nnti_debug_level, "using STRICT ordering");
        flags = GNI_MEM_READWRITE|GNI_MEM_STRICT_PI_ORDERING;
    } else if (config.pi_ordering==PI_ORDERING_RELAXED) {
        log_debug(nnti_debug_level, "using RELAXED ordering");
        flags = GNI_MEM_READWRITE|GNI_MEM_RELAXED_PI_ORDERING;
    } else if (config.pi_ordering==PI_ORDERING_DEFAULT) {
        log_debug(nnti_debug_level, "using DEFAULT ordering");
        flags = GNI_MEM_READWRITE;
    }

    if (need_mem_cq(type) == 1) {
        mem_cq_hdl=transport_global_data.mem_cq_hdl;
    }

    trios_start_timer(call_time);
    gni_rc=GNI_MemRegister (transport_global_data.nic_hdl,
            (uint64_t)buf,
            len+extra,
            mem_cq_hdl,
            flags,
            (uint32_t)-1,
            &mem_hdl);
    if (gni_rc!=GNI_RC_SUCCESS) {
        log_error(nnti_debug_level, "MemRegister(mem_hdl) failed: gni_rc=%d, %s", gni_rc, strerror(errno));
        goto cleanup;
    }
    trios_stop_timer("buf register", call_time);

    log_debug(nnti_debug_level, "register mem_cq_hdl  =%llu", (uint64_t)mem_cq_hdl);
    log_debug(nnti_debug_level, "register hdl->mem_hdl=(%llu,%llu)", (uint64_t)mem_hdl.qword1, (uint64_t)mem_hdl.qword2);

cleanup:
    switch(gni_rc) {
        case GNI_RC_SUCCESS:
            nnti_rc=NNTI_OK;
            break;
        default:
            nnti_rc=NNTI_EIO;
            break;
    }

    log_debug(nnti_debug_level, "exit  type(%d) buffer(%p) gni_rc(%d) nnti_rc(%d)", type, buf, gni_rc, nnti_rc);

    return (mem_hdl);
}

static NNTI_result_t register_memory(
        const NNTI_transport_t *trans_hdl,
        char                   *buffer,
        const uint64_t          element_size,
        const uint64_t          extra,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf)
{
	NNTI_result_t rc=NNTI_OK;
    trios_declare_timer(call_time);

    NNTI_buffer_t     *old_buf=NULL;
    nnti_gni_memory_handle_t *gni_mem_hdl=NULL;

    assert(trans_hdl);
    assert(buffer);
    assert(element_size>0);
    assert(extra>=0);
    assert(num_elements>0);
    assert(ops>0);
    assert(reg_buf);

    log_debug(nnti_ee_debug_level, "enter (reg_buf=%p)", reg_buf);

nthread_lock(&nnti_mem_lock);

    old_buf=get_buf_bufhash(hash6432shift((uint64_t)buffer));
    if (old_buf==NULL) {
        gni_mem_hdl=(nnti_gni_memory_handle_t*)calloc(1, sizeof(nnti_gni_memory_handle_t));
        assert(gni_mem_hdl);
        gni_mem_hdl->wr_queue =new wr_queue_t;
        gni_mem_hdl->ref_count=1;
        nthread_lock_init(&gni_mem_hdl->wr_queue_lock);
    } else {
        gni_mem_hdl=(nnti_gni_memory_handle_t*)old_buf->transport_private;
        gni_mem_hdl->ref_count++;
    }

    log_debug(nnti_ee_debug_level, "gni_mem_hdl->ref_count==%lu", gni_mem_hdl->ref_count);

    reg_buf->transport_id      = trans_hdl->id;
    reg_buf->buffer_owner      = trans_hdl->me;
    reg_buf->ops               = ops;
    reg_buf->payload_size      = element_size;
    reg_buf->payload           = (uint64_t)buffer;
    reg_buf->transport_private = (uint64_t)gni_mem_hdl;

    log_debug(nnti_debug_level, "rpc_buffer->payload_size=%ld",
            reg_buf->payload_size);

    if (gni_mem_hdl->ref_count==1) {

    	gni_mem_hdl->extra=extra;

        if (ops == NNTI_RECV_QUEUE) {
            nnti_gni_request_queue_handle_t *q_hdl=&transport_global_data.req_queue;

            gni_mem_hdl->type   =REQUEST_BUFFER;

            memset(q_hdl, 0, sizeof(nnti_gni_request_queue_handle_t));

            q_hdl->reg_buf=reg_buf;

            server_req_queue_init(
                    q_hdl,
                    buffer,
                    element_size,
                    extra,
                    num_elements);

            reg_buf->payload_size=q_hdl->req_size;

            post_recv_queue_work_requests(reg_buf);

        } else if (ops == NNTI_RECV_DST) {
            gni_mem_hdl->type    =RECEIVE_BUFFER;
        } else if (ops == NNTI_SEND_SRC) {
            gni_mem_hdl->type=SEND_BUFFER;
        } else if (ops == NNTI_GET_DST) {
            gni_mem_hdl->type    =GET_DST_BUFFER;
        } else if (ops == NNTI_GET_SRC) {
            gni_mem_hdl->type    =GET_SRC_BUFFER;
        } else if (ops == NNTI_PUT_SRC) {
            gni_mem_hdl->type    =PUT_SRC_BUFFER;
        } else if (ops == NNTI_PUT_DST) {
            gni_mem_hdl->type    =PUT_DST_BUFFER;
        } else if (ops == (NNTI_GET_SRC|NNTI_PUT_DST)) {
            gni_mem_hdl->type    =RDMA_TARGET_BUFFER;
        } else {
            gni_mem_hdl->type=UNKNOWN_BUFFER;
        }

        if (ops != NNTI_RECV_QUEUE) {
        	gni_mem_hdl->mem_hdl_list=&gni_mem_hdl->mem_hdl;
        	gni_mem_hdl->mem_hdl_count=1;

        	gni_mem_hdl->mem_hdl=register_memory_segment(gni_mem_hdl->type, buffer, element_size, extra);
        }
    }

    if (ops == NNTI_RECV_DST) {
        post_recv_work_request(reg_buf);
    }

    if (config.use_rdma_target_ack) {
        if ((gni_mem_hdl->type == RDMA_TARGET_BUFFER) ||
                (gni_mem_hdl->type == GET_SRC_BUFFER) ||
                (gni_mem_hdl->type == PUT_DST_BUFFER)) {
            post_recv_work_request(reg_buf);
        }
    }

    if (rc==NNTI_OK) {
        reg_buf->buffer_segments.NNTI_remote_addr_array_t_val=(NNTI_remote_addr_t *)calloc(1, sizeof(NNTI_remote_addr_t));
        reg_buf->buffer_segments.NNTI_remote_addr_array_t_len=1;
        assert(reg_buf->buffer_segments.NNTI_remote_addr_array_t_val);

        reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].transport_id                            = NNTI_TRANSPORT_GEMINI;
        reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.mem_hdl.qword1 = gni_mem_hdl->mem_hdl.qword1;
        reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.mem_hdl.qword2 = gni_mem_hdl->mem_hdl.qword2;

        if (gni_mem_hdl->type==REQUEST_BUFFER) {
            reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.size = transport_global_data.req_queue.req_size;
            reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.buf  = (uint64_t)transport_global_data.req_queue.req_buffer;
            reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.type = NNTI_GNI_REQUEST_BUFFER;
        } else {
            reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.size = reg_buf->payload_size;
            reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.buf  = (uint64_t)reg_buf->payload;
            reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.type = NNTI_GNI_SEND_SRC;
        }

        if (gni_mem_hdl->ref_count==1) {
            insert_buf_bufhash(reg_buf);
            log_debug(nnti_debug_level, "gni_mem_hdl->type==%llu",
                    (uint64_t)gni_mem_hdl->type);
            log_debug(nnti_debug_level, "reg_buf.buf.hash==%llu",
                    (uint64_t)hash6432shift(reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.buf));
            log_debug(nnti_debug_level, "ref_count==1 called insert_buf_bufhash() (reg_buf=%p, reg_buf.hash6432=%llu)",
                    reg_buf, (uint64_t)hash6432shift(reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.buf));
            NNTI_buffer_t *tmp_buf=get_buf_bufhash((uint64_t)hash6432shift(reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.buf));
            log_debug(nnti_debug_level, "immediate get_buf_bufhash() says tmp_buf=%p", tmp_buf);
        }
    }

nthread_unlock(&nnti_mem_lock);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "reg_buf",
                "end of NNTI_gni_register_memory", reg_buf);
    }

    log_debug(nnti_ee_debug_level, "exit (reg_buf=%p, reg_buf.hash6432=%llu)", reg_buf, (uint64_t)hash6432shift(reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.buf));

    return(rc);
}

static NNTI_result_t unregister_memory(gni_mem_handle_t mem_hdl)
{
    int rc=GNI_RC_SUCCESS; /* return code */
    trios_declare_timer(call_time);

    log_debug(nnti_debug_level, "enter mem_hdl(%p)", mem_hdl);

    log_debug(nnti_debug_level, "unregister mem_hdl=(%llu,%llu)", (uint64_t)mem_hdl.qword1, (uint64_t)mem_hdl.qword2);

    trios_start_timer(call_time);
    rc=GNI_MemDeregister (transport_global_data.nic_hdl, &mem_hdl);
    if (rc!=GNI_RC_SUCCESS) {
        log_error(nnti_debug_level, "MemDeregister(mem_hdl) failed: %d", rc);
    }
    trios_stop_timer("buf deregister", call_time);

    switch(rc) {
        case GNI_RC_SUCCESS:
            rc=NNTI_OK;
        default:
            rc=NNTI_EIO;
    }

    log_debug(nnti_debug_level, "exit mem_hdl(%p)", mem_hdl);

    return ((NNTI_result_t)rc);
}

static int need_mem_cq(nnti_gni_buffer_type_t type)
{
    int need_cq=0;

    log_debug(nnti_ee_debug_level, "enter");

    if (!config.use_rdma_events) {
        switch (type) {
            case RDMA_TARGET_BUFFER:
            case GET_SRC_BUFFER:
            case PUT_DST_BUFFER:
            case GET_DST_BUFFER:
            case PUT_SRC_BUFFER:
            case SEND_BUFFER:
                need_cq=0;
                break;
            case RECEIVE_BUFFER:
            case REQUEST_BUFFER:
                need_cq=1;
                break;
            case UNKNOWN_BUFFER:
            default:
                need_cq=0;
                break;
        }
    } else {
        switch (type) {
        case GET_DST_BUFFER:
        case PUT_SRC_BUFFER:
        case SEND_BUFFER:
            need_cq=0;
            break;
        case GET_SRC_BUFFER:
            if (config.rdma_mode==RDMA_CROSSOVER) {
                need_cq=0;
            } else {
                need_cq=1;
            }
            break;
        case PUT_DST_BUFFER:
        case RECEIVE_BUFFER:
        case REQUEST_BUFFER:
            need_cq=1;
            break;
        case RDMA_TARGET_BUFFER:
            if (config.rdma_mode==RDMA_CROSSOVER) {
                need_cq=0;
            } else {
                need_cq=1;
            }
            break;
        case UNKNOWN_BUFFER:
        default:
            need_cq=0;
            break;
        }
    }

    log_debug(nnti_ee_debug_level, "exit(%d)", need_cq);

    return(need_cq);
}

static int need_wc_mem_cq(nnti_gni_memory_handle_t *gni_mem_hdl)
{
    int need_cq=0;

    log_debug(nnti_ee_debug_level, "enter");

    assert(gni_mem_hdl);

    if (!config.use_rdma_target_ack) {
        switch (gni_mem_hdl->type) {
            case RDMA_TARGET_BUFFER:
            case GET_SRC_BUFFER:
            case PUT_DST_BUFFER:
            case GET_DST_BUFFER:
            case PUT_SRC_BUFFER:
            case SEND_BUFFER:
                need_cq=0;
                break;
            case RECEIVE_BUFFER:
            case REQUEST_BUFFER:
                need_cq=1;
                break;
            case UNKNOWN_BUFFER:
            default:
                need_cq=0;
                break;
        }
    } else {
        switch (gni_mem_hdl->type) {
            case GET_DST_BUFFER:
            case PUT_SRC_BUFFER:
            case SEND_BUFFER:
                need_cq=0;
                break;
            case RDMA_TARGET_BUFFER:
            case GET_SRC_BUFFER:
            case PUT_DST_BUFFER:
            case RECEIVE_BUFFER:
            case REQUEST_BUFFER:
                need_cq=1;
                break;
            case UNKNOWN_BUFFER:
            default:
                need_cq=0;
                break;
        }
    }

    log_debug(nnti_ee_debug_level, "exit(%d)", need_cq);

    return(need_cq);
}

static nnti_gni_sge_t *decode_sge(
        gni_cq_entry_t      *ev_data)
{
    const NNTI_buffer_t     *event_buf=NULL;
    nnti_gni_sge_t          *gni_sge  =NULL;
    nnti_gni_work_request_t *gni_wr   =NULL;

    log_debug(nnti_debug_level, "enter");

    if (gni_cq_get_source(*ev_data) == 6) {
    	// this is a target buffer.  determine if it is REQUEST_BUFFER or RECEIVE_BUFFER.
    	event_buf = get_buf_bufhash((uint32_t)gni_cq_get_inst_id(*ev_data));
    	if (event_buf == NULL) {
    		// ev_data.inst_id is not a buffer hash.  assume this is a REQUEST_BUFFER.
            log_debug(nnti_debug_level, "ev_data.source==6 and ev_data.inst_id is NOT a hash, "
            		"so the event buffer is a REQUEST BUFFER and ev_data.inst_id is the request buffer index.");

            uint64_t index=(uint64_t)gni_cq_get_inst_id(*ev_data);
            nnti_gni_memory_handle_t *gni_mem_hdl=GNI_MEM_HDL(transport_global_data.req_queue.reg_buf);

            gni_wr =gni_mem_hdl->wr_queue->at(index);
            gni_sge=&gni_wr->sge_list[0];
    	} else {
            log_debug(nnti_debug_level, "ev_data.source==6 and ev_data.inst_id is a hash, so the event buffer is a RECEIVE BUFFER.");

            nnti_gni_memory_handle_t *gni_mem_hdl=GNI_MEM_HDL(event_buf);

            gni_wr =gni_mem_hdl->wr_queue->front();
            gni_sge=&gni_wr->sge_list[0];
    	}
    } else {
    	gni_sge=get_sge_sgehash((uint32_t)gni_cq_get_inst_id(*ev_data));
        if (gni_sge==NULL) {
			event_buf = get_buf_bufhash((uint32_t)gni_cq_get_inst_id(*ev_data));
			if (event_buf == NULL) {
                log_error(nnti_debug_level, "ev_data.inst_id does NOT match either a bufhash or a sgehash.  what happened?  where did this event come from?");
                print_cq_event(ev_data, true);

                assert(gni_sge);
            } else {
                nnti_gni_memory_handle_t *gni_mem_hdl=GNI_MEM_HDL(event_buf);

                gni_wr =gni_mem_hdl->wr_queue->front();
                gni_sge=&gni_wr->sge_list[0];
            }
        }
        assert(gni_sge);
    }

    log_debug(nnti_debug_level, "exit (gni_sge==%p)", gni_sge);

    return(gni_sge);
}

static nnti_gni_work_request_t *decode_work_request(
        gni_cq_entry_t      *ev_data)
{
    const NNTI_buffer_t *event_buf=NULL;
    nnti_gni_work_request_t  *gni_wr=NULL;

    log_debug(nnti_debug_level, "enter");

    if (gni_cq_get_source(*ev_data) == 6) {
    	// this is a target buffer.  determine if it is REQUEST_BUFFER or RECEIVE_BUFFER.
    	event_buf = get_buf_bufhash((uint32_t)gni_cq_get_inst_id(*ev_data));
    	if (event_buf == NULL) {
    		// ev_data.inst_id is not a buffer hash.  assume this is a REQUEST_BUFFER.
            log_debug(nnti_debug_level, "ev_data.source==6 and ev_data.inst_id is NOT a hash, "
            		"so the event buffer is a REQUEST BUFFER and ev_data.inst_id is the request buffer index.");

            uint64_t index=(uint64_t)gni_cq_get_inst_id(*ev_data);
            nnti_gni_memory_handle_t *gni_mem_hdl=GNI_MEM_HDL(transport_global_data.req_queue.reg_buf);

            gni_wr=gni_mem_hdl->wr_queue->at(index);
    	} else {
            log_debug(nnti_debug_level, "ev_data.source==6 and ev_data.inst_id is a hash, so the event buffer is a RECEIVE BUFFER.");

            nnti_gni_memory_handle_t *gni_mem_hdl=GNI_MEM_HDL(event_buf);

            gni_wr=gni_mem_hdl->wr_queue->front();
    	}
    } else {
    	nnti_gni_sge_t *sge=get_sge_sgehash((uint32_t)gni_cq_get_inst_id(*ev_data));
        if (sge==NULL) {
			event_buf = get_buf_bufhash((uint32_t)gni_cq_get_inst_id(*ev_data));
			if (event_buf == NULL) {
                log_error(nnti_debug_level, "ev_data.inst_id does NOT match either a bufhash or a wrhash.  what happened?  where did this event come from?");
                print_cq_event(ev_data, true);

                assert(gni_wr);
            } else {
                nnti_gni_memory_handle_t *gni_mem_hdl=GNI_MEM_HDL(event_buf);

                gni_wr=gni_mem_hdl->wr_queue->front();
            }
        } else {
        	gni_wr=sge->gni_wr;
        }
        assert(gni_wr);
    }

    log_debug(nnti_debug_level, "exit (gni_wr==%p)", gni_wr);

    return(gni_wr);
}

static int cancel_wr(
        nnti_gni_work_request_t *gni_wr)
{
    NNTI_result_t rc=NNTI_OK;

    log_debug(nnti_debug_level, "enter (gni_wr=%p)", gni_wr);

    gni_wr->state=NNTI_GNI_WR_STATE_WAIT_COMPLETE;
    gni_wr->nnti_wr->result=NNTI_ECANCELED;

    log_debug(nnti_debug_level, "exit (gni_wr==%p)", gni_wr);

    return(rc);
}

static int process_event(
        const NNTI_buffer_t     *reg_buf,
        nnti_gni_sge_t          *gni_sge,
        gni_cq_handle_t          cq_hdl,
        gni_cq_entry_t          *ev_data,
        gni_post_descriptor_t   *post_desc_ptr)
{
    int rc=NNTI_OK;
    nnti_gni_memory_handle_t *gni_mem_hdl=NULL;
    nnti_gni_work_request_t  *gni_wr=NULL;
    const NNTI_buffer_t      *event_buf=NULL;

    log_level debug_level=nnti_debug_level;

    log_debug(debug_level, "enter (reg_buf=%p, gni_sge=%p, cq_hdl=%llu, ev_data=%llu)", reg_buf, gni_sge, (uint64_t)cq_hdl, (uint64_t)ev_data);

    if (!GNI_CQ_STATUS_OK(*ev_data)) {
        return NNTI_EIO;
    }

    assert(gni_sge);

    gni_wr=gni_sge->gni_wr;
    assert(gni_wr);

    if ((gni_wr->nnti_wr) && (gni_wr->nnti_wr->ops == NNTI_ATOMICS)) {
        gni_wr->state=NNTI_GNI_WR_STATE_RDMA_COMPLETE;
        gni_wr->nnti_wr->result=NNTI_OK;
        return NNTI_OK;
    }

    log_debug(debug_level, "ev_data.inst_id==%llu, reg_buf.buf.hash==%llu",
            (uint64_t)gni_cq_get_inst_id(*ev_data),
            (uint64_t)hash6432shift(reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.buf));

    event_buf=reg_buf;
    gni_mem_hdl=(nnti_gni_memory_handle_t *)event_buf->transport_private;
    assert(gni_mem_hdl);

    log_debug(debug_level, "event_buf=%p; gni_sge=%p; gni_wr=%p; gni_wr->last_op=%lld; gni_wr.hash==%llu",
    		event_buf, gni_sge, gni_wr, (int64_t)gni_wr->last_op, (uint64_t)hash6432shift((uint64_t)gni_wr));
    log_debug(debug_level, "gni_mem_hdl->type==%llu", (uint64_t)gni_mem_hdl->type);

    debug_level=nnti_debug_level;
    switch (gni_mem_hdl->type) {
        case SEND_BUFFER:
            if ((gni_wr->last_op==GNI_OP_SEND_REQUEST) &&
                (gni_sge->state==NNTI_GNI_SGE_STATE_STARTED)) {

                log_debug(debug_level, "SEND request event - event_buf==%p, state==%d", event_buf, gni_sge->state);

                if (post_desc_ptr == &gni_sge->post_desc) {
                    log_debug(debug_level, "SEND request completion - event_buf==%p", event_buf);
                    gni_sge->state=NNTI_GNI_SGE_STATE_COMPLETE;
                } else {
                    log_debug(debug_level, "SEND request - unknown post descriptor - post_desc_ptr(%p) != &gni_sge->post_desc(%p)",
                            post_desc_ptr, &gni_sge->post_desc);
                }
                gni_wr->wc->byte_offset=gni_wr->wc->src_offset;
            } else if ((gni_wr->last_op==GNI_OP_SEND_BUFFER) &&
                    (gni_sge->state==NNTI_GNI_SGE_STATE_STARTED)) {

                log_debug(debug_level, "SEND buffer event - event_buf==%p, state==%d", event_buf, gni_sge->state);

                if (post_desc_ptr == &gni_sge->post_desc) {
                    log_debug(debug_level, "SEND buffer completion - event_buf==%p", event_buf);
                    gni_sge->state=NNTI_GNI_SGE_STATE_COMPLETE;
                } else {
                    log_debug(debug_level, "SEND buffer - unknown post descriptor - post_desc_ptr(%p) != &gni_sge->post_desc(%p)",
                            post_desc_ptr, &gni_sge->post_desc);
                }
                gni_wr->wc->byte_offset=gni_wr->wc->src_offset;
            }
            break;
        case PUT_SRC_BUFFER:
            gni_wr->last_op=GNI_OP_PUT_INITIATOR;
            if (gni_sge->state==NNTI_GNI_SGE_STATE_STARTED) {
                log_debug(debug_level, "RDMA write event - event_buf==%p, gni_sge->state==%d", event_buf, gni_sge->state);
                if (post_desc_ptr == &gni_sge->post_desc) {
                    log_debug(debug_level, "RDMA write (initiator) completion - event_buf==%p", event_buf);
                    gni_sge->state=NNTI_GNI_SGE_STATE_COMPLETE;
                } else {
                    log_debug(debug_level, "RDMA write - unknown post descriptor - post_desc_ptr(%p) != &gni_wr->post_desc(%p)",
                            post_desc_ptr, &gni_sge->post_desc);
                }
                gni_wr->wc->byte_offset=gni_wr->wc->src_offset;
            }
            break;
        case GET_DST_BUFFER:
            gni_wr->last_op=GNI_OP_GET_INITIATOR;
            if (gni_sge->state==NNTI_GNI_SGE_STATE_STARTED) {
                log_debug(debug_level, "RDMA read event - event_buf==%p, state==%d", event_buf, gni_sge->state);
                if (post_desc_ptr == &gni_sge->post_desc) {
                    log_debug(debug_level, "RDMA read (initiator) completion - event_buf==%p", event_buf);
                    gni_sge->state=NNTI_GNI_SGE_STATE_COMPLETE;
                } else {
                    log_debug(debug_level, "RDMA read - unknown post descriptor - post_desc_ptr(%p) != &gni_wr->post_desc(%p)",
                            post_desc_ptr, &gni_sge->post_desc);
                }
                gni_wr->wc->byte_offset=gni_wr->wc->dest_offset;
            }
            break;
        case REQUEST_BUFFER:
            {
            uint64_t  index=0;
            nnti_gni_request_queue_handle_t *q=&transport_global_data.req_queue;

            gni_wr->last_op=GNI_OP_NEW_REQUEST;

            index = (uint64_t)gni_cq_get_inst_id(*ev_data);
            log_debug(debug_level, "wc_index(%llu)", index);

            nnti_gni_work_request_t *tmp_wr=gni_mem_hdl->wr_queue->at(index);

            assert(tmp_wr == gni_wr);

            gni_wr->wc->ack_received=1;

            log_debug(debug_level, "recv completion - event_buf=%p index=%llu", event_buf, index);

            gni_wr->state=NNTI_GNI_WR_STATE_RDMA_COMPLETE;

//            if ((q->req_index >= q->req_count) &&
//           	    (q->tail == 0)) {
//                reset_req_index(q);
//                log_debug(nnti_event_debug_level, "resetting req_index(%llu)", q->req_index);
//            }
            }
            break;
        case RECEIVE_BUFFER:
            gni_wr->last_op=GNI_OP_PUT_TARGET;
            log_debug(debug_level, "RDMA write event - event_buf==%p, state==%d", event_buf, gni_wr->state);
            if (config.use_rdma_events) {
//                if ((gni_wr->op_state.rdma_init    ==true)   &&
//                    (gni_wr->op_state.rdma_complete==false)) {
//                    log_debug(debug_level, "RDMA write (receive buffer) completion - event_buf==%p", event_buf);
//                    gni_wr->op_state.rdma_complete=true;
//                } else if ((gni_wr->op_state.rdma_init    ==true)   &&
//                           (gni_wr->op_state.rdma_complete==true)   &&
//                           (gni_wr->op_state.wc_complete ==false)) {
//                    log_debug(debug_level, "RDMA write ACK (receive buffer) completion - event_buf==%p", event_buf);
//                    gni_wr->op_state.wc_complete=true;
//                    gni_wr->wc->byte_offset=gni_wr->wc->dest_offset;
//                }
            } else {
                log_debug(debug_level, "RDMA write ACK (receive buffer) completion - event_buf==%p", event_buf);
                gni_wr->state=NNTI_GNI_WR_STATE_RDMA_COMPLETE;
                gni_wr->wc->byte_offset=gni_wr->wc->dest_offset;
            }
            break;
        case PUT_DST_BUFFER:
            gni_wr->last_op=GNI_OP_PUT_TARGET;
            log_debug(debug_level, "RDMA write event - event_buf==%p, state==%d", event_buf, gni_wr->state);
            if (config.use_rdma_events) {
//                if ((gni_wr->op_state.rdma_init    ==true)   &&
//                    (gni_wr->op_state.rdma_complete==false)) {
//                    log_debug(debug_level, "RDMA write (target) completion - event_buf==%p", event_buf);
//                    if (config.use_rdma_target_ack) {
//                        gni_wr->op_state.rdma_complete=true;
//                    } else {
//                        gni_wr->op_state.rdma_complete=true;
//                        gni_wr->wc->byte_offset=gni_wr->wc->dest_offset;
//                    }
//                } else if ((gni_wr->op_state.rdma_init    ==true)   &&
//                           (gni_wr->op_state.rdma_complete==true)   &&
//                           (gni_wr->op_state.wc_complete ==false)) {
//                    log_debug(debug_level, "RDMA write ACK (target) completion - event_buf==%p", event_buf);
//                    gni_wr->op_state.wc_complete=true;
//                    gni_wr->wc->byte_offset=gni_wr->wc->dest_offset;
//                }
            } else {
                log_debug(debug_level, "RDMA write ACK (target) completion - event_buf==%p", event_buf);
                gni_wr->state=NNTI_GNI_WR_STATE_RDMA_COMPLETE;
                gni_wr->wc->byte_offset=gni_wr->wc->dest_offset;
            }
            break;
        case GET_SRC_BUFFER:
            gni_wr->last_op=GNI_OP_GET_TARGET;
            log_debug(debug_level, "RDMA read event - event_buf==%p, state==%d", event_buf, gni_wr->state);
            if (config.use_rdma_events) {
//                if ((gni_wr->op_state.rdma_init    ==true)   &&
//                    (gni_wr->op_state.rdma_complete==false)) {
//                    log_debug(debug_level, "RDMA read (target) completion - event_buf==%p", event_buf);
//                    if (config.use_rdma_target_ack) {
//                        gni_wr->op_state.rdma_complete=true;
//                    } else {
//                        gni_wr->op_state.rdma_complete=true;
//                        gni_wr->wc->byte_offset=gni_wr->wc->src_offset;
//                    }
//                } else if ((gni_wr->op_state.rdma_init    ==true)   &&
//                           (gni_wr->op_state.rdma_complete==true)   &&
//                           (gni_wr->op_state.wc_complete ==false)) {
//                    log_debug(debug_level, "RDMA read ACK (target) completion - event_buf==%p", event_buf);
//                    gni_wr->op_state.wc_complete=true;
//                    gni_wr->wc->byte_offset=gni_wr->wc->src_offset;
//                }
            } else {
                log_debug(debug_level, "RDMA read ACK (target) completion - event_buf==%p", event_buf);
                gni_wr->state=NNTI_GNI_WR_STATE_RDMA_COMPLETE;
                gni_wr->wc->byte_offset=gni_wr->wc->src_offset;
            }
            break;
        case RDMA_TARGET_BUFFER:
            log_debug(debug_level, "RDMA target event - event_buf=%p, state=%d, last_op=%d",
                    event_buf, gni_wr->state, gni_wr->last_op);
            if ((gni_wr->last_op==GNI_OP_GET_INITIATOR) ||
                (gni_wr->last_op==GNI_OP_PUT_INITIATOR)) {

                if (gni_sge->state==NNTI_GNI_SGE_STATE_STARTED) {
                    if (post_desc_ptr == &gni_sge->post_desc) {
                        log_debug(debug_level, "RDMA target (initiator) completion - event_buf==%p", event_buf);
                        gni_sge->state=NNTI_GNI_SGE_STATE_COMPLETE;
//                    } else if (post_desc_ptr == &gni_wr->wc_post_desc) {
//                        log_debug(debug_level, "RDMA target ACK (initiator) completion - event_buf==%p", event_buf);
//                        gni_wr->op_state.wc_complete=true;
                    } else {
//                        log_debug(debug_level, "RDMA target (initiator) - unknown post descriptor - post_desc_ptr(%p) != &gni_wr->post_desc(%p) != &gni_wr->wc_post_desc(%p)",
//                                post_desc_ptr, &gni_wr->post_desc, &gni_wr->wc_post_desc);
                        log_debug(debug_level, "RDMA target (initiator) - unknown post descriptor - post_desc_ptr(%p) != &gni_sge->post_desc(%p)",
                                post_desc_ptr, &gni_sge->post_desc);
                    }
                    gni_wr->wc->byte_offset=gni_wr->wc->dest_offset;
                }

            } else {
                if (config.use_rdma_events) {
//                    if ((gni_wr->op_state.rdma_init==true)       &&
//                        (gni_wr->op_state.rdma_complete==false)) {
//                        log_debug(debug_level, "RDMA target (target) completion - event_buf==%p", event_buf);
//                        if (config.use_rdma_target_ack) {
//                            gni_wr->op_state.rdma_complete=true;
//                        } else {
//                            gni_wr->op_state.rdma_complete=true;
//                            gni_wr->wc->byte_offset=gni_wr->wc->src_offset;
//                        }
//                    } else if ((gni_wr->op_state.rdma_init==true)     &&
//                               (gni_wr->op_state.rdma_complete==true) &&
//                               (gni_wr->op_state.wc_complete==false)) {
//                        log_debug(debug_level, "RDMA target ACK (target) completion - event_buf==%p", event_buf);
//                        gni_wr->op_state.wc_complete=true;
//                        gni_wr->wc->byte_offset=gni_wr->wc->src_offset;
//                    }
                } else {
                    log_debug(debug_level, "RDMA target ACK (target) completion - event_buf==%p", event_buf);
                    gni_wr->state=NNTI_GNI_WR_STATE_RDMA_COMPLETE;
                    gni_wr->wc->byte_offset=gni_wr->wc->src_offset;
                }
            }
            break;
        case UNKNOWN_BUFFER:
        default:
            break;
    }

    print_wc(gni_wr->wc);

    log_debug(nnti_ee_debug_level, "exit");
    return (rc);
}

static NNTI_result_t post_recv_queue_work_requests(
        NNTI_buffer_t *reg_buf)
{
    nnti_gni_work_request_t  *gni_wr=NULL;
    nnti_gni_memory_handle_t *gni_mem_hdl=NULL;

    nnti_gni_request_queue_handle_t *q=&transport_global_data.req_queue;

    log_debug(nnti_ee_debug_level, "enter (reg_buf=%p)", reg_buf);

    gni_mem_hdl=(nnti_gni_memory_handle_t *)reg_buf->transport_private;
    assert(gni_mem_hdl);

    for (uint64_t index=0;index < q->req_count;index++) {
    	gni_wr=(nnti_gni_work_request_t *)malloc(sizeof(nnti_gni_work_request_t));
    	memset(gni_wr, 0, sizeof(nnti_gni_work_request_t));
    	nthread_lock_init(&gni_wr->lock);
    	assert(gni_wr);

    	gni_wr->sge_list=&gni_wr->sge;
    	gni_wr->sge_count=1;
    	gni_wr->sge_list[0].gni_wr=gni_wr;

    	gni_wr->reg_buf=reg_buf;
    	gni_wr->state  =NNTI_GNI_WR_STATE_POSTED;

    	gni_wr->wc = GNI_WC_ADDRESS(reg_buf, index);

    	gni_mem_hdl->wr_queue->push_back(gni_wr);
    }

    log_debug(nnti_ee_debug_level, "exit (reg_buf=%p)", reg_buf);

    return(NNTI_OK);
}

static NNTI_result_t post_recv_work_request(
        NNTI_buffer_t *reg_buf)
{
    nnti_gni_work_request_t  *gni_wr=NULL;
    nnti_gni_memory_handle_t *gni_mem_hdl=NULL;

    log_debug(nnti_ee_debug_level, "enter (reg_buf=%p)", reg_buf);

    gni_mem_hdl=(nnti_gni_memory_handle_t *)reg_buf->transport_private;
    assert(gni_mem_hdl);

    if (config.use_wr_pool) {
        if (need_wc_mem_cq(gni_mem_hdl) == 1) {
            gni_wr=wr_pool_target_pop();
        } else {
            gni_wr=wr_pool_initiator_pop();
        }
    } else {
        gni_wr=(nnti_gni_work_request_t *)malloc(sizeof(nnti_gni_work_request_t));
        memset(gni_wr, 0, sizeof(nnti_gni_work_request_t));
        nthread_lock_init(&gni_wr->lock);
    }
    assert(gni_wr);

	gni_wr->sge_list=&gni_wr->sge;
	gni_wr->sge_count=1;
	gni_wr->sge_list[0].gni_wr=gni_wr;

    gni_wr->reg_buf=reg_buf;
	gni_wr->state  =NNTI_GNI_WR_STATE_POSTED;

	gni_wr->wc = GNI_WC_ADDRESS(reg_buf, 0);

    nthread_lock(&gni_mem_hdl->wr_queue_lock);
    gni_mem_hdl->wr_queue->push_back(gni_wr);
    nthread_unlock(&gni_mem_hdl->wr_queue_lock);

    log_debug(nnti_ee_debug_level, "exit (reg_buf=%p)", reg_buf);

    return(NNTI_OK);
}

static NNTI_result_t repost_recv_work_request(
        NNTI_buffer_t    *reg_buf,
        nnti_gni_work_request_t *gni_wr)
{
    nnti_gni_memory_handle_t *gni_mem_hdl=NULL;

    log_debug(nnti_ee_debug_level, "enter (reg_buf=%p)", reg_buf);

    gni_mem_hdl=(nnti_gni_memory_handle_t *)reg_buf->transport_private;
    assert(gni_mem_hdl);

    gni_wr->last_op=0;

	gni_wr->state=NNTI_GNI_WR_STATE_POSTED;

    gni_mem_hdl->wr_queue->push_back(gni_wr);

    log_debug(nnti_ee_debug_level, "exit (reg_buf=%p)", reg_buf);

    return(NNTI_OK);
}

static int8_t is_wr_canceling(
        nnti_gni_work_request_t *gni_wr)
{
    int8_t rc=FALSE;

    if (gni_wr->state==NNTI_GNI_WR_STATE_CANCELING) {
        rc=TRUE;
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);

    return(rc);
}

static int8_t is_wr_complete(
        nnti_gni_work_request_t *gni_wr)
{
    int8_t rc=FALSE;

    if ((gni_wr->state==NNTI_GNI_WR_STATE_RDMA_COMPLETE) ||
    	(gni_wr->state==NNTI_GNI_WR_STATE_WAIT_COMPLETE)) {
    	rc=TRUE;
    } else {
    	rc=TRUE;  /* assume complete.  flip to incomplete if any SGE is not complete. */
    	for (int i=0;i<gni_wr->sge_count;i++) {
    		if (gni_wr->sge_list[i].state != NNTI_GNI_SGE_STATE_COMPLETE) {
    			rc=FALSE;
    			break;
    		}
    	}
    	if (rc==TRUE) {
    		// all SGEs are complete, so the WR is too.
    		gni_wr->state=NNTI_GNI_WR_STATE_RDMA_COMPLETE;
    	}
    }

    log_debug(nnti_ee_debug_level, "exit (rc=%d)", rc);
    return(rc);
}

static int8_t is_wr_complete(
        NNTI_work_request_t *wr)
{
    nnti_gni_work_request_t *gni_wr=NULL;

    log_debug(nnti_debug_level, "enter (wr=%p)", wr);

    gni_wr=GNI_WORK_REQUEST(wr);
    assert(gni_wr);

    return(is_wr_complete(gni_wr));
}

static int8_t is_any_wr_complete(
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count,
        uint32_t             *which)
{
    int8_t rc=FALSE;

    log_debug(nnti_debug_level, "enter");

    for (uint32_t i=0;i<wr_count;i++) {
        if ((wr_list[i] != NULL) &&
            (is_wr_complete(wr_list[i]) == TRUE)) {

            *which=i;
            rc = TRUE;
            break;
        }
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);

    return(rc);
}

static int8_t is_all_wr_complete(
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count)
{
    int8_t rc=TRUE;

    log_debug(nnti_debug_level, "enter");

    for (uint32_t i=0;i<wr_count;i++) {
        if ((wr_list[i] != NULL) &&
            (is_wr_complete(wr_list[i]) == FALSE)) {

            rc = FALSE;
            break;
        }
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);

    return(rc);
}

//static void create_wr_array(
//        const NNTI_buffer_t **buf_list,
//        const uint32_t        buf_count,
//        nnti_gni_work_request_t   ***wr_list,
//        uint32_t             *wr_count)
//{
//    *wr_list=(nnti_gni_work_request_t **)calloc(sizeof(nnti_gni_work_request_t *), buf_count);
//    assert(*wr_list);
//    *wr_count=buf_count;
//    for(uint32_t i=0;i<buf_count;i++) {
//        if (buf_list[i]==NULL) {
//            (*wr_list)[i]=NULL;
//            continue;
//        }
//        nnti_gni_memory_handle_t *gni_mem_hdl=(nnti_gni_memory_handle_t *)buf_list[i]->transport_private;
//        assert(gni_mem_hdl);
//
//        nthread_lock(&gni_mem_hdl->wr_queue_lock);
//        if (gni_mem_hdl->wr_queue->empty()) {
//            (*wr_list)[i]=NULL;
//        } else {
//            (*wr_list)[i]=gni_mem_hdl->wr_queue->front();
//            assert((*wr_list)[i]);
//        }
//        nthread_unlock(&gni_mem_hdl->wr_queue_lock);
//    }
//
//    return;
//}


//static nnti_gni_work_request_t *first_incomplete_wr(
//        nnti_gni_memory_handle_t *gni_mem_hdl)
//{
//    nnti_gni_work_request_t  *wr=NULL;
//
//    log_debug(nnti_ee_debug_level, "enter");
//
//    assert(gni_mem_hdl);
//
////    nthread_lock(&gni_mem_hdl->wr_queue_lock);
//
//    if (gni_mem_hdl->wr_queue->empty()) {
//        log_debug(nnti_ee_debug_level, "work request queue is empty");
//    } else {
//        wr_queue_iter_t i;
//        for (i=gni_mem_hdl->wr_queue->begin(); i != gni_mem_hdl->wr_queue->end(); i++) {
//            wr=*i;
//            assert(wr);
//            if (is_wr_complete(wr) == FALSE) {
//                log_debug(nnti_ee_debug_level, "incomplete work request found");
//                break;
//            }
//        }
//        if (wr==NULL) {
//            log_debug(nnti_ee_debug_level, "all work requests complete");
//        }
//    }
//
////    nthread_unlock(&gni_mem_hdl->wr_queue_lock);
//
//    log_debug(nnti_ee_debug_level, "exit (wr=%p)", wr);
//    return(wr);
//}

//static int8_t is_wr_queue_empty(
//        const NNTI_buffer_t *reg_buf)
//{
//    int8_t rc=FALSE;
//    nnti_gni_memory_handle_t *gni_mem_hdl=NULL;
//
//    log_debug(nnti_ee_debug_level, "enter");
//
//    gni_mem_hdl=(nnti_gni_memory_handle_t *)reg_buf->transport_private;
//    assert(gni_mem_hdl);
//
//    if (gni_mem_hdl->wr_queue->empty()) {
//        rc=TRUE;
//    }
//
//    log_debug(nnti_ee_debug_level, "exit (rc=%d)", rc);
//    return(rc);
//}

//static int8_t is_buf_op_complete(
//        const NNTI_buffer_t *reg_buf)
//{
//    int8_t rc=FALSE;
//    nnti_gni_memory_handle_t *gni_mem_hdl=NULL;
//    nnti_gni_work_request_t  *wr=NULL;
//
//    log_debug(nnti_ee_debug_level, "enter");
//
//    gni_mem_hdl=(nnti_gni_memory_handle_t *)reg_buf->transport_private;
//    assert(gni_mem_hdl);
//
//    nthread_lock(&gni_mem_hdl->wr_queue_lock);
//
//    if (is_wr_queue_empty(reg_buf) == TRUE) {
//        log_debug(nnti_ee_debug_level, "work request queue is empty - return FALSE");
//        rc=FALSE;
//    } else {
//        wr=gni_mem_hdl->wr_queue->front();
//        assert(wr);
//
//        rc = is_wr_complete(wr);
//    }
//
//    nthread_unlock(&gni_mem_hdl->wr_queue_lock);
//
//    log_debug(nnti_ee_debug_level, "exit (rc=%d)", rc);
//    return(rc);
//}
//
//static int8_t is_any_buf_op_complete(
//        const NNTI_buffer_t **buf_list,
//        const uint32_t        buf_count,
//        uint32_t             *which)
//{
//    int8_t rc=FALSE;
//
//    log_debug(nnti_debug_level, "enter");
//
//    for (int i=0;i<buf_count;i++) {
//        if ((buf_list[i] != NULL) &&
//            (is_wr_queue_empty(buf_list[i]) == FALSE) &&
//            (is_buf_op_complete(buf_list[i]) == TRUE)) {
//
//            *which=i;
//            rc = TRUE;
//            break;
//        }
//    }
//
//    log_debug(nnti_debug_level, "exit (rc=%d)", rc);
//
//    return(rc);
//}
//
//static int8_t is_all_buf_ops_complete(
//        const NNTI_buffer_t **buf_list,
//        const uint32_t        buf_count)
//{
//    int8_t rc=TRUE;
//
//    log_debug(nnti_debug_level, "enter");
//
//    for (int i=0;i<buf_count;i++) {
//        if ((buf_list[i] != NULL) &&
//            (is_wr_queue_empty(buf_list[i]) == FALSE) &&
//            (is_buf_op_complete(buf_list[i]) == FALSE)) {
//
//            rc = FALSE;
//            break;
//        }
//    }
//
//    log_debug(nnti_debug_level, "exit (rc=%d)", rc);
//
//    return(rc);
//}

static void create_status(
        NNTI_work_request_t     *wr,
        nnti_gni_work_request_t *gni_wr,
        NNTI_result_t            nnti_rc,
        gni_cq_entry_t          *ev_data,
        NNTI_status_t           *status)
{
    nnti_gni_connection_t    *conn       =NULL;

    memset(status, 0, sizeof(NNTI_status_t));
    status->op = wr->ops;
    if (is_wr_complete(gni_wr)) {
        status->result = wr->result;
    } else {
        status->result = nnti_rc;
    }
    if (status->result==NNTI_OK) {
        if (gni_wr->peer_instance) {
        	conn = get_conn_instance(gni_wr->peer_instance);
        } else {
        	conn = get_conn_instance(gni_wr->wc->inst_id);
        }
        assert(conn);

        if (gni_wr->nnti_wr->ops != NNTI_ATOMICS) {
            status->start  = (uint64_t)gni_wr->reg_buf->payload;
            status->offset = gni_wr->wc->byte_offset;
            status->length = gni_wr->wc->byte_len;
        }
        switch (gni_wr->last_op) {
            case GNI_OP_PUT_INITIATOR:
            case GNI_OP_GET_TARGET:
            case GNI_OP_SEND_REQUEST:
            case GNI_OP_SEND_BUFFER:
                create_peer(&status->src,
                        transport_global_data.listen_name,
                        transport_global_data.listen_addr,
                        transport_global_data.listen_port,
                        transport_global_data.alps_info.ptag,
                        transport_global_data.alps_info.cookie,
                        transport_global_data.instance);
                create_peer(&status->dest,
                        conn->peer_name,
                        conn->peer_addr,
                        conn->peer_port,
                        conn->peer_ptag,
                        conn->peer_cookie,
                        conn->peer_instance);
                break;
            case GNI_OP_GET_INITIATOR:
            case GNI_OP_PUT_TARGET:
            case GNI_OP_NEW_REQUEST:
            case GNI_OP_RECEIVE:
                create_peer(&status->src,
                        conn->peer_name,
                        conn->peer_addr,
                        conn->peer_port,
                        conn->peer_ptag,
                        conn->peer_cookie,
                        conn->peer_instance);
                create_peer(&status->dest,
                        transport_global_data.listen_name,
                        transport_global_data.listen_addr,
                        transport_global_data.listen_port,
                        transport_global_data.alps_info.ptag,
                        transport_global_data.alps_info.cookie,
                        transport_global_data.instance);
                break;
            }
    }
}

static void create_peer(NNTI_peer_t *peer, char *name, NNTI_ip_addr addr, NNTI_tcp_port port, uint32_t ptag, uint32_t cookie, NNTI_instance_id instance)
{
    log_debug(nnti_ee_debug_level, "enter");

    sprintf(peer->url, "gni://%s:%u/?ptag=%lu&cookie=%lu", name, ntohs(port), (uint64_t)ptag, (uint64_t)cookie);

    peer->peer.transport_id                       =NNTI_TRANSPORT_GEMINI;
    peer->peer.NNTI_remote_process_t_u.gni.addr   =addr;
    peer->peer.NNTI_remote_process_t_u.gni.port   =port;
    peer->peer.NNTI_remote_process_t_u.gni.inst_id=instance;

    log_debug(nnti_ee_debug_level, "exit");
}

static NNTI_result_t insert_conn_peer(const NNTI_peer_t *peer, nnti_gni_connection_t *conn)
{
    NNTI_result_t  rc=NNTI_OK;
    addrport_key key;

    key.addr = peer->peer.NNTI_remote_process_t_u.gni.addr;
    key.port = peer->peer.NNTI_remote_process_t_u.gni.port;

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_peer(logger_get_file(), "peer",
                "insert_conn_peer", peer);
    }

    if (nthread_lock(&nnti_conn_peer_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    if (connections_by_peer.find(key) != connections_by_peer.end()) {
        print_peer_map();
        assert(connections_by_peer.find(key) == connections_by_peer.end());
    }
    connections_by_peer[key] = conn;   // add to connection map
    nthread_unlock(&nnti_conn_peer_lock);

    log_debug(nnti_debug_level, "peer connection added (conn=%p)", conn);

    return(rc);
}
static NNTI_result_t insert_conn_instance(const NNTI_instance_id instance, nnti_gni_connection_t *conn)
{
    NNTI_result_t  rc=NNTI_OK;

    if (nthread_lock(&nnti_conn_instance_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    if (connections_by_instance.find(instance) != connections_by_instance.end()) {
        print_instance_map();
        assert(connections_by_instance.find(instance) == connections_by_instance.end());
    }
    connections_by_instance[instance] = conn;
    nthread_unlock(&nnti_conn_instance_lock);

    log_debug(nnti_debug_level, "instance connection added (conn=%p)", conn);

    return(rc);
}
static nnti_gni_connection_t *get_conn_peer(const NNTI_peer_t *peer)
{
    nnti_gni_connection_t *conn = NULL;
    addrport_key   key;

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_peer(logger_get_file(), "looking for peer",
                "get_conn_peer", peer);
    }

    memset(&key, 0, sizeof(addrport_key));
    key.addr=peer->peer.NNTI_remote_process_t_u.gni.addr;
    key.port=peer->peer.NNTI_remote_process_t_u.gni.port;

    if (nthread_lock(&nnti_conn_peer_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    if (connections_by_peer.find(key) != connections_by_peer.end()) {
        conn = connections_by_peer[key];
    }
    nthread_unlock(&nnti_conn_peer_lock);

    if (conn != NULL) {
        log_debug(nnti_debug_level, "connection found");
        return conn;
    }

    log_debug(nnti_debug_level, "connection NOT found");
    print_peer_map();

    return(NULL);
}
static nnti_gni_connection_t *get_conn_instance(const NNTI_instance_id instance)
{
    nnti_gni_connection_t *conn=NULL;

    log_debug(nnti_debug_level, "looking for instance=%llu", (unsigned long long)instance);
    if (nthread_lock(&nnti_conn_instance_lock)) log_warn(nnti_debug_level, "failed to get thread lock");;
    if (connections_by_instance.find(instance) != connections_by_instance.end()) {
        conn = connections_by_instance[instance];
    }
    nthread_unlock(&nnti_conn_instance_lock);

    if (conn != NULL) {
        log_debug(nnti_debug_level, "connection found");
        return conn;
    }

    log_debug(nnti_debug_level, "connection NOT found (instance==%llu)", (uint64_t)instance);
    print_instance_map();

    return(NULL);
}
static NNTI_peer_t *get_peer_by_url(const char *url)
{
    nnti_gni_connection_t *conn = NULL;

    log_debug(nnti_debug_level, "looking for url=%s", url);

    conn_by_peer_iter_t i;
    if (nthread_lock(&nnti_conn_peer_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    for (i=connections_by_peer.begin(); i != connections_by_peer.end(); i++) {
        log_debug(nnti_debug_level, "peer_map key=(%llu,%llu) conn=%p, url=%s",
                (uint64_t)i->first.addr, (uint64_t)i->first.port, i->second, i->second->peer.url);
        if (strcmp(i->second->peer.url, url) == 0) {
            conn=i->second;
            break;
        }
    }
    nthread_unlock(&nnti_conn_peer_lock);

    if (conn != NULL) {
        log_debug(nnti_debug_level, "peer found");
        NNTI_peer_t *peer=(NNTI_peer_t*)malloc(sizeof(NNTI_peer_t));
        *peer=conn->peer;
        return peer;
    }

    return(NULL);
}
static nnti_gni_connection_t *del_conn_peer(const NNTI_peer_t *peer)
{
    nnti_gni_connection_t *conn=NULL;
    addrport_key    key;

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_peer(logger_get_file(), "peer",
                "get_conn_peer", peer);
    }

    memset(&key, 0, sizeof(addrport_key));
    key.addr=peer->peer.NNTI_remote_process_t_u.gni.addr;
    key.port=peer->peer.NNTI_remote_process_t_u.gni.port;

    if (nthread_lock(&nnti_conn_peer_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    if (connections_by_peer.find(key) != connections_by_peer.end()) {
        conn = connections_by_peer[key];
    }

    if (conn != NULL) {
        log_debug(nnti_debug_level, "connection found");
        connections_by_peer.erase(key);
        del_conn_instance(conn->peer_instance);
    } else {
        log_debug(nnti_debug_level, "connection NOT found");
    }
    nthread_unlock(&nnti_conn_peer_lock);

    return(conn);
}
static nnti_gni_connection_t *del_conn_instance(const NNTI_instance_id instance)
{
    nnti_gni_connection_t *conn=NULL;
    log_level debug_level = nnti_debug_level;

    if (nthread_lock(&nnti_conn_instance_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    if (connections_by_instance.find(instance) != connections_by_instance.end()) {
        conn = connections_by_instance[instance];
    }

    if (conn != NULL) {
        log_debug(debug_level, "connection found");
        connections_by_instance.erase(instance);
    } else {
        log_debug(debug_level, "connection NOT found");
    }
    nthread_unlock(&nnti_conn_instance_lock);

    return(conn);
}
static void print_peer_map()
{
    if (!logging_debug(nnti_debug_level)) {
        return;
    }

    conn_by_peer_iter_t i;
    for (i=connections_by_peer.begin(); i != connections_by_peer.end(); i++) {
        log_debug(nnti_debug_level, "peer_map key=(%llu,%llu) conn=%p",
                (uint64_t)i->first.addr, (uint64_t)i->first.port, i->second);
    }
}
static void print_instance_map()
{
    if (!logging_debug(nnti_debug_level)) {
        return;
    }

    conn_by_inst_iter_t i;
    for (i=connections_by_instance.begin(); i != connections_by_instance.end(); i++) {
        log_debug(nnti_debug_level, "instance_map key=%llu conn=%p", i->first, i->second);
    }
}


static NNTI_result_t insert_buf_bufhash(NNTI_buffer_t *buf)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t h=hash6432shift((uint64_t)buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.buf);

    if (nthread_lock(&nnti_buf_bufhash_lock)) log_warn(nnti_debug_level, "failed to get thread lock");

    log_debug(nnti_debug_level, "adding buf=%p ; bufhash=%llu", buf, (uint64_t)h);

    assert(buffers_by_bufhash.find(h) == buffers_by_bufhash.end());
    buffers_by_bufhash[h] = buf;

    log_debug(nnti_debug_level, "bufhash buffer added (buf=%p ; buf.hash=%llu)", buf, (uint64_t)h);

    nthread_unlock(&nnti_buf_bufhash_lock);

    print_bufhash_map();

    return(rc);
}
static NNTI_buffer_t *get_buf_bufhash(const uint32_t bufhash)
{
    NNTI_buffer_t *buf=NULL;

    if (nthread_lock(&nnti_buf_bufhash_lock)) log_warn(nnti_debug_level, "failed to get thread lock");

    log_debug(nnti_debug_level, "looking for bufhash=%llu", (uint64_t)bufhash);

    if (buffers_by_bufhash.find(bufhash) != buffers_by_bufhash.end()) {
        buf = buffers_by_bufhash[bufhash];
    }

    log_debug(nnti_debug_level, "found buffer (buf=%p ; buf.hash=%llu)", buf, (uint64_t)bufhash);

    nthread_unlock(&nnti_buf_bufhash_lock);

    if (buf != NULL) {
        log_debug(nnti_debug_level, "buffer found (buf=%p ; buf.hash=%llu)", buf, (uint64_t)bufhash);
        return buf;
    }

    log_debug(nnti_debug_level, "buffer NOT found");

    print_bufhash_map();

    return(NULL);
}
static NNTI_buffer_t *del_buf_bufhash(NNTI_buffer_t *victim)
{
    NNTI_buffer_t *buf=NULL;
    uint32_t h=hash6432shift((uint64_t)victim->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.buf);
    log_level debug_level = nnti_debug_level;

    if (nthread_lock(&nnti_buf_bufhash_lock)) log_warn(nnti_debug_level, "failed to get thread lock");

    log_debug(debug_level, "deleting (victim=%p ; bufhash=%llu)", victim, (uint64_t)h);

    if (buffers_by_bufhash.find(h) != buffers_by_bufhash.end()) {
        buf = buffers_by_bufhash[h];
    }

    if (buf != NULL) {
        log_debug(debug_level, "buffer found and deleted (victim=%p ; buf=%p ; bufhash=%llu)", victim, buf, (uint64_t)h);
        buffers_by_bufhash.erase(h);
    } else {
        log_debug(debug_level, "buffer NOT found");
    }

    log_debug(debug_level, "deleted (victim=%p ; buf=%p ; bufhash=%llu)", victim, buf, (uint64_t)h);

    nthread_unlock(&nnti_buf_bufhash_lock);

    print_bufhash_map();

    return(buf);
}
static void print_bufhash_map()
{
    log_level debug_level=nnti_debug_level;

    if (!logging_debug(debug_level)) {
        return;
    }

    buf_by_bufhash_iter_t i;
    if (nthread_lock(&nnti_buf_bufhash_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    for (i=buffers_by_bufhash.begin(); i != buffers_by_bufhash.end(); i++) {
        log_debug(debug_level, "bufhash_map key=%llu buf=%p", i->first, i->second);
    }
    nthread_unlock(&nnti_buf_bufhash_lock);
}

static NNTI_result_t insert_wr_wrhash(nnti_gni_work_request_t *wr)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t h=hash6432shift((uint64_t)wr);

    if (nthread_lock(&nnti_wr_wrhash_lock)) log_warn(nnti_debug_level, "failed to get thread lock");

    log_debug(nnti_debug_level, "adding wrhash work request (wr=%p ; wr.hash=%llu)", wr, (uint64_t)h);

    assert(wr_by_wrhash.find(h) == wr_by_wrhash.end());
    wr_by_wrhash[h] = wr;

    log_debug(nnti_debug_level, "added wrhash work request (wr=%p ; wr.hash=%llu)", wr, (uint64_t)h);

    nthread_unlock(&nnti_wr_wrhash_lock);

    return(rc);

}
static nnti_gni_work_request_t *get_wr_wrhash(const uint32_t wrhash)
{
    nnti_gni_work_request_t *wr=NULL;

    if (nthread_lock(&nnti_wr_wrhash_lock)) log_warn(nnti_debug_level, "failed to get thread lock");

    log_debug(nnti_debug_level, "looking for wrhash=%llu", (uint64_t)wrhash);

    if (wr_by_wrhash.find(wrhash) != wr_by_wrhash.end()) {
        wr = wr_by_wrhash[wrhash];
    }

    log_debug(nnti_debug_level, "found work request (wr=%p ; wrhash=%llu)", wr, (uint64_t)wrhash);

    nthread_unlock(&nnti_wr_wrhash_lock);

    if (wr != NULL) {
        log_debug(nnti_debug_level, "work request found (wr=%p ; wrhash=%llu)", wr, (uint64_t)wrhash);
        return wr;
    }

    log_debug(nnti_debug_level, "work request NOT found");

    print_wrhash_map();

    return(NULL);
}
static nnti_gni_work_request_t *del_wr_wrhash(nnti_gni_work_request_t *victim)
{
    nnti_gni_work_request_t *wr=NULL;

    uint32_t h=hash6432shift((uint64_t)victim);
    log_level debug_level = nnti_debug_level;

    if (nthread_lock(&nnti_wr_wrhash_lock)) log_warn(nnti_debug_level, "failed to get thread lock");

    log_debug(debug_level, "deleting wrhash=%llu", (uint64_t)h);

    if (wr_by_wrhash.find(h) != wr_by_wrhash.end()) {
        wr = wr_by_wrhash[h];
    }

    if (wr != NULL) {
        log_debug(nnti_debug_level, "work request found and deleted (victim=%p ; wr=%p ; wrhash=%llu)", victim, wr, (uint64_t)h);
        wr_by_wrhash.erase(h);
    } else {
        log_debug(debug_level, "work request NOT found");
    }

    log_debug(debug_level, "deleted (wr=%p ; wrhash=%llu)", wr, (uint64_t)h);

    nthread_unlock(&nnti_wr_wrhash_lock);

    print_wrhash_map();

    return(wr);
}
static void print_wrhash_map()
{
    log_level debug_level=nnti_debug_level;

    if (!logging_debug(debug_level)) {
        return;
    }

    wr_by_wrhash_iter_t i;
    if (nthread_lock(&nnti_wr_wrhash_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    for (i=wr_by_wrhash.begin(); i != wr_by_wrhash.end(); i++) {
        log_debug(debug_level, "wrhash_map key=%llu wr=%p", i->first, i->second);
    }
    nthread_unlock(&nnti_wr_wrhash_lock);
}

static NNTI_result_t insert_sge_sgehash(nnti_gni_sge_t *sge)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t h=hash6432shift((uint64_t)sge);

    if (nthread_lock(&nnti_sge_sgehash_lock)) log_warn(nnti_debug_level, "failed to get thread lock");

    log_debug(nnti_debug_level, "adding sgehash (sge=%p ; sge.hash=%llu)", sge, (uint64_t)h);

    assert(sge_by_sgehash.find(h) == sge_by_sgehash.end());
    sge_by_sgehash[h] = sge;

    log_debug(nnti_debug_level, "added sgehash (sge=%p ; sge.hash=%llu)", sge, (uint64_t)h);

    nthread_unlock(&nnti_sge_sgehash_lock);

    return(rc);

}
static nnti_gni_sge_t *get_sge_sgehash(const uint32_t sgehash)
{
    nnti_gni_sge_t *sge=NULL;

    if (nthread_lock(&nnti_sge_sgehash_lock)) log_warn(nnti_debug_level, "failed to get thread lock");

    log_debug(nnti_debug_level, "looking for sgehash=%llu", (uint64_t)sgehash);

    if (sge_by_sgehash.find(sgehash) != sge_by_sgehash.end()) {
        sge = sge_by_sgehash[sgehash];
    }

    log_debug(nnti_debug_level, "found sge (sge=%p ; sgehash=%llu)", sge, (uint64_t)sgehash);

    nthread_unlock(&nnti_sge_sgehash_lock);

    if (sge != NULL) {
        log_debug(nnti_debug_level, "sge found (sge=%p ; sgehash=%llu)", sge, (uint64_t)sgehash);
        return sge;
    }

    log_debug(nnti_debug_level, "sge NOT found");

    print_sgehash_map();

    return(NULL);
}
static nnti_gni_sge_t *del_sge_sgehash(nnti_gni_sge_t *victim)
{
    nnti_gni_sge_t *sge=NULL;

    uint32_t h=hash6432shift((uint64_t)victim);
    log_level debug_level = nnti_debug_level;

    if (nthread_lock(&nnti_sge_sgehash_lock)) log_warn(nnti_debug_level, "failed to get thread lock");

    log_debug(debug_level, "deleting sgehash=%llu", (uint64_t)h);

    if (sge_by_sgehash.find(h) != sge_by_sgehash.end()) {
        sge = sge_by_sgehash[h];
    }

    if (sge != NULL) {
        log_debug(nnti_debug_level, "sge found and deleted (victim=%p ; sge=%p ; sgehash=%llu)", victim, sge, (uint64_t)h);
        sge_by_sgehash.erase(h);
    } else {
        log_debug(debug_level, "sge NOT found");
    }

    log_debug(debug_level, "deleted (sge=%p ; sgehash=%llu)", sge, (uint64_t)h);

    nthread_unlock(&nnti_sge_sgehash_lock);

    print_sgehash_map();

    return(sge);
}
static void print_sgehash_map()
{
    log_level debug_level=nnti_debug_level;

    if (!logging_debug(debug_level)) {
        return;
    }

    sge_by_sgehash_iter_t i;
    if (nthread_lock(&nnti_sge_sgehash_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    for (i=sge_by_sgehash.begin(); i != sge_by_sgehash.end(); i++) {
        log_debug(debug_level, "sgehash_map key=%llu sge=%p", i->first, i->second);
    }
    nthread_unlock(&nnti_sge_sgehash_lock);
}

static NNTI_result_t wr_pool_register(
        nnti_gni_work_request_t *wr,
        gni_cq_handle_t   cq_hdl)
{
    int rc=GNI_RC_SUCCESS;

    log_debug(nnti_debug_level, "enter");

//    rc=GNI_MemRegister (transport_global_data.nic_hdl,
//            (uint64_t)&wr->wc,
//            sizeof(nnti_gni_work_completion_t),
//            cq_hdl,
//            GNI_MEM_READWRITE,
//            (uint32_t)-1,
//            &wr->wc_mem_hdl);
//    if (rc!=GNI_RC_SUCCESS) {
//        log_error(nnti_debug_level, "MemRegister(wc_mem_hdl) failed: rc=%d, %s", rc, strerror(errno));
//        return(NNTI_EIO);
//    }
//
//    wr->wc_registered=FALSE;

    log_debug(nnti_debug_level, "exit");

    return(NNTI_OK);
}
static NNTI_result_t wr_pool_deregister(
        nnti_gni_work_request_t *wr)
{
    int rc=GNI_RC_SUCCESS;

    log_debug(nnti_debug_level, "enter");

//    if (wr->wc_registered==FALSE) {
//        log_debug(nnti_debug_level, "exit wr(%p) - not registered", wr);
//        return(NNTI_OK);
//    }
//
//    rc=GNI_MemDeregister(transport_global_data.nic_hdl, &wr->wc_mem_hdl);
//    if (rc!=GNI_RC_SUCCESS) {
//        log_error(nnti_debug_level, "MemDeregister(wc_mem_hdl) failed: rc=%d", rc);
//        return(NNTI_EIO);
//    }
//
//    wr->wc_registered=FALSE;

    log_debug(nnti_debug_level, "exit");

    return(NNTI_OK);
}
static NNTI_result_t wr_pool_init(void)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t i;
    nnti_gni_work_request_t *wr=NULL;

    log_debug(nnti_debug_level, "enter");

    for (i=0;i<config.wr_pool_initial_size;i++) {
        wr=(nnti_gni_work_request_t *)malloc(sizeof(nnti_gni_work_request_t));
        memset(wr, 0, sizeof(nnti_gni_work_request_t));
        assert(wr);
        nthread_lock_init(&wr->lock);
        rc=wr_pool_register(wr, transport_global_data.mem_cq_hdl);
        if (rc!=NNTI_OK) {
            log_error(nnti_debug_level, "failed to register target work request: rc=%d", rc);
            goto cleanup;
        }
        wr->is_initiator=FALSE;
        wr_pool_target_push(wr);

        wr=(nnti_gni_work_request_t *)malloc(sizeof(nnti_gni_work_request_t));
        memset(wr, 0, sizeof(nnti_gni_work_request_t));
        assert(wr);
        nthread_lock_init(&wr->lock);
        rc=wr_pool_register(wr, NULL);
        if (rc!=NNTI_OK) {
            log_error(nnti_debug_level, "failed to register initiator work request: rc=%d", rc);
            goto cleanup;
        }
        wr->is_initiator=TRUE;
        wr_pool_initiator_push(wr);
    }

cleanup:
    log_debug(nnti_debug_level, "exit");

    return(rc);
}
static nnti_gni_work_request_t *wr_pool_target_pop(void)
{
    nnti_gni_work_request_t *wr=NULL;

    log_debug(nnti_debug_level, "enter");

    if (nthread_lock(&nnti_wr_pool_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    if (!target_wr_pool.empty()) {
        wr=target_wr_pool.front();
        target_wr_pool.pop_front();
    }
    nthread_unlock(&nnti_wr_pool_lock);

    log_debug(nnti_debug_level, "exit");

    return(wr);
}
static nnti_gni_work_request_t *wr_pool_initiator_pop(void)
{
    nnti_gni_work_request_t *wr=NULL;

    log_debug(nnti_debug_level, "enter");

    if (nthread_lock(&nnti_wr_pool_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    if (!initiator_wr_pool.empty()) {
        wr=initiator_wr_pool.front();
        initiator_wr_pool.pop_front();
    }
    nthread_unlock(&nnti_wr_pool_lock);

    log_debug(nnti_debug_level, "exit");

    return(wr);
}
static void wr_pool_target_push(nnti_gni_work_request_t *wr)
{
    log_debug(nnti_debug_level, "enter");

    memset(&wr->wc, 0, sizeof(nnti_gni_work_completion_t));
    wr->state  =NNTI_GNI_WR_STATE_RESET;
    wr->last_op       =0;

    if (nthread_lock(&nnti_wr_pool_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    target_wr_pool.push_front(wr);
    nthread_unlock(&nnti_wr_pool_lock);

    log_debug(nnti_debug_level, "exit");

    return;
}
static void wr_pool_initiator_push(nnti_gni_work_request_t *wr)
{
    log_debug(nnti_debug_level, "enter");

    memset(&wr->wc, 0, sizeof(nnti_gni_work_completion_t));
    wr->state  =NNTI_GNI_WR_STATE_RESET;
    wr->last_op       =0;

    if (nthread_lock(&nnti_wr_pool_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    initiator_wr_pool.push_front(wr);
    nthread_unlock(&nnti_wr_pool_lock);

    log_debug(nnti_debug_level, "exit");

    return;
}
static NNTI_result_t wr_pool_fini(void)
{
    NNTI_result_t  rc=NNTI_OK;
    nnti_gni_work_request_t *wr=NULL;

    log_debug(nnti_debug_level, "enter");

    if (nthread_lock(&nnti_wr_pool_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    while (!target_wr_pool.empty()) {
        wr=target_wr_pool.front();
        target_wr_pool.pop_front();
        assert(wr);
        rc=wr_pool_deregister(wr);
        if (rc!=NNTI_OK) {
            log_error(nnti_debug_level, "failed to deregister target work request: rc=%d", rc);
            goto cleanup;
        }
        free(wr);
    }
    while (!initiator_wr_pool.empty()) {
        wr=initiator_wr_pool.front();
        initiator_wr_pool.pop_front();
        assert(wr);
        rc=wr_pool_deregister(wr);
        if (rc!=NNTI_OK) {
            log_error(nnti_debug_level, "failed to deregister initiator work request: rc=%d", rc);
            goto cleanup;
        }
        free(wr);
    }

cleanup:
    nthread_unlock(&nnti_wr_pool_lock);

    log_debug(nnti_debug_level, "exit");

    return(rc);
}

/*
 * Try hard to read the whole buffer.  Abort on read error.
 */
static int tcp_read(int sock, void *incoming, size_t len)
{
    int bytes_this_read=0;
    int bytes_left=len;
    int bytes_read=0;

    while (bytes_left > 0) {
        bytes_this_read = read(sock, (char *)incoming + bytes_read, bytes_left);
        if (bytes_this_read < 0) {
            return bytes_this_read;
        }
        if (bytes_this_read == 0) {
            break;
        }
        bytes_left -= bytes_this_read;
        bytes_read += bytes_this_read;
    }
    return bytes_read;
}

/*
 * Try hard to write the whole buffer.  Abort on write error.
 */
static int tcp_write(int sock, const void *outgoing, size_t len)
{
    int bytes_this_write=0;
    int bytes_left=len;
    int bytes_written=0;

    while (bytes_left > 0) {
        bytes_this_write = write(sock, (const char *)outgoing + bytes_written, bytes_left);
        if (bytes_this_write < 0) {
            return bytes_this_write;
        }
        bytes_left    -= bytes_this_write;
        bytes_written += bytes_this_write;
    }
    return bytes_written;
}

/*
 * Two processes exchange data over a TCP socket.  Both sides send and receive the
 * same amount of data.  Only one process can declare itself the server (is_server!=0),
 * otherwise this will hang because both will wait for the read to complete.
 *
 * Server receives, then sends.
 * Client sends, then receives.
 */
static int tcp_exchange(int sock, int is_server, void *incoming, void *outgoing, size_t len)
{
    int rc=0;

    if (is_server) {
        trios_declare_timer(callTime);
        trios_start_timer(callTime);
        rc = tcp_read(sock, incoming, len);
        trios_stop_timer("tcp_read", callTime);
        if (rc < 0) {
            log_warn(nnti_debug_level, "server failed to read GNI connection info: errno=%d", errno);
            goto out;
        }
        if (rc != (int) len) {
            log_error(nnti_debug_level, "partial read, %d/%d bytes", rc, (int) len);
            rc = 1;
            goto out;
        }
    } else {
        trios_declare_timer(callTime);
        trios_start_timer(callTime);
        rc = tcp_write(sock, outgoing, len);
        trios_stop_timer("tcp_write", callTime);
        if (rc < 0) {
            log_warn(nnti_debug_level, "client failed to write GNI connection info: errno=%d", errno);
            goto out;
        }
    }

    if (is_server) {
        trios_declare_timer(callTime);
        trios_start_timer(callTime);
        rc = tcp_write(sock, outgoing, len);
        trios_stop_timer("tcp_write", callTime);
        if (rc < 0) {
            log_warn(nnti_debug_level, "server failed to write GNI connection info: errno=%d", errno);
            goto out;
        }
    } else {
        trios_declare_timer(callTime);
        trios_start_timer(callTime);
        rc = tcp_read(sock, incoming, len);
        trios_stop_timer("tcp_read", callTime);
        if (rc < 0) {
            log_warn(nnti_debug_level, "client failed to read GNI connection info: errno=%d", errno);
            goto out;
        }
        if (rc != (int) len) {
            log_error(nnti_debug_level, "partial read, %d/%d bytes", rc, (int) len);
            rc = 1;
            goto out;
        }
    }

    rc = 0;

out:
    return rc;
}

static void transition_connection_to_ready(
        int sock,
        nnti_gni_connection_t *conn)
{
    int rc=NNTI_OK;
    trios_declare_timer(callTime);

    trios_start_timer(callTime);
    /* final sychronization to ensure both sides have posted RTRs */
    rc = tcp_exchange(sock, 0, &rc, &rc, sizeof(rc));
    trios_stop_timer("transition tcp_exchange", callTime);
}

static int new_client_connection(
        nnti_gni_connection_t *c,
        int sock)
{
    int rc;

    /*
     * Values passed through TCP to permit Gemini connection
     */
    struct {
        char             name[NNTI_HOSTNAME_LEN];
        uint32_t         addr;
        uint32_t         port;
        NNTI_instance_id instance;
        alpsAppGni_t     alps_info;

        uint64_t         atomics_addr;
        gni_mem_handle_t atomics_mem_hdl;

        uint32_t         client_has_recv_queue;
    } instance_in, instance_out;
    struct {
        nnti_gni_recv_queue_attrs_t server_attrs;
    } sa_in, sa_out;

    trios_declare_timer(call_time);

    c->connection_type=CLIENT_CONNECTION;

    /*
     * Prepare to exchange TCP and ALPS parameters with the server
     */
    memset(&instance_out, 0, sizeof(instance_out));
    strcpy(instance_out.name, transport_global_data.listen_name);
    instance_out.addr      = htonl((uint32_t)transport_global_data.listen_addr);
    instance_out.port      = htonl((uint32_t)transport_global_data.listen_port);
    instance_out.instance  = htonl(transport_global_data.instance);
    instance_out.alps_info = transport_global_data.alps_info;

    instance_out.atomics_addr   =(uint64_t)transport_global_data.atomics;
    instance_out.atomics_mem_hdl=transport_global_data.atomics_mem_hdl;

    if ((transport_global_data.req_queue.req_size  > 0) &&
        (transport_global_data.req_queue.req_count > 0)) {
        instance_out.client_has_recv_queue = 1;
    }

    /*
     * Exchange TCP and ALPS parameters with the server
     */
    trios_start_timer(call_time);
    rc = tcp_exchange(sock, 0, &instance_in, &instance_out, sizeof(instance_in));
    trios_stop_timer("tcp_exchange", call_time);
    if (rc)
        goto out;

    /*
     * Record the server's parameters
     */
    c->peer_name       = strdup(instance_in.name);
    c->peer_addr       = ntohl(instance_in.addr);
    c->peer_port       = ntohl(instance_in.port);
    c->peer_instance   = ntohl(instance_in.instance);
    c->peer_alps_info  = instance_in.alps_info;
    c->atomics_addr    = instance_in.atomics_addr;
    c->atomics_mem_hdl = instance_in.atomics_mem_hdl;

    /*
     * Read the receive queue attributes from the server
     */
    memset(&sa_in, 0, sizeof(sa_in));
    trios_start_timer(call_time);
    rc = tcp_read(sock, &sa_in, sizeof(sa_in));
    trios_stop_timer("read server queue attrs", call_time);
    if (rc == sizeof(sa_in)) {
        rc=0;
    }
    if (rc)
        goto out;

    /*
     * Record the server's attributes
     */
    c->queue_remote_attrs=sa_in.server_attrs;

    /*
     * Setup flow control attributes
     */
    client_req_queue_init(c);



    /*
     * If the client registered a receive queue (bidirectional connection)...
     */
    if (instance_out.client_has_recv_queue == 1) {
        memset(&sa_out, 0, sizeof(sa_out));

        nnti_gni_memory_handle_t *gni_mem_hdl=(nnti_gni_memory_handle_t *)transport_global_data.req_queue.reg_buf->transport_private;

        /*
         * Then write the receive queue attributes to the server
         */
        sa_out.server_attrs.req_index_addr   =transport_global_data.req_queue.req_index_addr;
        sa_out.server_attrs.req_index_mem_hdl=transport_global_data.req_queue.req_index_mem_hdl;
        sa_out.server_attrs.req_buffer_addr  =(uint64_t)transport_global_data.req_queue.req_buffer;
        sa_out.server_attrs.req_size         =transport_global_data.req_queue.req_size;
        sa_out.server_attrs.req_count        =transport_global_data.req_queue.req_count;
        sa_out.server_attrs.req_mem_hdl      =gni_mem_hdl->mem_hdl;

        trios_start_timer(call_time);
        rc = tcp_write(sock, &sa_out, sizeof(sa_out));
        trios_stop_timer("write server queue attrs", call_time);
        if (rc == sizeof(sa_out)) {
            rc=0;
        }
        if (rc)
            goto out;

    }


    rc=GNI_EpCreate (c->nic_hdl, transport_global_data.ep_cq_hdl, &c->ep_hdl);
    if (rc!=GNI_RC_SUCCESS) {
        log_error(nnti_debug_level, "EpCreate(c->ep_hdl) failed: %d", rc);
        goto out;
    }
    rc=GNI_EpBind (c->ep_hdl, c->peer_alps_info.local_addr, c->peer_instance);
    if (rc!=GNI_RC_SUCCESS) {
        log_error(nnti_debug_level, "EpBind(c->ep_hdl) failed: %d", rc);
        goto out;
    }

    log_debug(nnti_debug_level, "new connection ep_hdl(%llu)", c->ep_hdl);

out:
    return rc;
}

static int new_server_connection(
        nnti_gni_connection_t *c,
        int sock)
{
    int rc;

    /*
     * Values passed through TCP to permit Gemini connection
     */
    struct {
        char             name[NNTI_HOSTNAME_LEN];
        uint32_t         addr;
        uint32_t         port;
        NNTI_instance_id instance;
        alpsAppGni_t     alps_info;

        uint64_t         atomics_addr;
        gni_mem_handle_t atomics_mem_hdl;

        uint32_t         client_has_recv_queue;
    } instance_in, instance_out;
    struct {
        nnti_gni_recv_queue_attrs_t server_attrs;
    } sa_out, sa_in;

    trios_declare_timer(call_time);

    assert(transport_global_data.req_queue.reg_buf);

    nnti_gni_memory_handle_t *gni_mem_hdl=(nnti_gni_memory_handle_t *)transport_global_data.req_queue.reg_buf->transport_private;


    c->connection_type=SERVER_CONNECTION;

    /*
     * Prepare to exchange TCP and ALPS parameters with the client
     */
    memset(&instance_out, 0, sizeof(instance_out));
    strcpy(instance_out.name, transport_global_data.listen_name);
    instance_out.addr      = htonl((uint32_t)transport_global_data.listen_addr);
    instance_out.port      = htonl((uint32_t)transport_global_data.listen_port);
    instance_out.instance  = htonl(transport_global_data.instance);
    instance_out.alps_info = transport_global_data.alps_info;

    instance_out.atomics_addr   =(uint64_t)transport_global_data.atomics;
    instance_out.atomics_mem_hdl=transport_global_data.atomics_mem_hdl;

    /*
     * Exchange TCP and ALPS parameters with the client
     */
    trios_start_timer(call_time);
    rc = tcp_exchange(sock, 1, &instance_in, &instance_out, sizeof(instance_in));
    trios_stop_timer("tcp_exchange", call_time);
    if (rc)
        goto out;

    /*
     * Record the client's parameters
     */
    c->peer_name       = strdup(instance_in.name);
    c->peer_addr       = ntohl(instance_in.addr);
    c->peer_port       = ntohl(instance_in.port);
    c->peer_instance   = ntohl(instance_in.instance);
    c->peer_alps_info  = instance_in.alps_info;
    c->peer_ptag       = instance_in.alps_info.ptag;
    c->peer_cookie     = instance_in.alps_info.cookie;
    c->atomics_addr    = instance_in.atomics_addr;
    c->atomics_mem_hdl = instance_in.atomics_mem_hdl;

    c->cdm_hdl = transport_global_data.cdm_hdl;
    c->nic_hdl = transport_global_data.nic_hdl;

    /*
     * Write the receive queue attributes to the client
     */
    memset(&sa_out, 0, sizeof(sa_out));
    sa_out.server_attrs.req_index_addr   =transport_global_data.req_queue.req_index_addr;
    sa_out.server_attrs.req_index_mem_hdl=transport_global_data.req_queue.req_index_mem_hdl;
    sa_out.server_attrs.req_buffer_addr  =(uint64_t)transport_global_data.req_queue.req_buffer;
    sa_out.server_attrs.req_size         =transport_global_data.req_queue.req_size;
    sa_out.server_attrs.req_count        =transport_global_data.req_queue.req_count;
    sa_out.server_attrs.req_mem_hdl      =gni_mem_hdl->mem_hdl;

    trios_start_timer(call_time);
    rc = tcp_write(sock, &sa_out, sizeof(sa_out));
    trios_stop_timer("write server queue attrs", call_time);
    if (rc == sizeof(sa_out)) {
        rc=0;
    }
    if (rc)
        goto out;



    /*
     * If the client registered a receive queue (bidirectional connection)...
     */
    if (instance_in.client_has_recv_queue == 1) {
        /*
         * Then read the receive queue attributes from the client
         */
        memset(&sa_in, 0, sizeof(sa_in));
        trios_start_timer(call_time);
        rc = tcp_read(sock, &sa_in, sizeof(sa_in));
        trios_stop_timer("read server queue attrs", call_time);
        if (rc == sizeof(sa_in)) {
            rc=0;
        }
        if (rc)
            goto out;

        /*
         * Record the client's attributes
         */
        c->queue_remote_attrs=sa_in.server_attrs;

        /*
         * Setup flow control attributes
         */
        client_req_queue_init(c);

    }


    rc=GNI_EpCreate (c->nic_hdl, transport_global_data.ep_cq_hdl, &c->ep_hdl);
    if (rc!=GNI_RC_SUCCESS) {
        log_error(nnti_debug_level, "EpCreate(c->ep_hdl) failed: %d", rc);
        goto out;
    }
    rc=GNI_EpBind (c->ep_hdl, c->peer_alps_info.local_addr, c->peer_instance);
    if (rc!=GNI_RC_SUCCESS) {
        log_error(nnti_debug_level, "EpBind(c->ep_hdl) failed: %d", rc);
        goto out;
    }

    log_debug(nnti_debug_level, "new connection ep_hdl(%llu)", c->ep_hdl);

out:
    return rc;
}

/**
 * @brief initialize
 */
static NNTI_result_t init_connection(
        nnti_gni_connection_t **conn,
        const int sock,
        const int is_server)
{
    int rc=NNTI_OK; /* return code */

    trios_declare_timer(call_time);

    log_debug(nnti_debug_level, "initializing gni connection");

    trios_start_timer(call_time);
    if (is_server) {
        rc = new_server_connection(*conn, sock);
    } else {
        rc = new_client_connection(*conn, sock);
    }
    trios_stop_timer("new connection", call_time);
    if (rc) {
        close_connection(*conn);
        goto out;
    }

    print_gni_conn(*conn);

out:
    return((NNTI_result_t)rc);
}

static int check_for_waiting_connection()
{
    NNTI_result_t rc = NNTI_OK;

    struct sockaddr_in ssin;
    socklen_t len;
    int s;
    nnti_gni_connection_t *conn = NULL;
    NNTI_peer_t peer;

    len = sizeof(ssin);
    s = accept(transport_global_data.listen_sock, (struct sockaddr *) &ssin, &len);
    if (s < 0) {
        if ((errno == EAGAIN) || (errno == EWOULDBLOCK)) {
            log_debug(nnti_debug_level, "no connections waiting to be accepted: %s", strerror(errno));
            rc = NNTI_EWOULDBLOCK;
            goto cleanup;
        } else {
            log_error(nnti_debug_level, "failed to accept tcp socket connection: %s", strerror(errno));
            rc = NNTI_EIO;
            goto cleanup;
        }
    } else {
        conn = (nnti_gni_connection_t *)calloc(1, sizeof(nnti_gni_connection_t));
        log_debug(nnti_debug_level, "calloc returned conn=%p.", conn);
        if (conn == NULL) {
            log_error(nnti_debug_level, "calloc returned NULL.  out of memory?: %s", strerror(errno));
            rc=NNTI_ENOMEM;
            goto cleanup;
        }

//        nthread_lock(&nnti_gni_lock);
        rc=init_connection(&conn, s, 1);
        if (rc!=NNTI_OK) {
            goto cleanup;
        }
        create_peer(
                &peer,
                conn->peer_name,
                conn->peer_addr,
                conn->peer_port,
                conn->peer_ptag,
                conn->peer_cookie,
                conn->peer_instance);

        conn->peer=peer;

        del_conn_peer(&peer);
        del_conn_instance(conn->peer_instance);
        insert_conn_peer(&peer, conn);
        insert_conn_instance(conn->peer_instance, conn);

        transition_connection_to_ready(s, conn);
//        nthread_unlock(&nnti_gni_lock);

        log_debug(nnti_debug_level, "accepted new connection from %s:%u", conn->peer_name, conn->peer_port);

        if (close(s) < 0) {
            log_error(nnti_debug_level, "failed to close new tcp socket");
        }

        if (logging_debug(nnti_debug_level)) {
            fprint_NNTI_peer(logger_get_file(), "peer",
                    "end of check_listen_socket_for_new_connections", &peer);
        }
    }

cleanup:
    if (rc != NNTI_OK) {
        if (conn!=NULL) free(conn);
    }
    return rc;
}

/**
 * Check for new connections.  The listening socket is left nonblocking
 * so this test can be quick; but accept is not really that quick compared
 * to polling an Gemini interface, for instance.  Returns >0 if an accept worked.
 */
static int check_listen_socket_for_new_connections()
{
    bool done=false;
    while(!done) {
        if (check_for_waiting_connection() != NNTI_OK) {
            done=true;
        }
    }

    return(NNTI_OK);
}

static int init_server_listen_socket()
{
    NNTI_result_t rc=NNTI_OK;
    trios_declare_timer(call_time);
    int flags;
    struct hostent *host_entry;
    struct sockaddr_in skin;
    socklen_t skin_size=sizeof(struct sockaddr_in);

    trios_start_timer(call_time);
    transport_global_data.listen_sock = socket(AF_INET, SOCK_STREAM, 0);
    trios_stop_timer("socket", call_time);
    if (transport_global_data.listen_sock < 0)
        log_error(nnti_debug_level, "failed to create tcp socket: %s", strerror(errno));

    flags = 1;
    trios_start_timer(call_time);
    if (setsockopt(transport_global_data.listen_sock, SOL_SOCKET, SO_REUSEADDR, &flags, sizeof(flags)) < 0)
        log_error(nnti_debug_level, "failed to set tcp socket REUSEADDR flag: %s", strerror(errno));
    trios_stop_timer("setsockopt", call_time);

    flags=fcntl(transport_global_data.listen_sock, F_GETFL, 0);
    fcntl(transport_global_data.listen_sock, F_SETFL, flags | O_NONBLOCK);

    log_debug(nnti_debug_level, "listen_name (%s).", transport_global_data.listen_name);
    if (transport_global_data.listen_name[0]!='\0') {
        log_debug(nnti_debug_level, "using hostname from command-line (%s).", transport_global_data.listen_name);
    } else {
        trios_start_timer(call_time);
        gethostname(transport_global_data.listen_name, NNTI_HOSTNAME_LEN);
        trios_stop_timer("gethostname", call_time);
        log_debug(nnti_debug_level, "hostname not given on command-line.  using gethostname() result (%s).", transport_global_data.listen_name);
    }

    /* lookup the host provided on the command line */
    trios_start_timer(call_time);
    host_entry = gethostbyname(transport_global_data.listen_name);
    trios_stop_timer("gethostbyname", call_time);
    if (!host_entry) {
        log_warn(nnti_debug_level, "failed to resolve server name (%s): %s", transport_global_data.listen_name, strerror(errno));
        return NNTI_ENOENT;
    }

    memset(&skin, 0, sizeof(skin));
    skin.sin_family = AF_INET;
    memcpy(&skin.sin_addr, host_entry->h_addr_list[0], (size_t) host_entry->h_length);
    /* 0 here means to bind to a random port assigned by the kernel */
    skin.sin_port = 0;

retry:
    trios_start_timer(call_time);
    if (bind(transport_global_data.listen_sock, (struct sockaddr *)&skin, skin_size) < 0) {
        if (errno == EINTR) {
            goto retry;
        } else {
            log_error(nnti_debug_level, "failed to bind tcp socket: %s", strerror(errno));
        }
    }
    trios_stop_timer("bind", call_time);


    /* after the bind, get the "name" for the socket.  the "name" contains the port assigned by the kernel. */
    trios_start_timer(call_time);
    getsockname(transport_global_data.listen_sock, (struct sockaddr *)&skin, &skin_size);
    trios_stop_timer("getsockname", call_time);
    transport_global_data.listen_addr = (uint32_t)skin.sin_addr.s_addr;
    transport_global_data.listen_port = (uint16_t)skin.sin_port;
    log_debug(nnti_debug_level, "listening on ip(%s) addr(%u) port(%u)",
            transport_global_data.listen_name,
            (unsigned int)ntohl(skin.sin_addr.s_addr),
            (unsigned int)ntohs(skin.sin_port));
    trios_start_timer(call_time);
    if (listen(transport_global_data.listen_sock, 1024) < 0)
        log_error(nnti_debug_level, "failed to listen on tcp socket: %s", strerror(errno));
    trios_stop_timer("listen", call_time);

    return rc;
}

/*
 * At an explicit BYE message, or at finalize time, shut down a connection.
 * If descriptors are posted, defer and clean up the connection structures
 * later.
 */
static void close_connection(nnti_gni_connection_t *c)
{
    log_level debug_level = nnti_debug_level;  // nnti_ee_debug_level;

    if (c==NULL) return;

    log_debug(debug_level, "close_connection: start");

    print_gni_conn(c);

    if (c->peer_name) {
        free(c->peer_name);
        c->peer_name = NULL;
    }

    if (c->connection_type == CLIENT_CONNECTION) {
        client_req_queue_destroy(c);
    }

    c->state=DISCONNECTED;

    log_debug(debug_level, "close_connection: exit");
}

static void close_all_conn(void)
{
    log_level debug_level = nnti_debug_level;

    log_debug(debug_level, "enter (%d instance connections, %d peer connections)",
            connections_by_instance.size(), connections_by_peer.size());

    if (nthread_lock(&nnti_conn_instance_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    conn_by_inst_iter_t inst_iter = connections_by_instance.begin();
    while (inst_iter != connections_by_instance.end()) {
        log_debug(debug_level, "close connection (instance=%llu)", inst_iter->first);
        close_connection(inst_iter->second);

        connections_by_instance.erase(inst_iter++);
    }
    nthread_unlock(&nnti_conn_instance_lock);

    if (nthread_lock(&nnti_conn_peer_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    conn_by_peer_iter_t peer_iter = connections_by_peer.begin();
    while (peer_iter != connections_by_peer.end()) {
        log_debug(debug_level, "close connection (peer.addr=%llu)", peer_iter->first.addr);
//        close_connection(peer_iter->second);

        connections_by_peer.erase(peer_iter++);
    }
    nthread_unlock(&nnti_conn_peer_lock);

    log_debug(debug_level, "exit (%d instance connections, %d peer connections)",
            connections_by_instance.size(), connections_by_peer.size());

    return;
}

static uint32_t get_cpunum(void)
{
  int i, j;
  uint32_t cpu_num;

  cpu_set_t coremask;

  (void)sched_getaffinity(0, sizeof(coremask), &coremask);

  for (i = 0; i < CPU_SETSIZE; i++) {
    if (CPU_ISSET(i, &coremask)) {
      int run = 0;
      for (j = i + 1; j < CPU_SETSIZE; j++) {
        if (CPU_ISSET(j, &coremask)) run++;
        else break;
      }
      if (!run) {
        cpu_num=i;
      } else {
        fprintf(stdout, "This thread is bound to multiple CPUs(%d).  Using lowest numbered CPU(%d).", run+1, cpu_num);
        cpu_num=i;
      }
    }
  }
  return(cpu_num);
}

static void get_alps_info(
        alpsAppGni_t *alps_info)
{
    int alps_rc=0;
    int req_rc=0;
    size_t rep_size=0;

    uint64_t apid=0;
    alpsAppLLIGni_t *alps_info_list;
    char buf[1024];

    char *ptag_env_str=NULL;
    char *cookie_env_str=NULL;

    alps_info_list=(alpsAppLLIGni_t *)&buf[0];

    alps_app_lli_lock();

    log_debug(nnti_debug_level, "sending ALPS request");
    alps_rc = alps_app_lli_put_request(ALPS_APP_LLI_ALPS_REQ_GNI, NULL, 0);
    if (alps_rc != 0) log_debug(nnti_debug_level, "alps_app_lli_put_request failed: %d", alps_rc);
    log_debug(nnti_debug_level, "waiting for ALPS reply");
    alps_rc = alps_app_lli_get_response(&req_rc, &rep_size);
    if (alps_rc != 0) log_debug(nnti_debug_level, "alps_app_lli_get_response failed: alps_rc=%d", alps_rc);
    if (req_rc != 0) log_debug(nnti_debug_level, "alps_app_lli_get_response failed: req_rc=%d", req_rc);
    if (rep_size != 0) {
        log_debug(nnti_debug_level, "waiting for ALPS reply bytes (%d) ; sizeof(alps_info)==%d ; sizeof(alps_info_list)==%d", rep_size, sizeof(alps_info), sizeof(alps_info_list));
        alps_rc = alps_app_lli_get_response_bytes(alps_info_list, rep_size);
        if (alps_rc != 0) log_debug(nnti_debug_level, "alps_app_lli_get_response_bytes failed: %d", alps_rc);
    }

    log_debug(nnti_debug_level, "sending ALPS request");
    alps_rc = alps_app_lli_put_request(ALPS_APP_LLI_ALPS_REQ_APID, NULL, 0);
    if (alps_rc != 0) log_debug(nnti_debug_level, "alps_app_lli_put_request failed: %d", alps_rc);
    log_debug(nnti_debug_level, "waiting for ALPS reply");
    alps_rc = alps_app_lli_get_response(&req_rc, &rep_size);
    if (alps_rc != 0) log_debug(nnti_debug_level, "alps_app_lli_get_response failed: alps_rc=%d", alps_rc);
    if (req_rc != 0) log_debug(nnti_debug_level, "alps_app_lli_get_response failed: req_rc=%d", req_rc);
    if (rep_size != 0) {
        log_debug(nnti_debug_level, "waiting for ALPS reply bytes (%d) ; sizeof(apid)==%d", rep_size, sizeof(apid));
        alps_rc = alps_app_lli_get_response_bytes(&apid, rep_size);
        if (alps_rc != 0) log_debug(nnti_debug_level, "alps_app_lli_get_response_bytes failed: %d", alps_rc);
    }

    alps_app_lli_unlock();

    memcpy(alps_info, (alpsAppGni_t*)alps_info_list->u.buf, sizeof(alpsAppGni_t));
    transport_global_data.apid=apid;

    log_debug(nnti_debug_level, "apid                 =%llu", (unsigned long long)apid);
    log_debug(nnti_debug_level, "alps_info->device_id =%llu", (unsigned long long)alps_info->device_id);
    log_debug(nnti_debug_level, "alps_info->local_addr=%lld", (long long)alps_info->local_addr);
    log_debug(nnti_debug_level, "alps_info->cookie    =%llu", (unsigned long long)alps_info->cookie);
    log_debug(nnti_debug_level, "alps_info->ptag      =%llu", (unsigned long long)alps_info->ptag);

    log_debug(nnti_debug_level, "ALPS response - apid(%llu) alps_info->device_id(%llu) alps_info->local_addr(%llu) "
            "alps_info->cookie(%llu) alps_info->ptag(%llu)",
            (unsigned long long)apid,
            (unsigned long long)alps_info->device_id,
            (long long)alps_info->local_addr,
            (unsigned long long)alps_info->cookie,
            (unsigned long long)alps_info->ptag);

    if ((ptag_env_str=getenv("TRIOS_NNTI_GNI_PTAG")) != NULL) {
        uint32_t ptag_env=strtoul(ptag_env_str, NULL, 0);
        log_debug(nnti_debug_level, "replacing ALPS ptag (%lu) with user supplied ptag (%lu)", alps_info->ptag, ptag_env);
        alps_info->ptag=ptag_env;
    }
    if ((cookie_env_str=getenv("TRIOS_NNTI_GNI_COOKIE")) != NULL) {
        uint32_t cookie_env=strtoul(cookie_env_str, NULL, 0);
        log_debug(nnti_debug_level, "replacing ALPS cookie (%lu) with user supplied cookie (%lu)", alps_info->cookie, cookie_env);
        alps_info->cookie=cookie_env;
    }
    if ((ptag_env_str!=NULL) && (cookie_env_str==NULL)) {
        log_warn(nnti_debug_level, "user supplied ptag (%s) found without a user supplied cookie.  this probably won't work.", ptag_env_str);
    }
    if ((ptag_env_str==NULL) && (cookie_env_str!=NULL)) {
        log_warn(nnti_debug_level, "user supplied cookie (%s) found without a user supplied ptag.  this probably won't work.", cookie_env_str);
    }

    return;
}

static void print_wc(const nnti_gni_work_completion_t *wc)
{
    if (!logging_debug(nnti_debug_level)) {
        return;
    }

    log_debug(nnti_debug_level, "wc=%p, wc.op=%d, wc.inst_id=%llu, wc.byte_len=%llu, wc.byte_offset=%llu, wc.src_index=%llu, wc.src_offset=%llu, wc.dest_index=%llu, wc.dest_offset=%llu",
            wc,
            wc->op,
            wc->inst_id,
            wc->byte_len,
            wc->byte_offset,
            wc->src_index,
            wc->src_offset,
            wc->dest_index,
            wc->dest_offset);
}

static void print_cq_event(
        const gni_cq_entry_t *event,
        const bool            force)
{
    if (force) {
        log_debug(LOG_ALL, "event=%p, event.data=%llu, event.source=%llu, event.status=%llu, "
                "event.info=%llu, event.overrun=%llu, event.inst_id=%llu, event.tid=%llu, event.msg_id=%llu, event.type=%llu",
                event,
                (uint64_t)gni_cq_get_data(*event),
                (uint64_t)gni_cq_get_source(*event),
                (uint64_t)gni_cq_get_status(*event),
                (uint64_t)gni_cq_get_info(*event),
                (uint64_t)gni_cq_overrun(*event),
                (uint64_t)gni_cq_get_inst_id(*event),
                (uint64_t)gni_cq_get_tid(*event),
                (uint64_t)gni_cq_get_msg_id(*event),
                (uint64_t)gni_cq_get_type(*event));
    } else if (gni_cq_get_status(*event) != 0) {
        log_error(nnti_debug_level, "event=%p, event.data=%llu, event.source=%llu, event.status=%llu, "
                "event.info=%llu, event.overrun=%llu, event.inst_id=%llu, event.tid=%llu, event.msg_id=%llu, event.type=%llu",
                event,
                (uint64_t)gni_cq_get_data(*event),
                (uint64_t)gni_cq_get_source(*event),
                (uint64_t)gni_cq_get_status(*event),
                (uint64_t)gni_cq_get_info(*event),
                (uint64_t)gni_cq_overrun(*event),
                (uint64_t)gni_cq_get_inst_id(*event),
                (uint64_t)gni_cq_get_tid(*event),
                (uint64_t)gni_cq_get_msg_id(*event),
                (uint64_t)gni_cq_get_type(*event));
    } else {
        log_debug(nnti_debug_level, "event=%p, event.data=%llu, event.source=%llu, event.status=%llu, "
                "event.info=%llu, event.overrun=%llu, event.inst_id=%llu, event.tid=%llu, event.msg_id=%llu, event.type=%llu",
                event,
                (uint64_t)gni_cq_get_data(*event),
                (uint64_t)gni_cq_get_source(*event),
                (uint64_t)gni_cq_get_status(*event),
                (uint64_t)gni_cq_get_info(*event),
                (uint64_t)gni_cq_overrun(*event),
                (uint64_t)gni_cq_get_inst_id(*event),
                (uint64_t)gni_cq_get_tid(*event),
                (uint64_t)gni_cq_get_msg_id(*event),
                (uint64_t)gni_cq_get_type(*event));
    }
}

static void print_post_desc(
        const gni_post_descriptor_t *post_desc_ptr)
{
    if (!logging_debug(nnti_debug_level)) {
        return;
    }

    if (post_desc_ptr != NULL) {
        log_debug(nnti_debug_level, "post_desc_ptr                  ==%p", (uint64_t)post_desc_ptr);
        log_debug(nnti_debug_level, "post_desc_ptr->next_descr      ==%p", (uint64_t)post_desc_ptr->next_descr);
        log_debug(nnti_debug_level, "post_desc_ptr->prev_descr      ==%p", (uint64_t)post_desc_ptr->prev_descr);

        log_debug(nnti_debug_level, "post_desc_ptr->post_id         ==%llu", (uint64_t)post_desc_ptr->post_id);
        log_debug(nnti_debug_level, "post_desc_ptr->status          ==%llu", (uint64_t)post_desc_ptr->status);
        log_debug(nnti_debug_level, "post_desc_ptr->cq_mode_complete==%llu", (uint64_t)post_desc_ptr->cq_mode_complete);

        log_debug(nnti_debug_level, "post_desc_ptr->type            ==%llu", (uint64_t)post_desc_ptr->type);
        log_debug(nnti_debug_level, "post_desc_ptr->cq_mode         ==%llu", (uint64_t)post_desc_ptr->cq_mode);
        log_debug(nnti_debug_level, "post_desc_ptr->dlvr_mode       ==%llu", (uint64_t)post_desc_ptr->dlvr_mode);
        log_debug(nnti_debug_level, "post_desc_ptr->local_addr      ==%llu", (uint64_t)post_desc_ptr->local_addr);
        log_debug(nnti_debug_level, "post_desc_ptr->remote_addr     ==%llu", (uint64_t)post_desc_ptr->remote_addr);
        log_debug(nnti_debug_level, "post_desc_ptr->length          ==%llu", (uint64_t)post_desc_ptr->length);
        log_debug(nnti_debug_level, "post_desc_ptr->rdma_mode       ==%llu", (uint64_t)post_desc_ptr->rdma_mode);
        log_debug(nnti_debug_level, "post_desc_ptr->src_cq_hndl     ==%llu", (uint64_t)post_desc_ptr->src_cq_hndl);
        log_debug(nnti_debug_level, "post_desc_ptr->sync_flag_value ==%llu", (uint64_t)post_desc_ptr->sync_flag_value);
        log_debug(nnti_debug_level, "post_desc_ptr->sync_flag_addr  ==%llu", (uint64_t)post_desc_ptr->sync_flag_addr);
        log_debug(nnti_debug_level, "post_desc_ptr->amo_cmd         ==%llu", (uint64_t)post_desc_ptr->amo_cmd);
        log_debug(nnti_debug_level, "post_desc_ptr->first_operand   ==%llu", (uint64_t)post_desc_ptr->first_operand);
        log_debug(nnti_debug_level, "post_desc_ptr->second_operand  ==%llu", (uint64_t)post_desc_ptr->second_operand);
        log_debug(nnti_debug_level, "post_desc_ptr->cqwrite_value   ==%llu", (uint64_t)post_desc_ptr->cqwrite_value);
    } else {
        log_debug(nnti_debug_level, "post_desc_ptr == NULL");
    }
}

static void print_gni_conn(nnti_gni_connection_t *c)
{
    log_level debug_level=nnti_debug_level;

    if (!logging_debug(debug_level)) {
        return;
    }

    log_debug(debug_level, "c->peer_name       =%s", c->peer_name);
    log_debug(debug_level, "c->peer_addr       =%u", c->peer_addr);
    log_debug(debug_level, "c->peer_port       =%u", (uint32_t)c->peer_port);
    log_debug(debug_level, "c->peer_instance   =%llu", (uint64_t)c->peer_instance);

    log_debug(debug_level, "c->state           =%d", c->state);
}

//static void print_put_buf(void *buf, uint32_t size)
//{
//    struct data_t {
//            uint32_t int_val;
//            float float_val;
//            double double_val;
//    };
//
//    struct  data_array_t {
//            u_int data_array_t_len;
//            struct data_t *data_array_t_val;
//    };
//
////    struct data_array_t *da=(struct data_array_t *)buf;
//    const struct data_t *array = (struct data_t *)buf;
//    const int len = size/sizeof(struct data_t);
//    int idx=0;
//
//    for (idx=0;idx<len;idx++) {
//        log_debug(nnti_select_debug_level, "array[%d].int_val=%u ; array[%d].float_val=%f ; array[%d].double_val=%f",
//                idx, array[idx].int_val,
//                idx, array[idx].float_val,
//                idx, array[idx].double_val);
//    }
//
//}

static void print_raw_buf(void *buf, uint32_t size)
{
    if (logging_debug(nnti_debug_level)) {
        FILE* f=logger_get_file();
        uint64_t print_limit=(size<90) ? size : 90;
        fprintf(f, "\nbuf (%p)\n", buf);
        fflush(f);
        if (buf != NULL) {
            uint32_t l=0;
            for (l=0;l<print_limit;l++) {
                if (l%30 == 0) fprintf(f, "\nbuf (%lu) (offset(%u)) => ", (uint64_t)buf, l);
                fprintf(f, "%02hhX", ((char *)buf)[l]);
            }
            fprintf(f, "\n");
        }
    }
}

static uint16_t get_dlvr_mode_from_env()
{
    char *mode=getenv("NSSI_GNI_DELIVERY_MODE");

    log_debug(nnti_debug_level, "NSSI_GNI_DELIVERY_MODE=%s", mode);
    if ((mode==NULL) || !strcmp(mode, "GNI_DLVMODE_PERFORMANCE")) {
        log_debug(nnti_debug_level, "setting delivery mode to GNI_DLVMODE_PERFORMANCE");
        return GNI_DLVMODE_PERFORMANCE;
    } else if (!strcmp(mode, "GNI_DLVMODE_IN_ORDER")) {
        log_debug(nnti_debug_level, "setting delivery mode to GNI_DLVMODE_IN_ORDER");
        return GNI_DLVMODE_IN_ORDER;
    } else if (!strcmp(mode, "GNI_DLVMODE_NO_ADAPT")) {
        log_debug(nnti_debug_level, "setting delivery mode to GNI_DLVMODE_NO_ADAPT");
        return GNI_DLVMODE_NO_ADAPT;
    } else if (!strcmp(mode, "GNI_DLVMODE_NO_HASH")) {
        log_debug(nnti_debug_level, "setting delivery mode to GNI_DLVMODE_NO_HASH");
        return GNI_DLVMODE_NO_HASH;
    } else if (!strcmp(mode, "GNI_DLVMODE_NO_RADAPT")) {
        log_debug(nnti_debug_level, "setting delivery mode to GNI_DLVMODE_NO_RADAPT");
        return GNI_DLVMODE_NO_RADAPT;
    } else {
        log_debug(nnti_debug_level, "defaulting delivery mode to GNI_DLVMODE_PERFORMANCE");
        return GNI_DLVMODE_PERFORMANCE;
    }
}
static void set_dlvr_mode(
        gni_post_descriptor_t *pd)
{
    pd->dlvr_mode=transport_global_data.delivery_mode;
}

static void set_rdma_mode(
        gni_post_descriptor_t *pd)
{
    if (config.use_rdma_fence==true) {
        pd->rdma_mode=GNI_RDMAMODE_FENCE;
    }
}

static void set_post_desc(gni_post_descriptor_t *pd, nnti_gni_buffer_type_t target_buf_type, uint32_t buf_length)
{
    if (target_buf_type == GET_SRC_BUFFER) {
        // set type
        switch (config.rdma_mode) {
            case RDMA_FMA:
                pd->type=GNI_POST_FMA_GET;
                break;
            case RDMA_CROSSOVER:
                if (buf_length > config.fma_bte_crossover_size) {
                    pd->type=GNI_POST_RDMA_GET;
                } else {
                    pd->type=GNI_POST_FMA_GET;
                }
                break;
            case RDMA_BTE:
            case RDMA_MIXED:
            default:
                pd->type=GNI_POST_RDMA_GET;
                break;
        }
        // set cq_mode
        if (config.use_rdma_events==true) {
            pd->cq_mode=GNI_CQMODE_GLOBAL_EVENT|GNI_CQMODE_REMOTE_EVENT;
        } else {
            pd->cq_mode=GNI_CQMODE_GLOBAL_EVENT;
        }

    } else if (target_buf_type == PUT_DST_BUFFER) {
        // set type
        switch (config.rdma_mode) {
            case RDMA_FMA:
                pd->type=GNI_POST_FMA_PUT;
                break;
            case RDMA_CROSSOVER:
                if (buf_length > config.fma_bte_crossover_size) {
                    pd->type=GNI_POST_RDMA_PUT;
                } else {
                    pd->type=GNI_POST_FMA_PUT;
                }
                break;
            case RDMA_BTE:
            case RDMA_MIXED:
            default:
                pd->type=GNI_POST_RDMA_PUT;
                break;
        }
        // set cq_mode
        if (config.use_rdma_events==true) {
            pd->cq_mode=GNI_CQMODE_GLOBAL_EVENT|GNI_CQMODE_REMOTE_EVENT;
        } else {
            pd->cq_mode=GNI_CQMODE_GLOBAL_EVENT;
        }

    } else if (target_buf_type == RECEIVE_BUFFER) {
        // set type
        switch (config.rdma_mode) {
            case RDMA_CROSSOVER:
                if (buf_length > config.fma_bte_crossover_size) {
                    pd->type=GNI_POST_RDMA_PUT;
                } else {
                    pd->type=GNI_POST_FMA_PUT;
                }
                break;
            case RDMA_BTE:
                pd->type=GNI_POST_RDMA_PUT;
                break;
            case RDMA_FMA:
            case RDMA_MIXED:
            default:
                pd->type=GNI_POST_FMA_PUT;
                break;
        }
        // set cq_mode
        pd->cq_mode=GNI_CQMODE_GLOBAL_EVENT|GNI_CQMODE_REMOTE_EVENT;

    } else if (target_buf_type == REQUEST_BUFFER) {
        // set type
        switch (config.rdma_mode) {
            case RDMA_CROSSOVER:
                if (buf_length > config.fma_bte_crossover_size) {
                    pd->type=GNI_POST_RDMA_PUT;
                } else {
                    pd->type=GNI_POST_FMA_PUT;
                }
                break;
            case RDMA_BTE:
                pd->type=GNI_POST_RDMA_PUT;
                break;
            case RDMA_FMA:
            case RDMA_MIXED:
            default:
                pd->type=GNI_POST_FMA_PUT;
                break;
        }
        // set cq_mode
        pd->cq_mode=GNI_CQMODE_GLOBAL_EVENT|GNI_CQMODE_REMOTE_EVENT;
    } else {
        log_error(nnti_debug_level, "unknown target buffer type(%llu).", (uint64_t)target_buf_type);
    }

    set_dlvr_mode(pd);
    set_rdma_mode(pd);
}

static int post_wait(
        gni_cq_handle_t cq_hdl,
        int             timeout)
{
    int rc=0;
    trios_declare_timer(call_time);
    gni_post_descriptor_t *post_desc_ptr=NULL;
    gni_cq_entry_t ev_data;

    log_debug(nnti_ee_debug_level, "enter");

retry:
    memset(&ev_data, 0, sizeof(ev_data));
    log_debug(nnti_ee_debug_level, "calling CqWaitEvent");
//nthread_lock(&nnti_gni_lock);
    rc=GNI_CqWaitEvent_wrapper (cq_hdl, timeout, &ev_data);
    if (rc==GNI_RC_SUCCESS) {
        print_cq_event(&ev_data, false);
        log_debug(nnti_debug_level, "calling GetComplete");
        rc=GNI_GetCompleted_wrapper (cq_hdl, ev_data, &post_desc_ptr);
        if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "GetCompleted failed: %d", rc);
        print_post_desc(post_desc_ptr);
    } else if (rc==GNI_RC_TIMEOUT) {
        log_debug(nnti_debug_level, "CqWaitEvent timed out: %d", rc);
        goto retry;
    } else {
        log_error(nnti_debug_level, "CqWaitEvent failed: %d", rc);
    }
//nthread_unlock(&nnti_gni_lock);

    log_debug(nnti_ee_debug_level, "exit");

    return(rc);
}

static int reset_req_index(
        nnti_gni_request_queue_handle_t  *req_queue_attrs)
{
    int rc=0;
    gni_post_descriptor_t  post_desc;

    uint64_t value_before_reset=0;
    uint64_t value_before_reset_addr=(uint64_t)&value_before_reset;
    gni_mem_handle_t value_before_reset_mem_hdl;
    gni_cq_handle_t reset_cq_hdl;
    gni_ep_handle_t reset_ep_hdl;

    log_debug(nnti_ee_debug_level, "enter");

nthread_lock(&nnti_gni_lock);

    rc=GNI_CqCreate (transport_global_data.nic_hdl, 2, 0, GNI_CQ_BLOCKING, NULL, NULL, &reset_cq_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "CqCreate(value_before_reset_cq_hdl) failed: %d", rc);

//    rc=GNI_MemRegister (transport_global_data.nic_hdl, value_before_reset_addr, sizeof(uint64_t), NULL, GNI_MEM_READWRITE, (uint32_t)-1, &value_before_reset_mem_hdl);
//    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "MemRegister(value_before_reset) failed: %d", rc);

    rc=GNI_EpCreate (transport_global_data.nic_hdl, reset_cq_hdl, &reset_ep_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "EpCreate(reset_ep_hdl) failed: %d", rc);
    rc=GNI_EpBind (reset_ep_hdl, transport_global_data.alps_info.local_addr, transport_global_data.instance);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "EpBind(reset_ep_hdl) failed: %d", rc);

    log_debug(nnti_debug_level, "index before reset(%llu).", req_queue_attrs->req_index);


    memset(&post_desc, 0, sizeof(gni_post_descriptor_t));
    post_desc.type           =GNI_POST_AMO;
    post_desc.cq_mode        =GNI_CQMODE_SILENT;

    set_dlvr_mode(&post_desc);
    set_rdma_mode(&post_desc);

    post_desc.local_addr     =value_before_reset_addr;
    post_desc.local_mem_hndl =value_before_reset_mem_hdl;
    post_desc.remote_addr    =req_queue_attrs->req_index_addr;
    post_desc.remote_mem_hndl=req_queue_attrs->req_index_mem_hdl;
    post_desc.length         =sizeof(uint64_t);
    post_desc.amo_cmd        =GNI_FMA_ATOMIC_FAND;
    post_desc.first_operand  =0;

    log_debug(nnti_debug_level, "calling PostFma(reset index ep_hdl(%llu), local_addr=%llu, remote_addr=%llu)",
            reset_ep_hdl, post_desc.local_addr, post_desc.remote_addr);
    rc=GNI_PostFma(reset_ep_hdl, &post_desc);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "PostFma(reset index) failed: %d", rc);

//    post_wait(reset_cq_hdl, 1000);

//    rc=GNI_MemDeregister (transport_global_data.nic_hdl, &value_before_reset_mem_hdl);
//    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "MemDeregister(1) failed: %d", rc);

    rc=GNI_EpUnbind (reset_ep_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "EpUnbind() failed: %d", rc);
    rc=GNI_EpDestroy (reset_ep_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "EpDestroy() failed: %d", rc);

    rc=GNI_CqDestroy (reset_cq_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "CqDestroy() failed: %d", rc);

    log_debug(nnti_debug_level, "index after reset(%llu).", req_queue_attrs->req_index);

nthread_unlock(&nnti_gni_lock);

    log_debug(nnti_ee_debug_level, "exit");

    return(0);
}

static int fetch_add_buffer_offset(
        nnti_gni_client_queue_t       *local_req_queue_attrs,
        nnti_gni_recv_queue_attrs_t *remote_req_queue_attrs,
        uint64_t                     addend,
        uint64_t                    *prev_offset)
{
    int rc=0;
    trios_declare_timer(call_time);
    gni_post_descriptor_t  post_desc;
    gni_cq_entry_t ev_data;
    int extras_absorbed=0;

    log_level debug_level=nnti_debug_level;

//    static uint64_t last_offset=0;

    log_debug(nnti_ee_debug_level, "enter");

    memset(&post_desc, 0, sizeof(gni_post_descriptor_t));
    post_desc.type           =GNI_POST_AMO;
    post_desc.cq_mode        =GNI_CQMODE_GLOBAL_EVENT;

    set_dlvr_mode(&post_desc);
    set_rdma_mode(&post_desc);

    post_desc.local_addr     =local_req_queue_attrs->req_index_addr;
    post_desc.local_mem_hndl =local_req_queue_attrs->req_index_mem_hdl;
    post_desc.remote_addr    =remote_req_queue_attrs->req_index_addr;
    post_desc.remote_mem_hndl=remote_req_queue_attrs->req_index_mem_hdl;
    post_desc.length         =sizeof(uint64_t);
    post_desc.amo_cmd        =GNI_FMA_ATOMIC_FADD;
    post_desc.first_operand  =addend;
    post_desc.second_operand =0;

    log_debug(debug_level, "calling PostFma(fetch add - req_index_ep_hdl(%llu), local_addr=%llu, remote_addr=%llu)",
    		local_req_queue_attrs->req_index_ep_hdl, post_desc.local_addr, (uint64_t)post_desc.remote_addr);
    nthread_lock(&nnti_gni_lock);
    rc=GNI_PostFma(local_req_queue_attrs->req_index_ep_hdl, &post_desc);
    if (rc!=GNI_RC_SUCCESS) {
    	log_error(debug_level, "PostFma(fetch add) failed: %d", rc);
    	rc=-1;
    }

    rc=post_wait(local_req_queue_attrs->req_index_cq_hdl, 1000);
    nthread_unlock(&nnti_gni_lock);
    if (rc!=GNI_RC_SUCCESS) {
    	log_error(debug_level, "post_wait(fetch add) failed: %d", rc);
    	if (rc==GNI_RC_TRANSACTION_ERROR) {
    		log_error(debug_level, "post_wait(fetch add) return GNI_RC_TRANSACTION_ERROR.  This is unrecoverable.  Aborting.");
    		abort();
    	}
    	rc=-1;
    }

    log_debug(debug_level, "fetched queue_index(%llu)", local_req_queue_attrs->req_index);
    *prev_offset=local_req_queue_attrs->req_index;
    if (*prev_offset >= remote_req_queue_attrs->req_count) {
    	/*
    	 * The queue offset I just fetched is beyond the end of the queue.
    	 * The server needs to reset the offset counter.  I will return -2
    	 * to indicate the fetch returned an out-of-bounds value.
    	 */
    	rc=-1;
    } else {
    	local_req_queue_attrs->last_offset=*prev_offset;
    }

    log_debug(nnti_ee_debug_level, "exit");

    return(rc);
}

static int send_req(
        const NNTI_peer_t           *peer_hdl,
        nnti_gni_client_queue_t     *local_req_queue_attrs,
        nnti_gni_recv_queue_attrs_t *remote_req_queue_attrs,
        uint64_t                     offset,
        const NNTI_buffer_t         *reg_buf,
        nnti_gni_work_request_t     *gni_wr)
{
    int rc=0;
    nnti_gni_memory_handle_t *gni_mem_hdl=NULL;

    trios_declare_timer(call_time);

    bool use_fma=true;

    gni_mem_hdl=(nnti_gni_memory_handle_t *)reg_buf->transport_private;
    assert(gni_mem_hdl);

    log_debug(nnti_ee_debug_level, "enter");

    gni_wr->last_op =GNI_OP_SEND_REQUEST;
    gni_wr->state=NNTI_GNI_WR_STATE_STARTED;

    gni_wr->sge_list =&gni_wr->sge;
    gni_wr->sge_count=1;

    gni_wr->sge_list[0].gni_wr=gni_wr;
    gni_wr->sge_list[0].state=NNTI_GNI_SGE_STATE_STARTED;

    memset(&gni_wr->sge_list[0].post_desc, 0, sizeof(gni_post_descriptor_t));
    set_post_desc(&gni_wr->sge_list[0].post_desc, REQUEST_BUFFER, reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.size);

    gni_wr->sge_list[0].post_desc.src_cq_hndl=transport_global_data.ep_cq_hdl;

    gni_wr->sge_list[0].post_desc.local_addr           =reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.buf;
    gni_wr->sge_list[0].post_desc.local_mem_hndl.qword1=reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.mem_hdl.qword1;
    gni_wr->sge_list[0].post_desc.local_mem_hndl.qword2=reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.mem_hdl.qword2;
    gni_wr->sge_list[0].post_desc.remote_addr          =remote_req_queue_attrs->req_buffer_addr+offset;
    gni_wr->sge_list[0].post_desc.remote_mem_hndl      =remote_req_queue_attrs->req_mem_hdl;
    gni_wr->sge_list[0].post_desc.length               =reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.size + gni_mem_hdl->extra;

//    gni_wr->sge_list[0].post_desc.rdma_mode=GNI_RDMAMODE_FENCE;
    gni_wr->sge_list[0].post_desc.dlvr_mode=GNI_DLVMODE_IN_ORDER;

    print_raw_buf((void *)gni_wr->sge_list[0].post_desc.local_addr, gni_wr->sge_list[0].post_desc.length);

    insert_sge_sgehash(&gni_wr->sge_list[0]);

    if (config.rdma_mode==RDMA_CROSSOVER) {
        if (gni_wr->sge_list[0].post_desc.length > config.fma_bte_crossover_size) {
            use_fma=false;
        } else {
            use_fma=true;
        }
    } else if (config.rdma_mode==RDMA_BTE) {
        use_fma=false;
    } else if ((config.rdma_mode==RDMA_FMA) || (config.rdma_mode==RDMA_MIXED)) {
        use_fma=true;
    }

nthread_lock(&nnti_gni_lock);
    GNI_EpSetEventData(
            local_req_queue_attrs->req_ep_hdl,
            hash6432shift((uint64_t)&gni_wr->sge_list[0]),
            offset/(remote_req_queue_attrs->req_size+ + gni_mem_hdl->extra));

    if (use_fma) {
        log_debug(nnti_debug_level, "calling PostFma(send req ep_hdl(%llu), cq_hdl(%llu), local_addr=%llu, remote_addr=%llu)",
                local_req_queue_attrs->req_ep_hdl, transport_global_data.ep_cq_hdl, gni_wr->sge_list[0].post_desc.local_addr, gni_wr->sge_list[0].post_desc.remote_addr);
        trios_start_timer(call_time);
        rc=GNI_PostFma(local_req_queue_attrs->req_ep_hdl, &gni_wr->sge_list[0].post_desc);
        trios_stop_timer("PostFma req", call_time);
        if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "PostFma(fma put) failed: %d", rc);
    } else {
        log_debug(nnti_debug_level, "calling PostRdma(send req ep_hdl(%llu), cq_hdl(%llu), local_addr=%llu, remote_addr=%llu)",
                local_req_queue_attrs->req_ep_hdl, transport_global_data.ep_cq_hdl, gni_wr->sge_list[0].post_desc.local_addr, gni_wr->sge_list[0].post_desc.remote_addr);
        trios_start_timer(call_time);
        rc=GNI_PostRdma(local_req_queue_attrs->req_ep_hdl, &gni_wr->sge_list[0].post_desc);
        trios_stop_timer("PostRdma req", call_time);
        if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "PostRdma(rdma put) failed: %d", rc);
    }
nthread_unlock(&nnti_gni_lock);

    log_debug(nnti_ee_debug_level, "exit");

    return(0);
}

static int request_send(
        const NNTI_peer_t           *peer_hdl,
        nnti_gni_client_queue_t     *client_q,
        nnti_gni_recv_queue_attrs_t *server_q,
        const NNTI_buffer_t         *reg_buf,
        int                          req_num,
        nnti_gni_work_request_t     *gni_wr)
{
    int rc=0;
    uint64_t offset=0;
    trios_declare_timer(call_time);

    nnti_gni_memory_handle_t *gni_mem_hdl=NULL;

    log_debug(nnti_ee_debug_level, "enter");

    gni_mem_hdl=(nnti_gni_memory_handle_t *)reg_buf->transport_private;
    assert(gni_mem_hdl);

    log_debug(nnti_debug_level, "calling fetch_add_buffer_offset()");
    trios_start_timer(call_time);
    rc=fetch_add_buffer_offset(client_q, server_q, 1, &offset);
    trios_stop_timer("fetch_add_buffer_offset", call_time);
    if ((rc != -1) && (rc!=GNI_RC_SUCCESS)) log_error(nnti_debug_level, "fetch_add_buffer_offset() failed: %d", rc);

    if (rc == -1) {
    	log_debug(nnti_debug_level, "fetch_add_buffer failed: rc=%d", rc);
    } else {
    	gni_wr->reg_buf         =reg_buf;
    	gni_wr->wc              =GNI_WC_ADDRESS(reg_buf, 0);
    	gni_wr->wc->ack_received=0;
    	gni_wr->wc->inst_id     =transport_global_data.instance;
    	gni_wr->wc->op          =GNI_OP_SEND_REQUEST;
    	gni_wr->wc->byte_len    =reg_buf->payload_size;
    	gni_wr->wc->byte_offset =offset*(server_q->req_size+gni_mem_hdl->extra);
    	gni_wr->wc->src_offset  =0;
    	gni_wr->wc->dest_offset =offset*(server_q->req_size+gni_mem_hdl->extra);

    	nthread_lock(&gni_mem_hdl->wr_queue_lock);
    	gni_mem_hdl->wr_queue->push_back(gni_wr);
    	nthread_unlock(&gni_mem_hdl->wr_queue_lock);

//    	insert_wr_wrhash(gni_wr);

    	log_debug(nnti_debug_level, "calling send_req()");
    	trios_start_timer(call_time);
    	rc=send_req(peer_hdl, client_q, server_q, offset*(server_q->req_size+gni_mem_hdl->extra), reg_buf, gni_wr);
    	trios_stop_timer("send_req", call_time);
    	if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "send_req() failed: %d", rc);

    	gni_wr->last_op      =GNI_OP_SEND_REQUEST;
    	gni_wr->peer_instance=peer_hdl->peer.NNTI_remote_process_t_u.gni.inst_id;
    	log_debug(nnti_event_debug_level, "gni_wr->last_peer==%lu", gni_wr->peer_instance);
    }

    log_debug(nnti_ee_debug_level, "exit");

    return(rc);
}

static int send_buffer(
        const NNTI_peer_t       *peer_hdl,
        const NNTI_buffer_t     *src_hdl,
        const NNTI_buffer_t     *dest_hdl,
        nnti_gni_work_request_t *gni_wr)
{
    int rc=0;

    gni_ep_handle_t        ep_hdl;

    nnti_gni_memory_handle_t *gni_mem_hdl=NULL;

    bool use_fma=true;

    trios_declare_timer(call_time);

    gni_mem_hdl=(nnti_gni_memory_handle_t *)src_hdl->transport_private;
    assert(gni_mem_hdl);

    log_debug(nnti_ee_debug_level, "enter");

    gni_wr->last_op=GNI_OP_SEND_BUFFER;
    gni_wr->state=NNTI_GNI_WR_STATE_STARTED;

    gni_wr->sge_list =&gni_wr->sge;
    gni_wr->sge_count=1;

    gni_wr->sge_list[0].gni_wr=gni_wr;
    gni_wr->sge_list[0].state=NNTI_GNI_SGE_STATE_STARTED;

    memset(&gni_wr->sge_list[0].post_desc, 0, sizeof(gni_post_descriptor_t));
    set_post_desc(&gni_wr->sge_list[0].post_desc, RECEIVE_BUFFER, src_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.size);

    gni_wr->sge_list[0].post_desc.src_cq_hndl=transport_global_data.ep_cq_hdl;

    gni_wr->sge_list[0].post_desc.local_addr            =src_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.buf;
    gni_wr->sge_list[0].post_desc.local_mem_hndl.qword1 =src_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.mem_hdl.qword1;
    gni_wr->sge_list[0].post_desc.local_mem_hndl.qword2 =src_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.mem_hdl.qword2;
    gni_wr->sge_list[0].post_desc.remote_addr           =dest_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.buf;
    gni_wr->sge_list[0].post_desc.remote_mem_hndl.qword1=dest_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.mem_hdl.qword1;
    gni_wr->sge_list[0].post_desc.remote_mem_hndl.qword2=dest_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.mem_hdl.qword2;
    gni_wr->sge_list[0].post_desc.length                =src_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.size + gni_mem_hdl->extra;

//    gni_wr->sge_list[0].post_desc.rdma_mode=GNI_RDMAMODE_FENCE;
    gni_wr->sge_list[0].post_desc.dlvr_mode=GNI_DLVMODE_IN_ORDER;

    print_raw_buf((void *)gni_wr->sge_list[0].post_desc.local_addr, gni_wr->sge_list[0].post_desc.length);

    insert_sge_sgehash(&gni_wr->sge_list[0]);

    nnti_gni_connection_t *conn=get_conn_peer(peer_hdl);
    assert(conn);
    ep_hdl=conn->ep_hdl;

    if (config.rdma_mode==RDMA_CROSSOVER) {
        if (gni_wr->sge_list[0].post_desc.length > config.fma_bte_crossover_size) {
            use_fma=false;
        } else {
            use_fma=true;
        }
    } else if (config.rdma_mode==RDMA_BTE) {
        use_fma=false;
    } else if ((config.rdma_mode==RDMA_FMA) || (config.rdma_mode==RDMA_MIXED)) {
        use_fma=true;
    }

nthread_lock(&nnti_gni_lock);
    GNI_EpSetEventData(
            ep_hdl,
            hash6432shift((uint64_t)&gni_wr->sge_list[0]),
            hash6432shift((uint64_t)dest_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.gni.buf));

    if (use_fma) {
        log_debug(nnti_debug_level, "calling PostFma(send buffer ep_hdl(%llu), ep_cq_hdl(%llu), local_addr=%llu, remote_addr=%llu)",
                ep_hdl, transport_global_data.ep_cq_hdl, gni_wr->sge_list[0].post_desc.local_addr, gni_wr->sge_list[0].post_desc.remote_addr);
        trios_start_timer(call_time);
        rc=GNI_PostFma(ep_hdl, &gni_wr->sge_list[0].post_desc);
        trios_stop_timer("PostFma req", call_time);
        if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "PostFma(fma put) failed: %d", rc);
    } else {
        log_debug(nnti_debug_level, "calling PostRdma(send buffer ep_hdl(%llu), ep_cq_hdl(%llu), local_addr=%llu, remote_addr=%llu)",
                ep_hdl, transport_global_data.ep_cq_hdl, gni_wr->sge_list[0].post_desc.local_addr, gni_wr->sge_list[0].post_desc.remote_addr);
        trios_start_timer(call_time);
        rc=GNI_PostRdma(ep_hdl, &gni_wr->sge_list[0].post_desc);
        trios_stop_timer("PostRdma req", call_time);
        if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "PostRdma(rdma put) failed: %d", rc);
    }
nthread_unlock(&nnti_gni_lock);

    log_debug(nnti_ee_debug_level, "exit");

    return(0);
}

static int buffer_send(
        const NNTI_peer_t       *peer_hdl,
        const NNTI_buffer_t     *src_hdl,
        const NNTI_buffer_t     *dest_hdl,
        nnti_gni_work_request_t *gni_wr)
{
    int rc=0;
    trios_declare_timer(call_time);

    nnti_gni_memory_handle_t *gni_mem_hdl=NULL;


    gni_mem_hdl=(nnti_gni_memory_handle_t *)src_hdl->transport_private;
    assert(gni_mem_hdl);

    log_debug(nnti_ee_debug_level, "enter");

    gni_wr->reg_buf         =src_hdl;
    gni_wr->wc              =GNI_WC_ADDRESS(src_hdl, 0);
    gni_wr->wc->ack_received=0;
    gni_wr->wc->inst_id     =transport_global_data.instance;
    gni_wr->wc->op          =GNI_OP_SEND_BUFFER;
    gni_wr->wc->byte_len    =src_hdl->payload_size;
    gni_wr->wc->byte_offset =0;
    gni_wr->wc->src_offset  =0;
    gni_wr->wc->dest_offset =0;


    nthread_lock(&gni_mem_hdl->wr_queue_lock);
    gni_mem_hdl->wr_queue->push_back(gni_wr);
    nthread_unlock(&gni_mem_hdl->wr_queue_lock);

//    insert_wr_wrhash(gni_wr);

    log_debug(nnti_debug_level, "calling send_buffer()");
    trios_start_timer(call_time);
    rc=send_buffer(peer_hdl, src_hdl, dest_hdl, gni_wr);
    trios_stop_timer("send_buffer", call_time);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "send_buffer() failed: %d", rc);

    gni_wr->last_op      =GNI_OP_SEND_BUFFER;
    gni_wr->peer_instance=dest_hdl->buffer_owner.peer.NNTI_remote_process_t_u.gni.inst_id;
    log_debug(nnti_event_debug_level, "gni_wr->last_peer==%lu", gni_wr->peer_instance);

    log_debug(nnti_ee_debug_level, "exit");

    return(0);
}

static void client_req_queue_init(
        nnti_gni_connection_t *c)
{
    int rc;

    nnti_gni_client_queue_t  *q              =&c->queue_local_attrs;
    alpsAppGni_t           *server_params  =&c->peer_alps_info;
    uint64_t                server_instance=c->peer_instance;

    log_debug(nnti_debug_level, "enter");

    q->last_offset=0;

    q->req_index=0;
    q->req_index_addr=(uint64_t)&q->req_index;

    log_debug(nnti_debug_level, "client_req_queue->req_index_addr=%llu", (uint64_t)q->req_index_addr);

    rc=GNI_CqCreate (
            c->nic_hdl,
            1,
            0,
            GNI_CQ_BLOCKING,
            NULL,
            NULL,
            &q->req_index_cq_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "CqCreate() failed: %d", rc);

    rc=GNI_MemRegister (
            c->nic_hdl,
            q->req_index_addr,
            sizeof(uint64_t),
            NULL,
            GNI_MEM_READWRITE,
            (uint32_t)-1,
            &q->req_index_mem_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "MemRegister(1) failed: %d", rc);

    rc=GNI_EpCreate (
            c->nic_hdl,
            transport_global_data.ep_cq_hdl,
            &q->req_ep_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "EpCreate(req_ep_hdl) failed: %d", rc);
    rc=GNI_EpBind (
            q->req_ep_hdl,
            server_params->local_addr,
            server_instance);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "EpBind(req_ep_hdl) failed: %d", rc);

    rc=GNI_EpCreate (
            c->nic_hdl,
            q->req_index_cq_hdl,
            &q->req_index_ep_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "EpCreate(req_index_ep_hdl) failed: %d", rc);
    rc=GNI_EpBind (
            q->req_index_ep_hdl,
            server_params->local_addr,
            server_instance);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "EpBind(req_index_ep_hdl) failed: %d", rc);

    log_debug(nnti_debug_level, "exit");
}

static void client_req_queue_destroy(
        nnti_gni_connection_t *c)
{
    int rc;
    log_level debug_level = nnti_debug_level;  // nnti_debug_level;

    log_debug(debug_level, "enter");

    nnti_gni_client_queue_t     *q=&c->queue_local_attrs;

    rc=GNI_MemDeregister (c->nic_hdl, &q->req_index_mem_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "MemDeregister(1) failed: %d", rc);

    rc=GNI_EpUnbind (q->req_ep_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "EpUnbind() failed: %d", rc);
    rc=GNI_EpDestroy (q->req_ep_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "EpDestroy() failed: %d", rc);

    rc=GNI_EpUnbind (q->req_index_ep_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "EpUnbind() failed: %d", rc);
    rc=GNI_EpDestroy (q->req_index_ep_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "EpDestroy() failed: %d", rc);

    rc=GNI_CqDestroy (q->req_index_cq_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "CqDestroy() failed: %d", rc);

    log_debug(debug_level, "exit");
}

static void server_req_queue_init(
        nnti_gni_request_queue_handle_t *q,
        char                     *buffer,
        uint64_t                  req_size,
        uint64_t                  extra,
        uint64_t                  req_count)
{
    int rc;
    nnti_gni_memory_handle_t *gni_mem_hdl=(nnti_gni_memory_handle_t *)q->reg_buf->transport_private;

    log_debug(nnti_debug_level, "enter");

    q->req_buffer     =buffer;
    q->req_size       =req_size;
    q->req_count      =req_count;
    q->req_buffer_size=req_count*(req_size+extra);

    q->req_index     =0;
    q->req_index_addr=(uint64_t)&q->req_index;

    log_debug(nnti_debug_level, "server_req_queue->req_index_addr=%llu", (uint64_t)q->req_index_addr);

    q->head=0;
    q->tail=0;

    rc=GNI_MemRegister (
            transport_global_data.nic_hdl,
            (uint64_t)q->req_buffer,
            q->req_buffer_size,
            transport_global_data.mem_cq_hdl,
            GNI_MEM_READWRITE,
            (uint32_t)-1,
            &gni_mem_hdl->mem_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "MemRegister(1) failed: %d", rc);
    rc=GNI_MemRegister (
            transport_global_data.nic_hdl,
            q->req_index_addr,
            sizeof(uint64_t),
            NULL,
            GNI_MEM_READWRITE,
            (uint32_t)-1,
            &q->req_index_mem_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "MemRegister(1) failed: %d", rc);

    log_debug(nnti_debug_level, "exit");
}

static void server_req_queue_destroy(
        nnti_gni_request_queue_handle_t *q)
{
    int rc;
    nnti_gni_memory_handle_t *gni_mem_hdl=(nnti_gni_memory_handle_t *)q->reg_buf->transport_private;

    log_debug(nnti_debug_level, "enter");

    rc=GNI_MemDeregister (
            transport_global_data.nic_hdl,
            &gni_mem_hdl->mem_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "MemDeregister(1) failed: %d", rc);
    rc=GNI_MemDeregister (
            transport_global_data.nic_hdl,
            &q->req_index_mem_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "MemDeregister(1) failed: %d", rc);

    log_debug(nnti_debug_level, "exit");
}


#define EP_CQ_INDEX        0
#define MEM_CQ_INDEX       1
#define INTERRUPT_CQ_INDEX 2
#define CQ_COUNT 3
static NNTI_result_t progress(int timeout)
{
    int           rc     =0;
    gni_return_t  gni_rc =GNI_RC_SUCCESS;
    NNTI_result_t nnti_rc=NNTI_OK;

    nnti_gni_memory_handle_t        *gni_mem_hdl=NULL;
    nnti_gni_sge_t                  *gni_sge=NULL;
    nnti_gni_request_queue_handle_t *q_hdl=NULL;

    const NNTI_buffer_t      *wait_buf=NULL;

    gni_cq_handle_t  cq_list[CQ_COUNT];
    gni_cq_handle_t *cq;
    uint32_t         which_cq=0;
    gni_cq_entry_t   ev_data;


    gni_post_descriptor_t   *post_desc_ptr=NULL;

    long entry_time  =trios_get_time_ms();
    long elapsed_time=0;
    long gni_timeout=0;

    log_level debug_level  =nnti_debug_level;
    log_level old_log_level=logger_get_default_level();

    static bool in_progress  =false;   // if true, another thread is already making progress.
    bool        made_progress=false;
    bool        interrupted  =false;

    trios_declare_timer(call_time);
    trios_declare_timer(total_time);


    log_debug(debug_level, "enter (timeout=%d)", timeout);

    trios_start_timer(total_time);

    cq_list[EP_CQ_INDEX]       =transport_global_data.ep_cq_hdl;
    cq_list[MEM_CQ_INDEX]      =transport_global_data.mem_cq_hdl;
    cq_list[INTERRUPT_CQ_INDEX]=transport_global_data.interrupt_mem_cq_hdl;

    /*
     * Only one thread is allowed to make progress at a time.  All others
     * wait for the progress maker to finish, then everyone returns at once.
     */
    nthread_lock(&nnti_progress_lock);
    if (!in_progress) {
        log_debug(debug_level, "making progress");
        // no other thread is making progress.  we'll do it.
        in_progress=true;
        log_debug(debug_level, "set in_progress=true");
        nthread_unlock(&nnti_progress_lock);
    } else {
        // another thread is making progress.  we'll wait until they are done.
        rc=0;
        elapsed_time=0;
        if (in_progress) {
            if (timeout > 0) {
                log_debug(debug_level, "waiting for progress with timeout=%d", timeout-elapsed_time);
                // wait for progress or until timeout
                rc=nthread_timedwait(&nnti_progress_cond, &nnti_progress_lock, timeout-elapsed_time);
            } else if (timeout < 0) {
                log_debug(debug_level, "waiting infinitely for progress");
                // infinite wait for progress
                rc=nthread_wait(&nnti_progress_cond, &nnti_progress_lock);
            } else {
                log_debug(debug_level, "waiting for progress with timeout=0.  immediate timeout.");
                // timeout == 0 and we are not the progress maker.  report a timeout.
                rc=ETIMEDOUT;
                nthread_lock(&nnti_progress_lock);
            }
            elapsed_time = (trios_get_time_ms() - entry_time);
            log_debug(debug_level, "rc=%d, elapsed_time=%d", rc, elapsed_time);
        }
        nthread_unlock(&nnti_progress_lock);
        if (rc == ETIMEDOUT) {
            log_debug(debug_level, "timed out waiting for progress");
            nnti_rc = NNTI_ETIMEDOUT;
        } else if (rc == 0) {
            log_debug(debug_level, "someone made progress");
            nnti_rc=NNTI_OK;
        }
        goto cleanup;
    }

    while ((!made_progress) && (!interrupted))   {

        if (trios_exit_now()) {
            log_debug(debug_level, "caught abort signal");
            nnti_rc=NNTI_ECANCELED;
            break;
        }

        check_listen_socket_for_new_connections();

        memset(&ev_data, 0, sizeof(gni_cq_entry_t));
        post_desc_ptr=NULL;


        if (timeout-elapsed_time < 0) {
        	gni_timeout=MIN_TIMEOUT;
        } else {
        	gni_timeout=timeout-elapsed_time;
        }

        log_debug(nnti_cq_debug_level, "checking for events on any CQ (timeout=%d)", gni_timeout);
        trios_start_timer(call_time);
//        nthread_lock(&nnti_gni_lock);
        gni_rc=GNI_CqVectorWaitEvent_wrapper(cq_list, CQ_COUNT, gni_timeout, &ev_data, &which_cq);
        if (gni_cq_get_source(ev_data) == 1) {
            DEQUEUE_POST_DESCRIPTOR(cq_list[which_cq], &ev_data, &post_desc_ptr);
        }
//        nthread_unlock(&nnti_gni_lock);
        trios_stop_timer("progress - CqVectorWaitEvent", call_time);

        print_cq_event(&ev_data, false);
        print_post_desc(post_desc_ptr);
        log_debug(nnti_event_debug_level, "CqWaitEvent(progress) complete");
        if (gni_rc!=GNI_RC_SUCCESS) log_debug(nnti_debug_level, "CqWaitEvent() on cq_list(%p) failed: %d", cq_list, gni_rc);

        /* case 1: success */
        if (gni_rc == GNI_RC_SUCCESS) {
        	if (which_cq == INTERRUPT_CQ_INDEX) {
                log_debug(nnti_debug_level, "GNI_CqVectorWaitEvent() interrupted by NNTI_gni_interrupt");
        		interrupted=true;
                nnti_rc = NNTI_EINTR;
        		continue;
        	}
            gni_sge=decode_sge(&ev_data);
            process_event(gni_sge->gni_wr->reg_buf, gni_sge, cq_list[which_cq], &ev_data, post_desc_ptr);

            nnti_rc = NNTI_OK;

            made_progress=true;
        }
        /* case 2: timed out */
        else if ((gni_rc==GNI_RC_TIMEOUT) || (gni_rc==GNI_RC_NOT_DONE)) {
            elapsed_time = (trios_get_time_ms() - entry_time);

            /* if the caller asked for a legitimate timeout, we need to exit */
            if (((timeout >= 0) && (elapsed_time >= timeout)) || trios_exit_now()) {
                log_debug(nnti_debug_level, "CqVectorWaitEvent timed out...timeout(%d) elapsed_time(%d) exit_now(%d)",
                        timeout, elapsed_time, trios_exit_now());
                nnti_rc = NNTI_ETIMEDOUT;
                break;
            }
            /* continue if the timeout has not expired */
            log_debug(nnti_event_debug_level, "CqWaitEvent timedout...retrying...timeout(%d) elapsed_time(%d) exit_now(%d)",
                    timeout, elapsed_time, trios_exit_now());

            continue;
        }
        /* case 3: failure */
        else {
            char errstr[1024];
            uint32_t recoverable=0;
            errstr[0]='\0';
            GNI_CqErrorStr(ev_data, errstr, 1023);
            GNI_CqErrorRecoverable(ev_data, &recoverable);

            log_error(nnti_debug_level, "CqWaitEvent failed (gni_rc=%d) (recoverable=%llu) : %s",
                    gni_rc, (uint64_t)recoverable, errstr);
            print_cq_event(&ev_data, false);

            gni_sge=decode_sge(&ev_data);
            process_event(gni_sge->gni_wr->reg_buf, gni_sge, cq_list[which_cq], &ev_data, post_desc_ptr);

            nnti_rc = NNTI_EIO;

            made_progress=true;
        }
    }

unlock:
    nthread_lock(&nnti_progress_lock);
    in_progress=false;
    log_debug(debug_level, "set in_progress=false");
    nthread_broadcast(&nnti_progress_cond);
    log_debug(debug_level, "broadcasted on nnti_progress_cond");
    nthread_unlock(&nnti_progress_lock);

cleanup:
    trios_stop_timer("progress", total_time);

    log_debug(debug_level, "exit");

    return(nnti_rc);
}



#define LISTEN_IFACE_BASENAME "ipog"
static NNTI_result_t get_ipaddr(
        char *ipaddr,
        int maxlen)
{
    struct ifaddrs * ifAddrStruct=NULL;
    struct ifaddrs * ifa=NULL;
    void * tmpAddrPtr=NULL;

    getifaddrs(&ifAddrStruct);
    for (ifa = ifAddrStruct; ifa != NULL; ifa = ifa->ifa_next) {
        ipaddr[0]='\0';
        if (ifa ->ifa_addr->sa_family==AF_INET) { // check it is IP4
            // is a valid IP4 Address
            tmpAddrPtr=&((struct sockaddr_in *)ifa->ifa_addr)->sin_addr;
            log_debug(nnti_debug_level, "checking iface (IPv4) name (%s)", ifa->ifa_name);
            inet_ntop(AF_INET, tmpAddrPtr, ipaddr, maxlen);
            log_debug(nnti_debug_level, "hostname(%s) has IP Address %s", ifa->ifa_name, ipaddr);
            if (0==strncmp(ifa->ifa_name, LISTEN_IFACE_BASENAME, strlen(LISTEN_IFACE_BASENAME))) {
                log_debug(nnti_debug_level, "hostname(%s) matches", ifa->ifa_name);
                break;
            }
        } else if (ifa->ifa_addr->sa_family==AF_INET6) { // check it is IP6
            // is a valid IP6 Address
            tmpAddrPtr=&((struct sockaddr_in *)ifa->ifa_addr)->sin_addr;
            log_debug(nnti_debug_level, "checking iface (IPv6) name (%s)", ifa->ifa_name);
            inet_ntop(AF_INET6, tmpAddrPtr, ipaddr, maxlen);
            log_debug(nnti_debug_level, "hostname(%s) has IP Address %s", ifa->ifa_name, ipaddr);
            if (0==strncmp(ifa->ifa_name, LISTEN_IFACE_BASENAME, strlen(LISTEN_IFACE_BASENAME))) {
                log_debug(nnti_debug_level, "hostname(%s) matches", ifa->ifa_name);
                break;
            }
        }
    }
    if (ifAddrStruct!=NULL) freeifaddrs(ifAddrStruct);

    log_debug(nnti_debug_level, "ipaddr(%s)", ipaddr);
    if (ipaddr[0]=='\0') {
        return NNTI_ENOENT;
    } else {
        return NNTI_OK;
    }
}


static void config_init(nnti_gni_config_t *c)
{
    c->min_atomics_vars       =512;
    c->use_alps_ptag          =true;
    c->use_wr_pool            =false;
    c->wr_pool_initial_size   =50;
    c->wr_pool_max_size       =1000;
    c->wr_pool_create_if_empty=false;
    c->use_rdma_target_ack    =false;
    c->use_rdma_events        =false;
    c->use_rdma_fence         =false;
    c->pi_ordering            =PI_ORDERING_DEFAULT; /* DEFAULT, STRICT, RELAXED */
    c->rdma_mode              =RDMA_MIXED;          /* FMA, BTE, MIXED, CROSSOVER */
    c->fma_bte_crossover_size =4096;
    c->max_timeout_ms         =-1;
}

static void config_get_from_env(nnti_gni_config_t *c)
{
    char *env_str=NULL;

    if ((env_str=getenv("TRIOS_NNTI_USE_ALPS_PTAG")) != NULL) {
        if ((!strcasecmp(env_str, "TRUE")) ||
            (!strcmp(env_str, "1"))) {
            log_debug(nnti_debug_level, "setting c->use_alps_ptag to TRUE");
            c->use_alps_ptag=true;
        } else {
            log_debug(nnti_debug_level, "setting c->use_alps_ptag to FALSE");
            c->use_alps_ptag=false;
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_USE_ALPS_PTAG is undefined.  using c->use_alps_ptag default");
    }
    if ((env_str=getenv("TRIOS_NNTI_USE_WR_POOL")) != NULL) {
        if ((!strcasecmp(env_str, "TRUE")) ||
            (!strcmp(env_str, "1"))) {
            log_debug(nnti_debug_level, "setting c->use_wr_pool to TRUE");
            c->use_wr_pool=true;
        } else {
            log_debug(nnti_debug_level, "setting c->use_wr_pool to FALSE");
            c->use_wr_pool=false;
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_USE_WR_POOL is undefined.  using c->use_wr_pool default");
    }
    if ((env_str=getenv("TRIOS_NNTI_WR_POOL_INITIAL_SIZE")) != NULL) {
        errno=0;
        uint32_t initial_size=strtoul(env_str, NULL, 0);
        if (errno == 0) {
            log_debug(nnti_debug_level, "setting c->wr_pool_initial_size to %lu", initial_size);
            c->wr_pool_initial_size=initial_size;
        } else {
            log_debug(nnti_debug_level, "TRIOS_NNTI_WR_POOL_INITIAL_SIZE value conversion failed (%s).  using c->wr_pool_initial_size default.", strerror(errno));
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_WR_POOL_INITIAL_SIZE is undefined.  using c->wr_pool_initial_size default");
    }
    if ((env_str=getenv("TRIOS_NNTI_WR_POOL_MAX_SIZE")) != NULL) {
        errno=0;
        uint32_t max_size=strtoul(env_str, NULL, 0);
        if (errno == 0) {
            log_debug(nnti_debug_level, "setting c->wr_pool_max_size to %lu", max_size);
            c->wr_pool_max_size=max_size;
        } else {
            log_debug(nnti_debug_level, "TRIOS_NNTI_WR_POOL_MAX_SIZE value conversion failed (%s).  using c->wr_pool_max_size default.", strerror(errno));
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_WR_POOL_MAX_SIZE is undefined.  using c->wr_pool_max_size default");
    }
    if ((env_str=getenv("TRIOS_NNTI_WR_POOL_CREATE_IF_EMPTY")) != NULL) {
        if ((!strcasecmp(env_str, "TRUE")) ||
            (!strcmp(env_str, "1"))) {
            log_debug(nnti_debug_level, "setting c->wr_pool_create_if_empty to TRUE");
            c->wr_pool_create_if_empty=true;
        } else {
            log_debug(nnti_debug_level, "setting c->wr_pool_create_if_empty to FALSE");
            c->wr_pool_create_if_empty=false;
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_WR_POOL_CREATE_IF_EMPTY is undefined.  using c->wr_pool_create_if_empty default");
    }
    if ((env_str=getenv("TRIOS_NNTI_USE_RDMA_TARGET_ACK")) != NULL) {
        if ((!strcasecmp(env_str, "TRUE")) ||
            (!strcmp(env_str, "1"))) {
            log_debug(nnti_debug_level, "setting c->use_rdma_target_ack to TRUE");
            c->use_rdma_target_ack=true;
        } else {
            log_debug(nnti_debug_level, "setting c->use_rdma_target_ack to FALSE");
            c->use_rdma_target_ack=false;
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_USE_RDMA_TARGET_ACK is undefined.  using c->use_rdma_target_ack default");
    }
    if ((env_str=getenv("TRIOS_NNTI_USE_RDMA_EVENTS")) != NULL) {
        if ((!strcasecmp(env_str, "TRUE")) ||
            (!strcmp(env_str, "1"))) {
            log_debug(nnti_debug_level, "setting c->use_rdma_events to TRUE");
            c->use_rdma_events=true;
        } else {
            log_debug(nnti_debug_level, "setting c->use_rdma_events to FALSE");
            c->use_rdma_events=false;
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_USE_RDMA_EVENTS is undefined.  using c->use_rdma_events default");
    }
    if ((env_str=getenv("TRIOS_NNTI_USE_RDMA_FENCE")) != NULL) {
        if ((!strcasecmp(env_str, "TRUE")) ||
            (!strcmp(env_str, "1"))) {
            log_debug(nnti_debug_level, "setting c->use_rdma_fence to TRUE");
            c->use_rdma_fence=true;
        } else {
            log_debug(nnti_debug_level, "setting c->use_rdma_fence to FALSE");
            c->use_rdma_fence=false;
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_USE_RDMA_FENCE is undefined.  using c->use_rdma_fence default");
    }
    if ((env_str=getenv("TRIOS_NNTI_PI_ORDERING")) != NULL) {
        if (!strcasecmp(env_str, "PI_ORDERING_DEFAULT")) {
            log_debug(nnti_debug_level, "setting c->pi_ordering to PI_ORDERING_DEFAULT");
            c->pi_ordering=PI_ORDERING_DEFAULT;
        } else if (!strcasecmp(env_str, "PI_ORDERING_STRICT")) {
            log_debug(nnti_debug_level, "setting c->pi_ordering to PI_ORDERING_STRICT");
            c->pi_ordering=PI_ORDERING_STRICT;
        } else if (!strcasecmp(env_str, "PI_ORDERING_RELAXED")) {
            log_debug(nnti_debug_level, "setting c->pi_ordering to PI_ORDERING_RELAXED");
            c->pi_ordering=PI_ORDERING_RELAXED;
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_PI_ORDERING is undefined.  using c->pi_ordering default");
    }
    if ((env_str=getenv("TRIOS_NNTI_RDMA_MODE")) != NULL) {
        if (!strcasecmp(env_str, "RDMA_FMA")) {
            log_debug(nnti_debug_level, "setting c->rdma_mode to RDMA_FMA");
            c->rdma_mode=RDMA_FMA;
        } else if (!strcasecmp(env_str, "RDMA_BTE")) {
            log_debug(nnti_debug_level, "setting c->rdma_mode to RDMA_BTE");
            c->rdma_mode=RDMA_BTE;
        } else if (!strcasecmp(env_str, "RDMA_MIXED")) {
            log_debug(nnti_debug_level, "setting c->rdma_mode to RDMA_MIXED");
            c->rdma_mode=RDMA_MIXED;
        } else if (!strcasecmp(env_str, "RDMA_CROSSOVER")) {
            log_debug(nnti_debug_level, "setting c->rdma_mode to RDMA_CROSSOVER");
            c->rdma_mode=RDMA_CROSSOVER;
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_RDMA_MODE is undefined.  using c->rdma_mode default");
    }
    if ((env_str=getenv("TRIOS_NNTI_FMA_BTE_CROSSOVER_SIZE")) != NULL) {
        errno=0;
        uint32_t crossover_size=strtoul(env_str, NULL, 0);
        if (errno == 0) {
            log_debug(nnti_debug_level, "setting c->fma_bte_crossover_size to %lu", crossover_size);
            c->fma_bte_crossover_size=crossover_size;
        } else {
            log_debug(nnti_debug_level, "TRIOS_NNTI_FMA_BTE_CROSSOVER_SIZE value conversion failed (%s).  using c->fma_bte_crossover_size default.", strerror(errno));
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_FMA_BTE_CROSSOVER_SIZE is undefined.  using c->fma_bte_crossover_size default");
    }
    if ((env_str=getenv("TRIOS_NNTI_MAX_TIMEOUT_MS")) != NULL) {
        errno=0;
        int32_t max_timeout_ms=strtoul(env_str, NULL, 0);
        if (errno == 0) {
            log_debug(nnti_debug_level, "setting c->max_timeout_ms to %lu", max_timeout_ms);
            c->max_timeout_ms=max_timeout_ms;
        } else {
            log_debug(nnti_debug_level, "TRIOS_NNTI_MAX_TIMEOUT_MS value conversion failed (%s).  using c->max_timeout_ms default.", strerror(errno));
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_MAX_TIMEOUT_MS is undefined.  using c->max_timeout_ms default");
    }
    if ((env_str=getenv("TRIOS_NNTI_MIN_ATOMIC_VARS")) != NULL) {
        errno=0;
        uint32_t min_vars=strtoul(env_str, NULL, 0);
        if (errno == 0) {
            log_debug(nnti_debug_level, "setting c->min_atomics_vars to %lu", min_vars);
            c->min_atomics_vars=min_vars;
        } else {
            log_debug(nnti_debug_level, "TRIOS_NNTI_MIN_ATOMIC_VARS value conversion failed (%s).  using c->min_atomics_vars default.", strerror(errno));
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_MIN_ATOMIC_VARS is undefined.  using c->min_atomics_vars default");
    }
}

inline int DEQUEUE_POST_DESCRIPTOR(
        gni_cq_handle_t           cq_hdl,
        gni_cq_entry_t           *ev_data,
        gni_post_descriptor_t   **post_desc_ptr)
{
    int rc=GNI_RC_SUCCESS;

    *post_desc_ptr=NULL;

    log_debug(nnti_debug_level, "calling GetComplete(DEQUEUE_POST_DESCRIPTOR)");
    rc=GNI_GetCompleted_wrapper (cq_hdl, *ev_data, post_desc_ptr);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "GetCompleted(DEQUEUE_POST_DESCRIPTOR post_desc_ptr(%p)) failed: %d", *post_desc_ptr, rc);
    print_post_desc(*post_desc_ptr);

    return(rc);
}

inline int IS_RDMA_INITIATOR_BUFFER(
        const nnti_gni_work_request_t *wr)
{
    int rc=FALSE;
    nnti_gni_memory_handle_t *hdl=(nnti_gni_memory_handle_t *)wr->reg_buf->transport_private;

    switch (hdl->type) {
        case REQUEST_BUFFER:
        case RECEIVE_BUFFER:
        case SEND_BUFFER:
        case GET_SRC_BUFFER:
        case PUT_DST_BUFFER:
            rc=FALSE;
            break;
        case GET_DST_BUFFER:
        case PUT_SRC_BUFFER:
            rc=TRUE;
            break;
        case RDMA_TARGET_BUFFER:
            if ((wr->last_op==GNI_OP_GET_INITIATOR) ||
                (wr->last_op==GNI_OP_PUT_INITIATOR)) {
                rc=TRUE;
            } else {
                rc=FALSE;
            }
            break;
        case UNKNOWN_BUFFER:
        default:
            log_error(nnti_debug_level, "Unknown buffer type: %d", hdl->type);
            rc=FALSE;
            break;
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);
    return(rc);
}

inline int IS_RDMA_TARGET_BUFFER(
        const nnti_gni_work_request_t *wr)
{
    int rc=FALSE;
    nnti_gni_memory_handle_t *hdl=(nnti_gni_memory_handle_t *)wr->reg_buf->transport_private;

    switch (hdl->type) {
        case REQUEST_BUFFER:
        case RECEIVE_BUFFER:
        case SEND_BUFFER:
        case GET_DST_BUFFER:
        case PUT_SRC_BUFFER:
            rc=FALSE;
            break;
        case GET_SRC_BUFFER:
        case PUT_DST_BUFFER:
            rc=TRUE;
            break;
        case RDMA_TARGET_BUFFER:
            if ((wr->last_op==GNI_OP_GET_TARGET) ||
                (wr->last_op==GNI_OP_PUT_TARGET)) {
                rc=TRUE;
            } else {
                rc=FALSE;
            }
            break;
        case UNKNOWN_BUFFER:
        default:
            log_error(nnti_debug_level, "Unknown buffer type: %d", hdl->type);
            rc=FALSE;
            break;
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);
    return(rc);
}
