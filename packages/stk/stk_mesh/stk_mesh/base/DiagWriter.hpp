#ifndef stk_mesh_DiagWriter_h
#define stk_mesh_DiagWriter_h

#include <sstream>

#include <stk_mesh/base/Types.hpp>      // for EntityProc
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey

#ifdef STK_MESH_TRACE_ENABLED
#include <stk_mesh/base/DiagWriter_fwd.hpp>
#include <stk_util/util/Writer.hpp>
#include <stk_util/util/Trace.hpp>
#include <stk_util/util/WriterParser.hpp>
#endif

namespace stk { namespace diag { class Writer; } }
namespace stk { namespace mesh { class Part; } }


// Note, this classes/functions in this header are for internal use only.
// The API for tracing is defined in Trace.hpp

namespace stk {
namespace mesh {

namespace impl {
class Partition;
}

#ifdef STK_MESH_TRACE_ENABLED

// Must be called before theDiagWriter/meshlog
void initDiagWriter(std::ostream& stream);

stk::diag::Writer &theDiagWriter();
#define meshlog stk::mesh::theDiagWriter()

class DiagWriterParser : public diag::WriterParser
{
public:
  DiagWriterParser();
};

DiagWriterParser &theDiagWriterParser();

typedef diag::Tracespec Tracespec;
typedef diag::Traceback Traceback;

/**
 * Defines the Trace objects that will be used in stk... this class
 * mostly just ensures that the correct DiagWriter is used for
 * tracing.
 */
class Trace : public diag::Trace
{
public:
  explicit Trace(const char *message)
    : diag::Trace(meshlog, message)
  {}

  Trace(const char *message, int print_mask)
    : diag::Trace(meshlog, message, print_mask)
  {}

  Trace(const char *message, int print_mask, bool do_trace)
    : diag::Trace(meshlog, message, print_mask, do_trace)
  {}
};

// If Writer does not know how to output an object you want to trace, you
// can address that here by defining an operator<< for that object. Note
// that Writer handles vectors and pointers automatically.

stk::diag::Writer& operator<<(stk::diag::Writer& writer, const Part& part);

stk::diag::Writer& operator<<(stk::diag::Writer& writer, const EntityKey& key);

stk::diag::Writer& operator<<(stk::diag::Writer& writer, const EntityProc& entity_proc);

stk::diag::Writer& operator<<(stk::diag::Writer& writer, const Bucket& bucket);

namespace impl {
stk::diag::Writer& operator<<(stk::diag::Writer& writer, const Partition& partition);
}

#endif // STKMESH_TRACE_ENABLED

} // namespace mesh
} // namespace stk

#endif // stk_mesh_DiagWriter_h
