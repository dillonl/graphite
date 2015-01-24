// Boost.Process
//
// Copyright (c) 2006 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

#ifndef BOOST_PROCESSES_WINDOWS_LAUNCH_HPP_122945
#define BOOST_PROCESSES_WINDOWS_LAUNCH_HPP_122945

#include "boost/processes/launch.hpp"
#include "boost/processes/detail/error.hpp"
#include "boost/processes/detail/handle.hpp"
#include "boost/processes/detail/pipe.hpp"
#include "boost/processes/detail/windows_api.hpp"
#include "boost/processes/detail/windows_util.hpp"
#include "boost/processes/environment.hpp"
#include "boost/processes/pipe_end.hpp"
#include "boost/processes/stream_behavior.hpp"
#include "boost/assert.hpp"
#include "boost/scoped_array.hpp"
#include "boost/variant/apply_visitor.hpp"
#include "boost/variant/static_visitor.hpp"
#include <cstring>

namespace boost {
namespace processes {
namespace detail {

class setup_stream: public boost::static_visitor<windows_api::HANDLE>
{
    handle& h_;
    pipe_end& pe_;
    const windows_api::HANDLE sout_;
    const windows_api::DWORD n_std_handle_;

public:
    setup_stream(handle& h,
                 pipe_end& pe,
                 windows_api::HANDLE sout,
                 windows_api::DWORD n_std_handle)
    : h_(h)
    , pe_(pe)
    , sout_(sout)
    , n_std_handle_(n_std_handle)
    {}

    windows_api::HANDLE
    operator()(const inherit_stream&) const
    {
        using namespace windows_api;
        HANDLE h = GetStdHandle(n_std_handle_);
        if (h == INVALID_HANDLE_VALUE_)
            BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("GetStdHandle() failed");
        // standard handle may be overwritten by SetStdHandle() with some
        // non-inheritable handle, so we duplicate it here
        h_ = windows_dup(h, true);
        return h_.native();
    }

    windows_api::HANDLE
    operator()(const capture_stream& c) const
    {
        pipe p;
        if (c.mode() == capture_stream::in) {
            h_ = p.rend();
            pe_ = p.wend();
        } else {
            h_ = p.wend();
            pe_ = p.rend();
        }
        windows_set_inheritable(h_.native(), true);
        return h_.native();
    }

    windows_api::HANDLE
    operator()(const redirect_to_stdout&) const
    {
        BOOST_ASSERT(sout_ != windows_api::INVALID_HANDLE_VALUE_);
        h_ = windows_dup(sout_, true);
        return h_.native();
    }

    windows_api::HANDLE
    operator()(const silence_stream& s) const
    {
        return (*this)(redirect_to_named_file(s.fn(),
                                              s.mode() == silence_stream::in
                                                  ? redirect_to_named_file::in
                                                  : redirect_to_named_file::out,
                                              "NUL"));
    }

    windows_api::HANDLE
    operator()(const close_stream&) const
    {
        // do nothing
        return windows_api::INVALID_HANDLE_VALUE_;
    }

    windows_api::HANDLE
    operator()(const redirect_to_named_file& r) const
    {
        using namespace windows_api;
        HANDLE h = CreateFileA(r.filename().c_str(),
                               r.mode() == redirect_to_named_file::in
                                   ? GENERIC_READ_ : GENERIC_WRITE_,
                               0, 0,
                               r.mode() == redirect_to_named_file::in
                                    ? OPEN_EXISTING_ : CREATE_NEW_,
                               FILE_ATTRIBUTE_NORMAL_, 0);
        if (h == INVALID_HANDLE_VALUE_)
            BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("CreateFile() failed");
        h_ = handle(h);
        return h_.native();
    }

    windows_api::HANDLE
    operator()(const redirect_to_native_file& r) const
    {
        windows_set_inheritable(r.native_file(), true);
        return r.native_file();
    }
};

} // namespace detail

BOOST_PROCESSES_INLINE_IF_HEADER_ONLY child
launch(const std::string& exe,
       const std::vector<std::string>& args,
       const context& ctx)
{
    BOOST_ASSERT(!args.empty());
    // Validate execution context.
    // XXX Should this be a 'validate()' method in it?
    BOOST_ASSERT(!ctx.work_directory.empty());

    using namespace detail::windows_api;
    using namespace detail;

    STARTUPINFOA* si;
    STARTUPINFOA si1 = {0};
    boost::scoped_array<char> si2;
    if (ctx.startupinfo) {
        si = static_cast<STARTUPINFOA*>(ctx.startupinfo);
        BOOST_ASSERT(si->cb >= sizeof(STARTUPINFOA));
        si2.reset(new char[si->cb]);
        std::memcpy(si2.get(), si, si->cb);
        si = static_cast<STARTUPINFOA*>(static_cast<void*>(si2.get()));
    } else {
        si = &si1;
        si->cb = sizeof(si1);
    }
    si->dwFlags |= STARTF_USESTDHANDLES_;

    handle sin_handle, sout_handle, serr_handle;
    pipe_end sin_pe, sout_pe, serr_pe;

    for (streams_behavior::const_iterator it = ctx.streams_behavior.begin(),
         end = ctx.streams_behavior.end(); it != end; ++it) {
        switch (get_fileno(*it)) {
        case stdin_fileno:
            si->hStdInput = boost::apply_visitor(
                setup_stream(
                  sin_handle, sin_pe, si->hStdOutput, STD_INPUT_HANDLE_),
                *it);
            break;
        case stdout_fileno:
            si->hStdOutput = boost::apply_visitor(
                setup_stream(
                  sout_handle, sout_pe,  si->hStdOutput, STD_OUTPUT_HANDLE_),
                *it);
            break;
        case stderr_fileno:
            si->hStdError = boost::apply_visitor(
                setup_stream(
                  serr_handle, serr_pe, si->hStdOutput, STD_ERROR_HANDLE_),
                *it);
            break;
        default:
            BOOST_ASSERT(false);
        }
    }

    boost::shared_array<char> cmdline(vector_to_windows_cmdline(args));
    boost::shared_array<char> envstrs(
        environment_to_windows_string(ctx.environment));
    PROCESS_INFORMATION pi;
	
    if (!CreateProcessA(exe.c_str(), cmdline.get(), 0, 0, true, 0,
                        envstrs.get(), ctx.work_directory.c_str(), si, &pi)) {
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("CreateProcess() failed");
    }
    CloseHandle(pi.hThread); // FIXME: check for errors
    return child(pi.hProcess, sin_pe, sout_pe, serr_pe);
}

} // namespace processes
} // namespace boost

#endif
