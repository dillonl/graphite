// Boost.Process
//
// Copyright (c) 2006 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

#ifndef BOOST_PROCESSES_POSIX_LAUNCH_HPP_122922
#define BOOST_PROCESSES_POSIX_LAUNCH_HPP_122922

#include "boost/processes/launch.hpp"
#include "boost/processes/detail/error.hpp"
#include "boost/processes/detail/pipe.hpp"
#include "boost/processes/detail/posix_util.hpp"
#include "boost/processes/stream_behavior.hpp"
#include "boost/assert.hpp"
#include "boost/scoped_array.hpp"
#include "boost/variant/apply_visitor.hpp"
#include "boost/variant/get.hpp"
#include "boost/variant/static_visitor.hpp"
#include <fcntl.h>
#include <unistd.h>
#include <cstdlib>
#include <cstring>
#include <map>
#include <utility>

namespace boost {
namespace processes {
namespace detail {

class setup_stream: public boost::static_visitor<>
{
    child::pe_map& pes_;

public:
    setup_stream(child::pe_map& pes)
    : pes_(pes)
    {}

    void operator()(const inherit_stream&) const
    {
        // do nothing
    }

    void operator()(const capture_stream& c) const
    {
        const fileno fn = c.fn();
        child::pe_map::iterator it = pes_.find(fn);
        BOOST_ASSERT(it != pes_.end());
        pipe_end& pe = it->second;
        BOOST_ASSERT(pe.native() != fn);
        posix_dup2(pe.native(), fn);
        pe.close();
    }

    void operator()(const redirect_to_stdout& r) const
    {
        posix_dup2(STDOUT_FILENO, r.fn());
    }

    void operator()(const silence_stream& s) const
    {
        if (s.mode() == silence_stream::in)
            (*this)(redirect_to_named_file(s.fn(), redirect_to_named_file::in,
                    "/dev/zero"));
        else
            (*this)(redirect_to_named_file(s.fn(), redirect_to_named_file::out,
                    "/dev/null"));
    }

    void operator()(const close_stream& c) const
    {
        ::close(c.fn());
    }

    void operator()(const redirect_to_named_file& r) const
    {
        int fd = ::open(r.filename().c_str(),
                        r.mode() == redirect_to_named_file::in
                            ? O_RDONLY: O_WRONLY);
        if (fd == -1)
            BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("open(2) failed");
        if (fd != r.fn()) {
            handle h(fd);
            posix_dup2(h.native(), r.fn());
        }
    }

    void operator()(const redirect_to_native_file& r) const
    {
        posix_dup2(r.native_file(), r.fn());
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

    typedef std::vector<std::string> arguments;
    using namespace detail;

    // occupy fd's before creating pipes
    handle tmp_handle(pipe().rend());
    typedef std::vector<handle> handle_vector;
    handle_vector tmp_handles;
    tmp_handles.reserve(ctx.streams_behavior.size());
    for (streams_behavior::const_iterator it = ctx.streams_behavior.begin(),
         end = ctx.streams_behavior.end(); it != end; ++it) {
        if (!boost::get<inherit_stream>(&*it)
           && !boost::get<close_stream>(&*it)) {
            const fileno fn = get_fileno(*it);
            if (fn != tmp_handle.native() && ::fcntl(fn, F_GETFD, 0) < 0) {
                if (errno != EBADF)
                    BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("fcntl(2) failed");
                posix_dup2(tmp_handle.native(), fn);
                handle h(fn);
                tmp_handles.push_back(h);
            }
        }
    }

    // create pipes
    child::pe_map parent_pes, child_pes;
    for (streams_behavior::const_iterator it = ctx.streams_behavior.begin(),
         end = ctx.streams_behavior.end(); it != end; ++it) {
        if (const capture_stream* c = boost::get<capture_stream>(&*it)) {
            fileno fn = c->fn();
            pipe p;
            child_pes.insert(std::make_pair(
                fn, c->mode() == capture_stream::in ? p.rend() : p.wend()));
            pipe_end& pe =
                c->mode() == capture_stream::in ? p.wend() : p.rend();
            posix_close_on_exec(pe.native(), true);
            parent_pes.insert(std::make_pair(fn, pe));
        }
    }

    // close unnecessary handles
    tmp_handle.close();
    for (handle_vector::iterator it = tmp_handles.begin(),
        end = tmp_handles.end(); it != end; ++it) {
            it->close();
    }

    pid_t pid = ::fork();
    if (pid == -1) {
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("fork(2) failed");
    } else if (pid == 0) {
        std::pair<std::size_t, char**> argv;
        char** envp;

        try {
            // close unnecessary handles
            for (child::pe_map::iterator it = parent_pes.begin(),
                 end = parent_pes.end(); it != end; ++it) {
                it->second.close();
            }

            for (streams_behavior::const_iterator it =
                 ctx.streams_behavior.begin(), end = ctx.streams_behavior.end();
                 it != end; ++it) {
                boost::apply_visitor(setup_stream(child_pes), *it);
            }

            ctx.apply();

            argv = vector_to_argv(args);
            envp = environment_to_envp(ctx.environment);

        } catch (const boost::system::system_error& e) {
            ::write(STDERR_FILENO, e.what(), std::strlen(e.what()));
            ::write(STDERR_FILENO, "\n", 1);
            ::_exit(127);
        }

        ::execve(exe.c_str(), argv.second, envp);

        for (arguments::size_type i = 0; i < argv.first; i++)
            delete [] argv.second[i];
        delete [] argv.second;

        for (arguments::size_type i = 0; i < ctx.environment.size(); i++)
            delete [] envp[i];
        delete [] envp;

        boost::system::system_error e(
            boost::system::error_code(errno,
                                      boost::system::get_system_category()),
            BOOST_PROCESSES_SOURCE_LOCATION "execve(2) failed");

        ::write(STDERR_FILENO, e.what(), std::strlen(e.what()));
        ::write(STDERR_FILENO, "\n", 1);
        ::_exit(127);
    }

    // close unnecessary handles
    for (child::pe_map::iterator it = child_pes.begin(),
        end = child_pes.end(); it != end; ++it) {
            it->second.close();
    }

    return child(pid, parent_pes);
}

} // namespace processes
} // namespace boost

#endif
