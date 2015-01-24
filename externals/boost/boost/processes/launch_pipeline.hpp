// Boost.Process
//
// Copyright (c) 2006 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

/// \file boost/process/launch_pipeline.hpp
///
/// Provides launch_pipeline() function.

#ifndef BOOST_PROCESSES_LAUNCH_PIPELINE_HPP_175117
#define BOOST_PROCESSES_LAUNCH_PIPELINE_HPP_175117

#include "boost/processes/config.hpp"
#include "boost/processes/child.hpp"
#include "boost/processes/context.hpp"
#include "boost/processes/detail/error.hpp"
#include "boost/processes/detail/pipe.hpp"
#include "boost/processes/launch.hpp"
#include "boost/assert.hpp"
#include "boost/scoped_array.hpp"
#include "boost/variant/get.hpp"

#ifdef BOOST_POSIX_API
#include "boost/processes/detail/posix_util.hpp"
#endif

namespace boost {
namespace processes {

/// \brief Launches a pipelined set of child processes.
///
/// Given a collection of pipeline_entry objects describing how to launch
/// a set of child processes, spawns them all and connects their inputs and
/// outputs in a way that permits pipelined communication.
///
/// \pre Let 1..N be the processes in the collection: the input behavior of
///      the 2..N processes must be set to close_stream().
/// \pre Let 1..N be the processes in the collection: the output behavior of
///      the 1..N-1 processes must be set to close_stream().
///
/// \remark <b>Blocking remarks</b>: This function may block if the
///         device holding the executable of one of the entries
///         blocks when loading the image. This might happen if, e.g.,
///         the binary is being loaded from a network share.
///
/// \return A set of Child objects that represent all the processes spawned
///         by this call. You should use wait_children() to wait for their
///         termination.
template<class Entries>
children
launch_pipeline(const Entries& entries)
{
    BOOST_ASSERT(entries.size() >= 2);
//  BOOST_ASSERT(boost::get<close_stream>(&entries[0].context.stdout_behavior));

    // Method's result value.
    children cs;

    // The pipes used to connect the pipeline's internal process.
    boost::scoped_array<detail::pipe> pipes
        (new detail::pipe[entries.size() - 1]);

#if defined(BOOST_POSIX_API)

    for (std::size_t i = 0; i < entries.size() - 1; ++i) {
        detail::posix_close_on_exec(pipes[i].rend().native(), true);
        detail::posix_close_on_exec(pipes[i].wend().native(), true);
    }

    // Configure and spawn the pipeline's first process.
    {
        typename Entries::value_type::context_type ctx(entries[0].context);

        pipe_end& we = pipes[0].wend();
        ctx.streams_behavior.push_front(
            redirect_to_native_file(stdout_fileno, we.native()));

        child ch(launch(entries[0].executable, entries[0].arguments, ctx));

        we.close();

        cs.push_back(ch);
    }

    // Configure and spawn the pipeline's internal processes.
    for (typename Entries::size_type i = 1; i < entries.size() - 1; i++) {
        typename Entries::value_type::context_type ctx(entries[i].context);

//        BOOST_ASSERT(boost::get<close_stream>(&ctx.stdin_behavior));
        pipe_end& re = pipes[i - 1].rend();
        ctx.streams_behavior.push_front(
            redirect_to_native_file(stdin_fileno, re.native()));
//        BOOST_ASSERT(boost::get<close_stream>(&ctx.stdout_behavior));
        pipe_end& we = pipes[i].wend();
        ctx.streams_behavior.push_front(
            redirect_to_native_file(stdout_fileno, we.native()));

        child ch(launch(entries[i].executable, entries[i].arguments, ctx));

        re.close();
        we.close();

        cs.push_back(ch);
    }

    // Configure and spawn the pipeline's last process.
    {
        const typename Entries::size_type last = entries.size() - 1;
        typename Entries::value_type::context_type ctx(entries[last].context);

//        BOOST_ASSERT(boost::get<close_stream>(&ctx.stdin_behavior));
        pipe_end& re = pipes[last - 1].rend();
        ctx.streams_behavior.push_front(
            redirect_to_native_file(stdin_fileno, re.native()));

        child ch(launch(entries[last].executable,
                        entries[last].arguments,
                        ctx));

        re.close();

        cs.push_back(ch);
    }
#elif defined(BOOST_WINDOWS_API)
    // Configure and spawn the pipeline's first process.
    {
        typename Entries::value_type::context_type ctx(entries[0].context);

        pipe_end& we = pipes[0].wend();
        ctx.streams_behavior.push_front(
            redirect_to_native_file(stdout_fileno, we.native()));

        child ch(launch(entries[0].executable, entries[0].arguments, ctx));

        we.close();

        cs.push_back(ch);
    }

    // Configure and spawn the pipeline's internal processes.
    for (typename Entries::size_type i = 1; i < entries.size() - 1; i++) {
        typename Entries::value_type::context_type ctx(entries[i].context);

//        BOOST_ASSERT(boost::get<close_stream>(&ctx.stdin_behavior));
        pipe_end& re = pipes[i - 1].rend();
        ctx.streams_behavior.push_front(
            redirect_to_native_file(stdin_fileno, re.native()));

//        BOOST_ASSERT(boost::get<close_stream>(&ctx.stdout_behavior));
        pipe_end& we = pipes[i].wend();
        ctx.streams_behavior.push_front(
            redirect_to_native_file(stdout_fileno, we.native()));

        child ch(launch(entries[i].executable, entries[i].arguments, ctx));

        re.close();
        we.close();

        cs.push_back(ch);
    }

    // Configure and spawn the pipeline's last process.
    {
        const typename Entries::size_type last = entries.size() - 1;
        typename Entries::value_type::context_type ctx(entries[last].context);

//        BOOST_ASSERT(boost::get<close_stream>(&ctx.stdin_behavior));
        pipe_end re = pipes[last - 1].rend();
        ctx.streams_behavior.push_front(
            redirect_to_native_file(stdin_fileno, re.native()));

        child ch(launch(entries[last].executable,
                        entries[last].arguments,
                        ctx));

        re.close();

        cs.push_back(ch);
    }
#endif

    return cs;
}

/// \brief Waits for a collection of children to terminate.
///
/// Given a collection of Child objects (such as std::vector<child> or
/// the convenience children type), waits for the termination of all of
/// them.
///
/// \remark <b>Blocking remarks</b>: This call blocks if any of the
/// children processes in the collection has not finalized execution and
/// waits until it terminates.
///
/// \return The exit status of the first process that returns an error
///         code or, if all of them executed correctly, the exit status
///         of the last process in the collection.
template<class Children>
const status
wait_children(Children& cs)
{
    BOOST_ASSERT(cs.size() >= 2);

    typename Children::iterator iter = cs.begin();
    while (iter != cs.end()) {
        const status s = iter->wait();
        iter++;
        if (iter == cs.end())
            return s;
        else if (!s.exited() || s.exit_status() != EXIT_SUCCESS) {
            while (iter != cs.end()) {
                iter->wait();
                iter++;
            }
            return s;
        }
    }

    BOOST_ASSERT(false);
    return cs.begin()->wait();
}

} // namespace processes
} // namespace boost

#endif
