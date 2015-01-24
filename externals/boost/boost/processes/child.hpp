// Boost.Process
//
// Copyright (c) 2006, 2007 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

/// \file boost/process/child.hpp
///
/// Includes the declaration of the child class.

#ifndef BOOST_PROCESSES_CHILD_HPP_171406
#define BOOST_PROCESSES_CHILD_HPP_171406

#include "boost/processes/detail/pipe.hpp"
#include "boost/processes/pipe_end.hpp"
#include "boost/processes/status.hpp"
#include "boost/assert.hpp"
#include "boost/shared_ptr.hpp"
#if defined(BOOST_POSIX_API)
    #include <sys/types.h>
    #include <sys/wait.h>
	#include <signal.h>	//added by rtao, testing!!!
	#include <unistd.h>	//added by rtao, testing!!!
#elif defined(BOOST_WINDOWS_API)
    #include "boost/processes/detail/windows_api.hpp"
#endif
#include <map>
#include <vector>

//added by rtao, for kill
#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>
#include <stdio.h> //for remove( ) 
//added end

#ifdef BOOST_MSVC
#pragma warning(push)
// 4251 * needs to have dll-interface to be used by clients of *
// 4275 * used as base for dll-interface class *
#pragma warning(disable: 4251 4275)
#endif

namespace boost {
namespace processes {

/// \brief Generic implementation of the Child concept.
///
/// The child class implements the Child concept in an operating system
/// agnostic way.
///
class BOOST_PROCESSES_DECL child
{
public:
#if defined(BOOST_PROCESS_DOXYGEN)
    /// \brief The native representaion of a child.
    typedef implementation_defined native_type;
#elif defined(BOOST_WINDOWS_API)
    typedef detail::windows_api::HANDLE native_type;
#elif defined(BOOST_POSIX_API)
    typedef pid_t native_type;
#endif

    /// \brief Constructs a new child object representing a just spawned
    ///        child process.
    ///
    /// Creates a new child object that represents the just spawned process
    /// \a h.
    ///
    /// The \a chin, \a chout and \a cherr file handles represent
    /// the parent's handles used to communicate with the corresponding
    /// data streams. They needn't be valid but their availability must
    /// match the redirections configured by the launcher that spawned this
    /// process.
    ///
    child(native_type native,
          const pipe_end& chin,
          const pipe_end& chout,
          const pipe_end& cherr)
    :
#if defined(BOOST_WINDOWS_API)
      handle_(native, detail::windows_api::CloseHandle)
#elif defined(BOOST_POSIX_API)
      pid_(native)
#endif
    , chin_(chin)
    , chout_(chout)
    , cherr_(cherr)
    {}

    // Default constructor
	child() {}	/*This is added by Jianrong and rtao*/

    native_type native()
    {
        return
#if defined(BOOST_WINDOWS_API)
            handle_.get();
#elif defined(BOOST_POSIX_API)
            pid_;
#endif
    }

    /// \brief Gets a reference to the child's standard input stream.
    ///
    /// Returns a reference to an object that represents the
    /// standard input communication channel with the child process.
    ///
    pipe_end get_stdin() const
    {
        return chin_;
    }

    /// \brief Gets a reference to the child's standard output stream.
    ///
    /// Returns a reference to an object that represents the
    /// standard output communication channel with the child process.
    ///
    pipe_end get_stdout() const
    {
        return chout_;
    }

    /// \brief Gets a reference to the child's standard error stream.
    ///
    /// Returns a reference to an object that represents the
    /// standard error communication channel with the child process.
    ///
    pipe_end get_stderr() const
    {
        return cherr_;
    }

    /// \brief Terminates the process execution.
    ///
    /// Forces the termination of the process execution. Some platforms
    /// allow processes to ignore some external termination notifications
    /// or to capture them for a proper exit cleanup. You can set the
    /// \a force flag to true in them to force their termination regardless
    /// of any exit handler.
    ///
    /// After this call, accessing this object can be dangerous because the
    /// process identifier may have been reused by a different process. It
    /// might still be valid, though, if the process has refused to die.
    ///
    /// \throw system_error If the system call used to terminate the
    ///                     process fails.
    void terminate(bool force = false) const;

	//This part is added to boost::processes by rtao for terminate the children
	//processes tree of this process.
	//bool terminate_tree(std::string chilldren_info_file) const;
	bool terminate_process_tree() const;


    /// \brief Blocks and waits for the child process to terminate.
    ///
    /// Returns a status object that represents the child process'
    /// finalization condition. The child process object ceases to be
    /// valid after this call.
    ///
    /// \remark <b>Blocking remarks</b>: This call blocks if the child
    /// process has not finalized execution and waits until it terminates.
    ///
    const status wait();

#if defined(BOOST_POSIX_API)

    /// \brief Maps child's file descriptors to file handles.
    ///
    typedef std::map<int, pipe_end> pe_map;

    /// \brief Constructs a new POSIX child object representing a just
    ///        spawned child process.
    ///
    child(native_type native, const pe_map& pes)
    : pid_(native)
    , pes_(pes)
    {
        chin_ = pes_[STDIN_FILENO];
        chout_ = pes_[STDOUT_FILENO];
        cherr_ = pes_[STDERR_FILENO];
    }


    /// \brief Gets a reference to the child's stream \a desc.
    ///
    /// Returns an object that represents one of the multiple communication
    /// channels with the child process. The parent can use the returned stream
    /// to send/recieve data to the child.
    ///
    /// Giving this function the STDIN_FILENO constant (defined in
    /// cstdlib) is a synonym for the get_stdin() call.
    ///
    /// Giving this function the STDOUT_FILENO or STDERR_FILENO constants
    /// (both defined in cstdlib) are synonyms for the get_stdout() and
    /// get_stderr() calls.
    ///
    pipe_end get_pipe_end(int desc) const;

#endif // #if defined(BOOST_POSIX_API)

private:
    /// \brief The native representaion of a child.
#if defined(BOOST_WINDOWS_API)
    shared_ptr<void> handle_;
#elif defined(BOOST_POSIX_API)
    pid_t pid_;
#endif

    /// \brief The standard input stream attached to the child process.
    ///
    /// This object holds the communication channel with the
    /// child's process standard input.
    ///
    pipe_end chin_;

    /// \brief The standard output stream attached to the child process.
    ///
    /// This object holds the communication channel with the
    /// child's process standard output.
    ///
    pipe_end chout_;

    /// \brief The standard error stream attached to the child process.
    ///
    /// This object holds the communication channel with the
    /// child's process standard error.
    ///
    pipe_end cherr_;

#if defined(BOOST_POSIX_API)

    /// \brief Contains all relationships between child's file
    ///        descriptors and their corresponding file descriptors.
    pe_map pes_;

#endif // #if defined(BOOST_POSIX_API)
};

/// \brief Collection of child objects.
///
/// This convenience type represents a collection of child objects backed
/// by a vector.
///
typedef std::vector<child> children;

} // namespace processes
} // namespace boost

#ifdef BOOST_MSVC
#pragma warning(pop)
#endif

#ifdef BOOST_PROCESSES_HEADER_ONLY
#include "boost/processes/detail/impl/child.hpp"
#endif

#endif
