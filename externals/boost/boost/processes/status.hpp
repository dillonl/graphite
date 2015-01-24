// Boost.Process
//
// Copyright (c) 2006 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

/// \file boost/process/status.hpp
///
/// Includes the declaration of the status class.

#ifndef BOOST_PROCESSES_STATUS_HPP_172235
#define BOOST_PROCESSES_STATUS_HPP_172235

#include "boost/processes/config.hpp"
#include "boost/assert.hpp"
#if defined(BOOST_POSIX_API)
    #include <sys/wait.h>
#endif

namespace boost {
namespace processes {

/// \brief Status returned by a finalized child process.
///
/// This class represents the %status returned by a child process after it
/// has terminated. It only provides that information available under all
/// supported platforms.
///
/// \see posix::status
///
class BOOST_PROCESSES_DECL status
{
public:
    /// \brief Creates a status object based on exit information.
    ///
    /// Creates a new status object representing the exit status of a
    /// child process.
    ///
    /// \param flags In a POSIX system this parameter contains the
    ///              flags returned by the ::waitpid() call. In a
    ///              Windows system it contains the exit code only.
    status(int flags);

    /// \brief Returns whether the process exited gracefully or not.
    ///
    /// Returns whether the process exited gracefully or not.
    ///
    bool exited() const;

    /// \brief If exited, returns the exit code.
    ///
    /// If the process exited, returns the exit code it returned.
    ///
    /// \pre exited() is true.
    ///
    int exit_status() const;

#ifdef BOOST_POSIX_API

    /// \brief Returns whether the process exited due to an external
    ///        signal.
    ///
    /// Returns whether the process exited due to an external signal.
    /// The result is always false in Windows systems.
    ///
    bool signaled() const;

    /// \brief If signaled, returns the terminating signal code.
    ///
    /// If the process was signaled, returns the terminating signal code.
    /// Cannnot be called under Windows because the preconditions will not
    /// ever be met.
    ///
    /// \pre signaled() is true.
    ///
    int term_signal() const;

    /// \brief If signaled, returns whether the process dumped core.
    ///
    /// If the process was signaled, returns whether the process
    /// produced a core dump.
    /// Cannnot be called under Windows because the preconditions will not
    /// ever be met.
    ///
    /// \pre signaled() is true.
    ///
    bool dumped_core() const;

    /// \brief Returns whether the process was stopped by an external
    ///        signal.
    ///
    /// Returns whether the process was stopped by an external signal.
    /// The result is always false in Windows systems.
    ///
    bool stopped() const;

    /// \brief If stpped, returns the stop signal code.
    ///
    /// If the process was stopped, returns the stop signal code.
    /// Cannnot be called under Windows because the preconditions will not
    /// ever be met.
    ///
    /// \pre signaled() is true.
    ///
    int stop_signal() const;

#endif // #ifdef BOOST_POSIX_API

protected:
    /// \brief OS-specific codification of exit status.
    ///
    int flags_;
};

inline
status::status(int flags)
: flags_(flags)
{}

inline
bool
status::exited() const
{
#if defined(BOOST_POSIX_API)
    return WIFEXITED(flags_);
#elif defined(BOOST_WINDOWS_API)
    return true;
#endif
}

inline
int
status::exit_status() const
{
    BOOST_ASSERT(exited());
#if defined(BOOST_POSIX_API)
    return WEXITSTATUS(flags_);
#elif defined(BOOST_WINDOWS_API)
    return flags_;
#endif
}

#ifdef BOOST_POSIX_API

inline bool status::signaled() const
{
    return WIFSIGNALED(flags_);
}

inline int status::term_signal() const
{
    BOOST_ASSERT(signaled());
    return WTERMSIG(flags_);
}

inline bool status::dumped_core() const
{
    BOOST_ASSERT(signaled());
    return WCOREDUMP(flags_);
}

inline bool status::stopped() const
{
    return WIFSTOPPED(flags_);
}

inline int status::stop_signal() const
{
    BOOST_ASSERT(stopped());
    return WSTOPSIG(flags_);
}

#endif // #ifdef BOOST_POSIX_API

} // namespace processes
} // namespace boost

#endif
