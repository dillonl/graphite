// Boost.Process
//
// Copyright (c) 2006 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

/// \file boost/process/detail/posix_util.hpp
///
/// Provides some convenience functions to start processes under POSIX
/// operating systems.

#ifndef BOOST_PROCESSES_POSIX_UTIL_HPP_170432
#define BOOST_PROCESSES_POSIX_UTIL_HPP_170432

#include "boost/processes/config.hpp"
#ifndef BOOST_POSIX_API
#error "Unsupported platform."
#endif
#include "boost/processes/environment.hpp"
#include <string>
#include <vector>
#include <utility>

namespace boost {
namespace processes {
namespace detail {

/// \brief Duplicates an open native file handle.
///
/// Given a native file handle \a h1, this routine duplicates it so
/// that it ends up being identified by the native file handle \a h2.
///
/// This operation is only available in POSIX systems.
///
/// \throw system_error If dup2() fails.
///
BOOST_PROCESSES_DECL void
posix_dup2(int fd1, int fd2);

BOOST_PROCESSES_DECL void
posix_close_on_exec(int fd, bool close_on_exec);

/// \brief Converts the command line to an array of C strings.
///
/// Converts the command line's list of arguments to the format expected
/// by the \a argv parameter in the POSIX execve() system call.
///
/// This operation is only available in POSIX systems.
///
/// \return The first argument of the pair is an integer that indicates
///         how many strings are stored in the second argument. The
///         second argument is a null-terminated, dynamically allocated
///         vector of dynamically allocated strings holding the arguments
///         to the executable. The caller is responsible of freeing them.
BOOST_PROCESSES_DECL std::pair<std::size_t, char**>
vector_to_argv(const std::vector<std::string>& args);

/// \brief Converts an environment to a char** table as used by execve().
///
/// Converts the environment's contents to the format used by the
/// execve() system call. The returned char** array is allocated
/// in dynamic memory and the caller must free it when not used any
/// more. Each entry is also allocated in dynamic memory and is a
/// null-terminated string of the form var=value; these must also be
/// released by the caller.
///
BOOST_PROCESSES_DECL char**
environment_to_envp(const environment& env);

} // namespace detail
} // namespace processes
} // namespace boost

#ifdef BOOST_PROCESSES_HEADER_ONLY
#include "boost/processes/detail/impl/posix_util.hpp"
#endif

#endif
