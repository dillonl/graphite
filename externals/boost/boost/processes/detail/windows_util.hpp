// Boost.Process
//
// Copyright (c) 2006 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

/// \file boost/process/detail/windows_util.hpp
///
/// Provides some convenience functions to start processes under Windows
/// operating systems.

#ifndef BOOST_PROCESSES_WINDOWS_UTIL_HPP_170710
#define BOOST_PROCESSES_WINDOWS_UTIL_HPP_170710

#include "boost/processes/config.hpp"
#ifndef BOOST_WINDOWS_API
#error "Unsupported platform."
#endif
#include "boost/processes/detail/handle.hpp"
#include "boost/processes/detail/windows_api.hpp"
#include "boost/processes/environment.hpp"
#include "boost/shared_array.hpp"
#include <string>
#include <vector>

namespace boost {
namespace processes {
namespace detail {

/// \brief Converts an environment to a string used by CreateProcess().
///
/// Converts the environment's contents to the format used by the
/// CreateProcess() system call. The returned char* string is
/// allocated in dynamic memory and the caller must free it when not
/// used any more. This is enforced by the use of a shared pointer.
/// The string is of the form var1=value1\\0var2=value2\\0\\0.
///
BOOST_PROCESSES_DECL boost::shared_array<char>
environment_to_windows_string(const environment& env);

/// \brief Duplicates the \a h native file handle.
///
/// Given a \a file handle \a h, this routine constructs a new
/// \a file_handle object that owns a new duplicate of \a h. The
/// duplicate's inheritable flag is set to the value of \a inheritable.
///
/// This operation is only available in Windows systems.
///
/// \return A file handle owning a duplicate of \a h.
/// \throw system_error If DuplicateHandle() fails.
///
BOOST_PROCESSES_DECL handle
windows_dup(windows_api::HANDLE h, bool inheritable);

/// \brief Changes the file handle's inheritable flag.
///
/// Changes the file handle's inheritable flag to \a i. It is not
/// necessary for the file handle's flag to be different than \a i.
///
/// This operation is only available in Windows systems.
///
/// \pre The file handle is valid. // FIXME?
/// \post The native file handle's inheritable flag is set to \a i.
/// \throw system_error If the property change fails.
///
BOOST_PROCESSES_DECL void
windows_set_inheritable(windows_api::HANDLE h, bool inheritable);

/// \brief Converts the command line to a plain string.
///
/// Converts the command line's list of arguments to the format
/// expected by the \a lpCommandLine parameter in the CreateProcess()
/// system call.
///
/// This operation is only available in Windows systems.
///
/// \return A dynamically allocated string holding the command line
///         to be passed to the executable. It is returned in a
///         shared_array object to ensure its release at some point.
BOOST_PROCESSES_DECL boost::shared_array<char>
vector_to_windows_cmdline(const std::vector<std::string>& args);

} // namespace detail
} // namespace processes
} // namespace boost

#ifdef BOOST_PROCESSES_HEADER_ONLY
#include "boost/processes/detail/impl/windows_util.hpp"
#endif

#endif
