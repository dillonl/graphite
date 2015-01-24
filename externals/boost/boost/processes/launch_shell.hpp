// Boost.Process
//
// Copyright (c) 2006 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

/// \file boost/process/launch_shell.hpp
///
/// Provides the launch_shell() function.

#ifndef BOOST_PROCESSES_LAUNCH_SHELL_HPP_124744
#define BOOST_PROCESSES_LAUNCH_SHELL_HPP_124744

#include "boost/processes/config.hpp"
#include "boost/processes/child.hpp"
#include "boost/processes/context.hpp"
#include "boost/processes/detail/error.hpp"
#include <string>

namespace boost {
namespace processes {

/// \brief Launches a shell-based command.
///
/// Executes the given command through the default system shell. The
/// command is subject to pattern expansion, redirection and pipelining.
/// The shell is launched as described by the parameters in the execution
/// context.
///
/// This function behaves similarly to the system(3) system call. In a
/// POSIX system, the command is fed to /bin/sh whereas under a Windows
/// system, it is fed to cmd.exe. It is difficult to write portable
/// commands as the first parameter, but this function comes in handy in
/// multiple situations.
///
BOOST_PROCESSES_DECL child
launch_shell(const std::string& command, const context& ctx);

} // namespace processes
} // namespace boost

#ifdef BOOST_PROCESSES_HEADER_ONLY
#include "boost/processes/detail/impl/launch_shell.hpp"
#endif

#endif
