// Boost.Process
//
// Copyright (c) 2006 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

/// \file boost/process/launch.hpp
///
/// Provides launch() function.

#ifndef BOOST_PROCESSES_LAUNCH_HPP_172625
#define BOOST_PROCESSES_LAUNCH_HPP_172625

#include "boost/processes/child.hpp"
#include "boost/processes/context.hpp"
#include <string>
#include <vector>

namespace boost {
namespace processes {

/// \brief Starts a new child process.
///
/// Launches a new process based on the binary image specified by the
/// executable, the set of arguments to be passed to it and several
/// parameters that describe the execution context.
///
/// \remark <b>Blocking remarks</b>: This function may block if the device
/// holding the executable blocks when loading the image. This might
/// happen if, e.g., the binary is being loaded from a network share.
///
/// \return A handle to the new child process.
///
BOOST_PROCESSES_DECL child
launch(const std::string& exe,
       const std::vector<std::string>& args,
       const context& ctx);

} // namespace processes
} // namespace boost

#ifdef BOOST_PROCESSES_HEADER_ONLY
#include "boost/processes/detail/impl/launch.hpp"
#endif

#endif
