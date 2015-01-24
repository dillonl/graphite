// Boost.Process
//
// Copyright (c) 2006 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

/// \file boost/process/utility.hpp
///
/// Provides miscellaneous utilities.

#ifndef BOOST_PROCESSES_UTILITY_HPP_205401
#define BOOST_PROCESSES_UTILITY_HPP_205401

#include "boost/processes/config.hpp"
#include <string>
#include <vector>

namespace boost {
namespace processes {

/// \brief Locates a program in the path.
///
/// Locates the executable program \a file in all the directory components
/// specified in \a path. If \a path is empty, the value of the PATH
/// environment variable is used.
///
/// The path variable is interpreted following the same conventions used
/// to parse the PATH environment variable in the underlying platform.
///
/// \throw system_error If the file cannot be found in the path.
///
BOOST_PROCESSES_DECL std::string
find_executable_in_path(const std::string& file, std::string path = "");

} // namespace processes
} // namespace boost

#ifdef BOOST_PROCESSES_HEADER_ONLY
#include "boost/processes/detail/impl/utility.hpp"
#endif

#endif
