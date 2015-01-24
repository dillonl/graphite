// Boost.Process
//
// Copyright (c) 2006 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

/// \file boost/process/environment.hpp
///
/// Includes the declaration of the environment class and related
/// free functions.

#ifndef BOOST_PROCESSES_ENVIRONMENT_HPP_172757
#define BOOST_PROCESSES_ENVIRONMENT_HPP_172757

#include "boost/processes/config.hpp"
#include <map>
#include <string>

namespace boost {
namespace processes {

/// \brief Representation of a process' environment variables.
///
/// The environment is a map that stablishes an unidirectional
/// association between variable names and their values and is
/// represented by a string to string map.
///
/// Variables may be defined to the empty string. Be aware that doing so
/// is not portable: POSIX systems will treat such variables as being
/// defined to the empty value, but Windows systems are not able to
/// distinguish them from undefined variables.
///
/// Similarly, Windows systems support a variable with no name that holds
/// the path to the current working directory; you may set it if you want
/// to, but the library will do the job for you if unset. Contrarywise
/// POSIX systems do not support variables without names.
///
/// It is worthy to note that the environment is sorted alphabetically.
/// This is provided for-free by the map container used to implement this
/// type, and this behavior is required by Windows systems.
///
typedef std::map<std::string, std::string> environment;

/// \brief Returns a snapshot of the current environment.
///
/// This function grabs a snapshot of the current environment and returns
/// it to the caller. The returned object can be modified later on but
/// changes to it do \b not modify the current environment.
///
/// XXX self.environment() does the same as this. One of the two has to
/// go away, most likely this one.
///
BOOST_PROCESSES_DECL environment current_environment();

} // namespace processes
} // namespace boost

#ifdef BOOST_PROCESSES_HEADER_ONLY
#include "boost/processes/detail/impl/environment.hpp"
#endif

#endif
