// Boost.Process
//
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

#ifndef BOOST_PROCESSES_LAUNCH_HPP_155059
#define BOOST_PROCESSES_LAUNCH_HPP_155059

#include "boost/processes/config.hpp"
#if defined(BOOST_WINDOWS_API)
#include "boost/processes/detail/impl/windows_launch.hpp"
#elif defined(BOOST_POSIX_API)
#include "boost/processes/detail/impl/posix_launch.hpp"
#endif

#endif
