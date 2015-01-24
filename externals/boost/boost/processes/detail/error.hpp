// Boost.Process
//
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

#ifndef BOOST_PROCESSES_ERROR_HPP_221915
#define BOOST_PROCESSES_ERROR_HPP_221915

#include "boost/processes/config.hpp"

#include "boost/preprocessor/stringize.hpp"
#include "boost/system/system_error.hpp"
#include "boost/throw_exception.hpp"

#if defined(BOOST_WINDOWS_API)
#include "boost/processes/detail/windows_api.hpp"
#elif defined(BOOST_POSIX_API)
#include <cerrno>
#endif

#define BOOST_PROCESSES_SOURCE_LOCATION \
    "in file '" __FILE__ "', line " BOOST_PP_STRINGIZE(__LINE__) ": "

#if defined(BOOST_WINDOWS_API)
    #define BOOST_PROCESSES_LAST_ERROR detail::windows_api::GetLastError()
#elif defined(BOOST_POSIX_API)
    #define BOOST_PROCESSES_LAST_ERROR errno
#endif

#define BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR(what) \
    boost::throw_exception(boost::system::system_error( \
        boost::system::error_code(BOOST_PROCESSES_LAST_ERROR, \
                                  boost::system::get_system_category()), \
        BOOST_PROCESSES_SOURCE_LOCATION what))

#endif
