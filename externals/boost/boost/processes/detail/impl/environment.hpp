// Boost.Process
//
// Copyright (c) 2006 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

#ifndef BOOST_PROCESSES_ENVIRONMENT_HPP_115014
#define BOOST_PROCESSES_ENVIRONMENT_HPP_115014

#include "boost/processes/environment.hpp"

#if defined(BOOST_POSIX_API)
extern "C" {
    extern char** environ;
}
#elif defined(BOOST_WINDOWS_API)
#include "boost/processes/detail/error.hpp"
#include "boost/processes/detail/windows_api.hpp"
// workaround a bug in WinBase.h
#ifdef GetEnvironmentStrings
#undef GetEnvironmentStrings
#endif
#endif

#define BOOST_PROCESSES_JOIN(a, b) a ## b

namespace boost {
namespace processes {

BOOST_PROCESSES_INLINE_IF_HEADER_ONLY
environment
current_environment()
{
    environment env;

#if defined(BOOST_POSIX_API)
    char** ptr = ::environ;
    while (*ptr != 0) {
        std::string str = *ptr;
        std::string::size_type pos = str.find('=');
        env.insert
            (environment::value_type(str.substr(0, pos),
                                     str.substr(pos + 1, str.length())));
        ptr++;
    }
#elif defined(BOOST_WINDOWS_API)
    using namespace detail::windows_api;
    char* es = GetEnvironmentStrings();
    if (es == 0)
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("GetEnvironmentStrings() failed");

    try {
        char* escp = es;
        while (*escp != '\0') {
            std::string str = escp;
            std::string::size_type pos = str.find('=');
            env.insert
                (environment::value_type(str.substr(0, pos),
                                         str.substr(pos + 1, str.length())));
            escp += str.length() + 1;
        }
    } catch (...) {
        FreeEnvironmentStringsA(es);
        throw;
    }

    FreeEnvironmentStringsA(es);
#endif

    return env;
}

} // namespace processes
} // namespace boost

#endif
