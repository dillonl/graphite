// Boost.Process
//
// Copyright (c) 2006 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

#ifndef BOOST_PROCESSES_OPERATIONS_HPP_114946
#define BOOST_PROCESSES_OPERATIONS_HPP_114946

#include "boost/processes/utility.hpp"
#include "boost/processes/detail/error.hpp"
#include "boost/assert.hpp"
#if defined(BOOST_POSIX_API)
    #include <unistd.h>
#elif defined(BOOST_WINDOWS_API)
    #include "boost/processes/detail/windows_api.hpp"
#endif

namespace boost {
namespace processes {

BOOST_PROCESSES_INLINE_IF_HEADER_ONLY std::string
find_executable_in_path(const std::string& file, std::string path)
{
#if defined(BOOST_POSIX_API)
    BOOST_ASSERT(file.find('/') == std::string::npos);
#elif defined(BOOST_WINDOWS_API)
    BOOST_ASSERT(file.find('\\') == std::string::npos);
#endif

    std::string result;

#if defined(BOOST_POSIX_API)
    if (path.empty()) {
        const char* envpath = ::getenv("PATH");
        if (envpath == 0)
            BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("getenv() failed");
        path = envpath;
    }
    BOOST_ASSERT(!path.empty());

    std::string::size_type pos1 = 0, pos2;
    do {
        pos2 = path.find(':', pos1);
        std::string dir = path.substr(pos1, pos2 - pos1);
        std::string f = dir + '/' + file;
        if (::access(f.c_str(), X_OK) == 0)
            result = f;
        pos1 = pos2 + 1;
    } while (pos2 != std::string::npos && result.empty());

    if (result.empty())
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("cannot locate file in path");

#elif defined(BOOST_WINDOWS_API)
    const char* exts[] = { "", ".exe", ".com", ".bat", 0 };
    const char** ext = exts;
    using namespace detail::windows_api;
    while (*ext != 0) {
        char buf[MAX_PATH_];
        char* dummy;
        DWORD len = SearchPathA(path.empty() ? 0 : path.c_str(),
                                file.c_str(), *ext, MAX_PATH_, buf, &dummy);
        BOOST_ASSERT(len < MAX_PATH_);
        if (len > 0) {
            result = buf;
            break;
        }
        ext++;
    }
    if (result.empty())
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("SearchPath() failed");
#endif

    return result;
}

} // namespace processes
} // namespace boost

#endif
