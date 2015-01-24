// Boost.Process
//
// Copyright (c) 2006 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

#ifndef BOOST_PROCESSES_POSIX_UTIL_HPP_114253
#define BOOST_PROCESSES_POSIX_UTIL_HPP_114253

#include "boost/processes/detail/posix_util.hpp"
#include "boost/processes/detail/error.hpp"
#include "boost/assert.hpp"
#include <fcntl.h>
#include <unistd.h>
#include <cstring>

namespace boost {
namespace processes {
namespace detail {

BOOST_PROCESSES_INLINE_IF_HEADER_ONLY void
posix_dup2(int fd1, int fd2)
{
    if (::dup2(fd1, fd2) == -1)
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("dup2(2) failed");
}

BOOST_PROCESSES_INLINE_IF_HEADER_ONLY void
posix_close_on_exec(int fd, bool c)
{
    int old_flags = ::fcntl(fd, F_GETFD, 0);
    if (old_flags < 0)
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("fcntl(2) failed");
    int new_flags = c ? old_flags|FD_CLOEXEC : old_flags&~FD_CLOEXEC;
    if (new_flags != old_flags && ::fcntl(fd, F_SETFD, new_flags) == -1)
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("fcntl(2) failed");
}

BOOST_PROCESSES_INLINE_IF_HEADER_ONLY std::pair<std::size_t, char**>
vector_to_argv(const std::vector<std::string>& args)
{
    typedef std::vector<std::string> arguments;
    arguments::size_type nargs = args.size();
    BOOST_ASSERT(nargs > 0);

    char** argv = new char*[nargs + 1];
    arguments::size_type i = 0;
    for (arguments::const_iterator iter = args.begin();
         iter != args.end(); iter++) {
        argv[i] = ::strdup(iter->c_str());
        i++;
    }
    argv[nargs] = 0;

    return std::pair<std::size_t, char **>(nargs, argv);
}

BOOST_PROCESSES_INLINE_IF_HEADER_ONLY char**
environment_to_envp(const environment& env)
{
    char** ep = new char*[env.size() + 1];
    environment::size_type i = 0;
    for (environment::const_iterator iter = env.begin();
         iter != env.end(); iter++) {
        std::string tmp = iter->first + "=" + iter->second;

        char* cstr = new char[tmp.length() + 1];
        std::strncpy(cstr, tmp.c_str(), tmp.length());
        cstr[tmp.length()] = '\0';

        ep[i++] = cstr;
    }
    ep[i] = 0;
    return ep;
}

} // namespace detail
} // namespace processes
} // namespace boost

#endif
