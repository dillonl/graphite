// Boost.Process
//
// Copyright (c) 2006 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

#ifndef BOOST_PROCESSES_LAUNCH_SHELL_HPP_125031
#define BOOST_PROCESSES_LAUNCH_SHELL_HPP_125031

#include "boost/processes/launch_shell.hpp"
#include "boost/processes/launch.hpp"

namespace boost {
namespace processes {

BOOST_PROCESSES_INLINE_IF_HEADER_ONLY child
launch_shell(const std::string& command, const context& ctx)
{
    std::string exe;
    std::vector<std::string> args;

#if defined(BOOST_POSIX_API)
    exe = "/bin/sh";
    args.push_back("sh");
    args.push_back("-c");
    args.push_back(command);
#elif defined(BOOST_WINDOWS_API)
    using namespace detail::windows_api;
    char buf[MAX_PATH_];
    unsigned int res =
        GetSystemDirectoryA(buf, MAX_PATH_);
    if (res == 0)
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("GetSystemDirectory() failed");
    BOOST_ASSERT(res < MAX_PATH_);

    exe = std::string(buf) + "\\cmd.exe";
    args.push_back("cmd");
    args.push_back("/c");
    args.push_back(command);
#endif

    return launch(exe, args, ctx);
}

} // namespace processes
} // namespace boost

#endif
