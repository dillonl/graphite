// Boost.Process
//
// Copyright (c) 2006 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

#ifndef BOOST_PROCESSES_WINDOWS_UTIL_HPP_183852
#define BOOST_PROCESSES_WINDOWS_UTIL_HPP_183852

#include "boost/processes/detail/windows_util.hpp"
#include "boost/processes/detail/error.hpp"
#include "boost/processes/detail/windows_api.hpp"
#include "boost/assert.hpp"
#include <cstring>

namespace boost {
namespace processes {
namespace detail {

BOOST_PROCESSES_INLINE_IF_HEADER_ONLY handle
windows_dup(windows_api::HANDLE h, bool inheritable)
{
    using namespace windows_api;
    HANDLE result;
    if (DuplicateHandle(GetCurrentProcess(), h,
                        GetCurrentProcess(), &result,
                        0, inheritable, DUPLICATE_SAME_ACCESS_) == 0)
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("DuplicateHandle() failed");
    return handle(result);
}

BOOST_PROCESSES_INLINE_IF_HEADER_ONLY void
windows_set_inheritable(windows_api::HANDLE h, bool inheritable)
{
    using namespace windows_api;
    BOOST_ASSERT(h != INVALID_HANDLE_VALUE_); // FIXME?
    if (SetHandleInformation(h, HANDLE_FLAG_INHERIT_,
                             inheritable ? HANDLE_FLAG_INHERIT_ : 0) == 0)
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("SetHandleInformation() failed");
}

BOOST_PROCESSES_INLINE_IF_HEADER_ONLY boost::shared_array<char>
vector_to_windows_cmdline(const std::vector<std::string>& args)
{
    typedef std::vector<std::string> arguments;
    arguments args2;

    arguments::size_type i = 0;
    std::size_t length = 0;
    for (arguments::const_iterator iter = args.begin();
         iter != args.end(); iter++) {
        std::string arg = (*iter);

        std::string::size_type pos = 0;
        while ((pos = arg.find('"', pos)) != std::string::npos) {
            arg.replace(pos, 1, "\\\"");
            pos += 2;
        }

        if (arg.find(' ') != std::string::npos)
            arg = '\"' + arg + '\"';

        if (i != args.size() - 1)
            arg += ' ';

        args2.push_back(arg);
        length += arg.size() + 1;

        i++;
    }

    boost::shared_array<char> cmdline(new char[length]);
    cmdline[0] = 0;
    for (arguments::size_type i = 0; i < args2.size(); i++)
        std::strncat(cmdline.get(), args2[i].c_str(), args2[i].size());

    return cmdline;
}

/// \brief Converts an environment to a string used by CreateProcess().
///
/// Converts the environment's contents to the format used by the
/// CreateProcess() system call. The returned char* string is
/// allocated in dynamic memory and the caller must free it when not
/// used any more. This is enforced by the use of a shared pointer.
/// The string is of the form var1=value1\\0var2=value2\\0\\0.
///
BOOST_PROCESSES_INLINE_IF_HEADER_ONLY boost::shared_array<char>
environment_to_windows_string(const environment& env)
{
    boost::shared_array<char> strs(0);

    // TODO: Add the "" variable to the returned string; it shouldn't
    // be in the environment if the user didn't add it.

    if (env.size() == 0) {
        strs.reset(new char[2]);
        strs[0] = 0;
        strs[1] = 0;
    } else {
        std::string::size_type len = sizeof(char);
        for (environment::const_iterator iter = env.begin();
             iter != env.end(); iter++)
            len += (iter->first.length() + 1 + iter->second.length() +
                    1) * sizeof(char);

        strs.reset(new char[len]);

        char* ptr = strs.get();
        for (environment::const_iterator iter = env.begin();
             iter != env.end(); iter++) {
            std::string tmp = iter->first + "=" + iter->second;
            std::strcpy(ptr, tmp.c_str());
            ptr += (tmp.length() + 1) * sizeof(char);

            BOOST_ASSERT(static_cast<std::string::size_type>
                (ptr - strs.get()) * sizeof(char) < len);
        }
        *ptr = '\0';
    }

    BOOST_ASSERT(strs.get() != 0);
    return strs;
}

} // namespace detail
} // namespace processes
} // namespace boost

#endif
