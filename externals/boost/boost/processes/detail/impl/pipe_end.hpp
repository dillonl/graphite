// Boost.Process
//
// Inspired by fdstream.hpp, (C) Copyright Nicolai M. Josuttis 2001,
// available at http://www.josuttis.com/cppcode/fdstream.html.
// (C) Copyright 2003-2007 Jonathan Turkanis
// (C) Copyright 2008 CodeRage, LLC (turkanis at coderage dot com)
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

#ifndef BOOST_PROCESSES_PIPE_END_HPP_144148
#define BOOST_PROCESSES_PIPE_END_HPP_144148

#include "boost/processes/pipe_end.hpp"

namespace boost {
namespace processes {

std::streamsize BOOST_PROCESSES_INLINE_IF_HEADER_ONLY
pipe_end::read(char* s, std::streamsize n)
{
#if defined(BOOST_WINDOWS_API)
    using namespace detail::windows_api;
    DWORD result;
    if (!ReadFile(native(), s, n, &result, 0))
    {
        DWORD le = GetLastError();
        if (ERROR_BROKEN_PIPE_ == le)
            result = -1; // eof
        else
            boost::throw_exception(boost::system::system_error(
                boost::system::error_code(
                    le, boost::system::get_system_category()),
                BOOST_PROCESSES_SOURCE_LOCATION "ReadFile() failed"));
    }
    return static_cast<std::streamsize>(result);
#elif defined(BOOST_POSIX_API)
    errno = 0;
    std::streamsize result = ::read(native(), s, n);
    if (errno != 0)
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("read(2) failed");
    return result == 0 ? -1 : result;
#endif
}

std::streamsize BOOST_PROCESSES_INLINE_IF_HEADER_ONLY
pipe_end::write(const char* s, std::streamsize n)
{
#if defined(BOOST_WINDOWS_API)
    detail::windows_api::DWORD ignore;
    if (!detail::windows_api::WriteFile(native(), s, n, &ignore, 0))
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("WriteFile() failed");
    return n;
#elif defined(BOOST_POSIX_API)
    int amt = ::write(native(), s, n);
    if (amt < n) // Handles blocking fd's only.
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("write(2) failed");
    return n;
#endif
}

} // namespace processes {
} // namespace boost {

#endif
