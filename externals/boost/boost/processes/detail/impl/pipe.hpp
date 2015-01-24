// Boost.Process
//
// Copyright (c) 2006 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

#ifndef BOOST_PROCESSES_PIPE_HPP_114312
#define BOOST_PROCESSES_PIPE_HPP_114312

#include "boost/processes/detail/pipe.hpp"
#include "boost/processes/detail/error.hpp"
#if defined(BOOST_POSIX_API)
    #include <unistd.h>
#elif defined(BOOST_WINDOWS_API)
    #include "boost/processes/detail/windows_api.hpp"
#endif

namespace boost {
namespace processes {
namespace detail {

BOOST_PROCESSES_INLINE_IF_HEADER_ONLY
pipe::pipe()
{
    pipe_end::native_type hs[2];

 #if defined(BOOST_WINDOWS_API)
    if (!detail::windows_api::CreatePipe(&hs[0], &hs[1], 0, 0))
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("CreatePipe() failed");
#else
    if (::pipe(hs) == -1)
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("pipe(2) failed");
#endif

    read_end_ = pipe_end(hs[0]);
    write_end_ = pipe_end(hs[1]);
}

} // namespace detail
} // namespace processes
} // namespace boost

#endif
