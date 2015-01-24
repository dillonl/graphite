// Boost.Process
//
// Copyright (c) 2006, 2007 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

/// \file boost/process/process.hpp
///
/// Includes the declaration of the process class.

#ifndef BOOST_PROCESSES_PROCESS_HPP_172335
#define BOOST_PROCESSES_PROCESS_HPP_172335

#include "boost/processes/config.hpp"
#if defined(BOOST_POSIX_API)
    #include <sys/types.h>
#elif defined(BOOST_WINDOWS_API)
    #include "boost/processes/detail/windows_api.hpp"
#endif

namespace boost {
namespace processes {

// reserve this name for future uses
class process { process(); public:
#if defined(BOOST_PROCESS_DOXYGEN)
    /// \brief The native process ID type.
    typedef implementation_defined native_id_type;
#elif defined(BOOST_WINDOWS_API)
    typedef detail::windows_api::DWORD native_id_type;
#elif defined(BOOST_POSIX_API)
    typedef pid_t native_id_type;
#endif
};

} // namespace processes
} // namespace boost

#endif
