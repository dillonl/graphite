// Boost.Process
//
// Copyright (c) 2006, 2007 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

/// \file boost/process/self.hpp
///
/// Includes the declaration of the self class.

#ifndef BOOST_PROCESSES_SELF_HPP_172258
#define BOOST_PROCESSES_SELF_HPP_172258

#include "boost/processes/environment.hpp"
#include "boost/processes/process.hpp"
#include "boost/noncopyable.hpp"
#if defined(BOOST_POSIX_API)
    #include <sys/types.h>
    #include <unistd.h>
#elif defined(BOOST_WINDOWS_API)
    #include "boost/processes/detail/windows_api.hpp"
#endif

namespace boost {
namespace processes {

/// \brief Generic implementation of the Process concept.
///
/// The self singleton provides access to the current process.
///
class BOOST_PROCESSES_DECL self: boost::noncopyable
{
    /// \brief Constructs a new self object.
    ///
    /// Creates a new self object that represents the current process.
    ///
    self()
    : native_id_
#if defined(BOOST_POSIX_API)
        (::getpid())
#elif defined(BOOST_WINDOWS_API)
        (detail::windows_api::GetCurrentProcessId())
#endif
    {}

    process::native_id_type native_id_;

public:
    /// \brief Returns the self instance representing the caller's process.
    ///
    /// Returns a reference to the self instance that represents the
    /// caller's process.
    ///
    static self& get_instance()
    {
        static self *instance = 0;
        if (instance == 0)
            instance = new self;
        return *instance;
    }

    /// \brief Returns the process ID.
    ///
    process::native_id_type native_id() const
    {
        return native_id_;
    }

    /// \brief Returns the current environment.
    ///
    /// Returns the current process' environment variables. Modifying the
    /// returned object has no effect on the current environment.
    ///
    const environment get_environment() const
    {
        return current_environment();
    }

    // TODO: Add methods to modify the current environment.
    // A preliminary approach could be:
    // has_environment(variable)
    // get_environment(variable)
    // set_environment(variable, value)
    // unset_environment(variable)
};

} // namespace processes
} // namespace boost

#endif
