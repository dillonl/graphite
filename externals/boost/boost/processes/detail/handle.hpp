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

#ifndef BOOST_PROCESSES_HANDLE_HPP_121403
#define BOOST_PROCESSES_HANDLE_HPP_121403

#include "boost/processes/config.hpp"
#include "boost/processes/detail/error.hpp"
#include "boost/noncopyable.hpp"
#include "boost/shared_ptr.hpp"

#if defined(BOOST_POSIX_API)
    #include <unistd.h>
#elif defined(BOOST_WINDOWS_API)
    #include "boost/processes/detail/windows_api.hpp"
#endif

#ifdef BOOST_MSVC
#pragma warning(push)
// 4251 * needs to have dll-interface to be used by clients of *
// 4275 * used as base for dll-interface class *
#pragma warning(disable: 4251 4275)
#endif

namespace boost {
namespace processes {
namespace detail {

class BOOST_PROCESSES_DECL handle {
public:
#if defined(BOOST_WINDOWS_API)
    typedef detail::windows_api::HANDLE native_type;
#elif defined(BOOST_POSIX_API)
    typedef int native_type;
#endif

    static native_type invalid_native()
    {
#if defined(BOOST_WINDOWS_API)
        return detail::windows_api::INVALID_HANDLE_VALUE_;
#elif defined(BOOST_POSIX_API)
        return -1;
#endif
    }

    handle() {} // never throws

    explicit handle(native_type native)
    : impl_(new impl(native))
    {}

    bool valid() const
    { return impl_ && impl_->valid(); }

    void close()
    { if (impl_) impl_->close(); }

    native_type native() const
    { return impl_ ? impl_->native() : invalid_native(); }

private:
    class BOOST_PROCESSES_DECL impl: boost::noncopyable {
    public:
        typedef handle::native_type native_type;

        impl(native_type native): native_(native) {}

        ~impl()
        {
            if (valid())
            {
#if defined(BOOST_WINDOWS_API)
                detail::windows_api::CloseHandle(native_);
#elif defined(BOOST_POSIX_API)
                ::close(native_);
#endif
            }
        }

        bool valid() const
        { return native_ != handle::invalid_native(); }

        // FIXME: throw on error?
        void close()
        {
            if (valid())
            {
#if defined(BOOST_WINDOWS_API)
                detail::windows_api::CloseHandle(native_);
#elif defined(BOOST_POSIX_API)
                ::close(native_);
#endif
                native_ = invalid_native();
            }
        }

        native_type native()
        { return native_; }

    private:
        native_type native_;
    }; // class impl

    // FIXME: use intrusive_ptr?
    shared_ptr<impl> impl_;
}; // class handle

} // namespace detail {
} // namespace processes {
} // namespace boost {

#ifdef BOOST_MSVC
#pragma warning(pop)
#endif

#endif
