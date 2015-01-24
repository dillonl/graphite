// Boost.Process
//
// Copyright (c) 2006 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

/// \file boost/process/stream_behavior.hpp
///
/// Includes the declaration of the stream_behavior class and associated
/// free functions.

#ifndef BOOST_PROCESSES_STREAM_BEHAVIOR_HPP_172214
#define BOOST_PROCESSES_STREAM_BEHAVIOR_HPP_172214

#include "boost/processes/config.hpp"
#include "boost/processes/detail/handle.hpp"
#include "boost/variant/get.hpp"
#ifdef BOOST_MSVC
#pragma warning(push)
// 4345 behavior change: an object of POD type constructed with a
//      an initializer of the form () will be default-initialized
#pragma warning(disable: 4345)
#endif
#include "boost/variant/variant.hpp"
#ifdef BOOST_MSVC
#pragma warning(pop)
#endif
#include <string>

namespace boost {
namespace processes {

enum std_fileno {
#ifdef BOOST_WINDOWS_API
    stdin_fileno,
    stdout_fileno,
    stderr_fileno
#elif defined(BOOST_POSIX_API)
    stdin_fileno = STDIN_FILENO,
    stdout_fileno = STDOUT_FILENO,
    stderr_fileno = STDERR_FILENO
#endif
};

#ifdef BOOST_PROCESS_DOXYGEN
typedef implementation_defined fileno;
#elif defined(BOOST_WINDOWS_API)
typedef std_fileno fileno;
#elif defined(BOOST_POSIX_API)
typedef int fileno;
#endif

class BOOST_PROCESSES_DECL fileno_holder
{
public:
    explicit fileno_holder(fileno fn): fileno_(fn) {}
    fileno fn() const { return fileno_; }

private:
    fileno fileno_;
};

/// The child's stream is connected to the same stream used by the parent.
///
class BOOST_PROCESSES_DECL inherit_stream: public fileno_holder
{
public:
    explicit inherit_stream(fileno f): fileno_holder(f) {}
};

/// The child's stream is connected to the parent by using an
/// anonymous pipe so that they can send and receive data to/from
/// each other.
///
class BOOST_PROCESSES_DECL capture_stream: public fileno_holder
{
public:
    enum mode_type { in, out };

    explicit capture_stream(std_fileno fn)
    : fileno_holder(fn)
    , mode_(default_mode(fn))
    {}

    capture_stream(fileno fn, mode_type mode)
    : fileno_holder(fn)
    , mode_(mode)
    { BOOST_ASSERT(mode == in || mode == out); }

    mode_type mode() const { return mode_; }

    static mode_type default_mode(std_fileno fn)
    {
        switch (fn) {
        case stdin_fileno:
            return in;
        case stdout_fileno:
        case stderr_fileno:
            return out;
        default:
            BOOST_ASSERT(false);
        }
        return mode_type(0); // supress warning
    }

private:
    mode_type mode_;
};

/// The child's stream is connected to child's standard out.
/// This is typically used when configuring the standard error
/// stream.
///
class BOOST_PROCESSES_DECL redirect_to_stdout: public fileno_holder
{
public:
    explicit redirect_to_stdout(fileno f): fileno_holder(f) {}
};

/// The child's stream is redirected to a null device so that its
/// in is always zero or its out is lost, depending on
/// whether the stream is an in or an out one. It is
/// important to notice that this is different from close because
/// the child is still able to write data. If we closed, e.g.
/// stdout, the child might not work at all!
///
class BOOST_PROCESSES_DECL silence_stream: public fileno_holder
{
public:
    enum mode_type { in, out };

    explicit silence_stream(std_fileno fn)
    : fileno_holder(fn)
    , mode_(default_mode(fn))
    {}

    silence_stream(fileno fn, mode_type mode)
    : fileno_holder(fn)
    , mode_(mode)
    { BOOST_ASSERT(mode == in || mode == out); }

    mode_type mode() const { return mode_; }

    static mode_type default_mode(std_fileno fn)
    {
        switch (fn) {
        case stdin_fileno:
            return in;
        case stdout_fileno:
        case stderr_fileno:
            return out;
        default:
            BOOST_ASSERT(false);
        }
        return mode_type(0); // supress warning
    }

private:
    mode_type mode_;
};

/// The child's stream is closed upon startup so that it will not
/// have any access to it.
///
class BOOST_PROCESSES_DECL close_stream: public fileno_holder
{
public:
    explicit close_stream(fileno f): fileno_holder(f) {}
};

class BOOST_PROCESSES_DECL redirect_to_named_file: public fileno_holder
{
public:
    enum mode_type { in, out };

    redirect_to_named_file(std_fileno fn, const std::string& filename)
    : fileno_holder(fn)
    , mode_(default_mode(fn))
    , filename_(filename)
    {}

    redirect_to_named_file(fileno fn,
                           mode_type mode,
                           const std::string& filename)
    : fileno_holder(fn)
    , mode_(mode)
    , filename_(filename)
    { BOOST_ASSERT(mode == in || mode == out); }

    mode_type mode() const { return mode_; }
    const std::string& filename() const { return filename_; }

    static mode_type default_mode(std_fileno fn)
    {
        switch (fn) {
        case stdin_fileno:
            return in;
        case stdout_fileno:
        case stderr_fileno:
            return out;
        default:
            BOOST_ASSERT(false);
        }
        return mode_type(0); // supress warning
    }

private:
    mode_type mode_;
    std::string filename_;
};

/// The child redirects the stream's out to the provided native file.
///
class BOOST_PROCESSES_DECL redirect_to_native_file: public fileno_holder
{
public:
    typedef detail::handle::native_type native_type;

    redirect_to_native_file(fileno fn, native_type native_file)
    : fileno_holder(fn)
    , native_file_(native_file)
    {}

    native_type native_file() const { return native_file_; }

private:
    native_type native_file_;
};

/// \brief Describes the possible states for a child's communication stream.
///
typedef boost::variant<
    inherit_stream
  , capture_stream
  , redirect_to_stdout
  , silence_stream
  , close_stream
  , redirect_to_named_file
  , redirect_to_native_file
  > stream_behavior;

inline fileno get_fileno(const stream_behavior& sb)
{
    boost::variant<fileno_holder> t(sb);
    return boost::get<fileno_holder>(t).fn();
}

} // namespace processes
} // namespace boost

#endif
