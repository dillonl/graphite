// Boost.Process
//
// Copyright (c) 2006 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

/// \file boost/process/context.hpp
///
/// Includes the declaration of the context class and several accessory
/// base classes.

#ifndef BOOST_PROCESSES_CONTEXT_HPP_171940
#define BOOST_PROCESSES_CONTEXT_HPP_171940

#include "boost/processes/environment.hpp"
#include "boost/processes/detail/error.hpp"
#include "boost/processes/stream_behavior.hpp"
#include "boost/assert.hpp"
#if defined(BOOST_WINDOWS_API)
    #include "boost/processes/detail/windows_api.hpp"
#endif
#include <deque>
#include <string>
#include <vector>

namespace boost {
namespace processes {

typedef std::deque<stream_behavior> streams_behavior;

/// \brief Process startup execution context.
///
/// The context structure groups all the parameters needed to configure a
/// process' environment during its creation.
///
template<class Path>
struct basic_context
{
    /// \brief Constructs a new context.
    ///
    /// Constructs a new context. It is configured as follows:
    /// * The initial work directory of the child processes is set to the
    ///   current working directory.
    /// * The environment variables table is empty.
    /// * All communcation channels with the child process are closed.
    /// * There are no channel mergings.
    ///
    basic_context();

    /// \brief The process' initial work directory.
    ///
    /// The work directory is the directory in which the process starts
    /// execution.
    ///
    Path work_directory;

    /// \brief The process' environment.
    ///
    /// Contains the list of environment variables, alongside with their
    /// values, that will be passed to the spawned child process.
    ///
    boost::processes::environment environment;

    boost::processes::streams_behavior streams_behavior;

    basic_context& add(const stream_behavior& sb)
    {
        streams_behavior.push_back(sb);
        return *this;
    }

#if defined(BOOST_WINDOWS_API)

    /// \brief Process startup information.
    ///
    void* startupinfo;

#elif defined(BOOST_POSIX_API)

    /// \brief The effective user credentials.
    ///
    /// EUID that specifies the effective user credentials to use to run
    /// the %child process. Defaults to the current EUID.
    ///
    uid_t euid;

    /// \brief The effective group credentials.
    ///
    /// EGID that specifies the effective group credentials to use to run
    /// the %child process. Defaults to the current EGID.
    ///
    uid_t egid;

    /// \brief The chroot directory, if any.
    ///
    /// Specifies the directory in which the %child process is chrooted
    /// before execution. Empty if this feature is not desired.
    ///
    Path chroot;

    /// \brief Apply properties to the current execution environment.
    ///
    /// Modifies the current execution environment so that the properties become
    /// effective.
    ///
    /// \throw system_error If any error ocurred during environment
    ///                     configuration.
    void apply() const;

#endif // #elif defined(BOOST_POSIX_API)
};

typedef basic_context<std::string> context;

/// \brief Represents a child process in a pipeline.
///
/// This convenience class is a triplet that holds all the data required
/// to spawn a new child process in a pipeline.
///
template<class Executable, class Arguments, class Context>
class basic_pipeline_entry
{
public:
    /// \brief The executable to launch.
    ///
    Executable executable;

    /// \brief The set of arguments to pass to the executable.
    ///
    Arguments arguments;

    /// \brief The child's execution context.
    ///
    Context context;

    /// \brief The type of the Executable concept used in this template
    ///        instantiation.
    typedef Executable executable_type;

    /// \brief The type of the Arguments concept used in this template
    ///        instantiation.
    typedef Arguments arguments_type;

    /// \brief The type of the Context concept used in this template
    ///        instantiation.
    typedef Context context_type;

    /// \brief Constructs a new pipeline_entry object.
    ///
    /// Given the executable, set of arguments and execution triplet,
    /// constructs a new pipeline_entry object that holds the three
    /// values.
    ///
    basic_pipeline_entry(const Executable& e, const Arguments& a,
                         const Context& c) :
        executable(e),
        arguments(a),
        context(c)
    {}
};

/// \brief Default instantiation of basic_pipeline_entry.
///
typedef basic_pipeline_entry<std::string, std::vector<std::string>,
                             context> pipeline_entry;

template<class Path>
basic_context<Path>::basic_context()
// FIXME: temporary workaround
: environment(current_environment())
#if defined(BOOST_WINDOWS_API)
,  startupinfo(0)
#elif defined(BOOST_POSIX_API)
, euid(::geteuid())
, egid(::getegid())
#endif // #elif defined(BOOST_POSIX_API)
{
#if defined(BOOST_POSIX_API)
    const char* buf = ::getcwd(0, 0);
    if (buf == 0)
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("getcwd(2) failed");
    work_directory = buf;
#elif defined(BOOST_WINDOWS_API)
    using namespace detail::windows_api;
    DWORD length = GetCurrentDirectoryA(0, 0);
    char* buf = new char[length * sizeof(char)];
    if (GetCurrentDirectoryA(length, buf) == 0) {
        delete buf;
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("GetCurrentDirectory() failed");
    }
    work_directory = buf;
    delete buf;
#endif
    BOOST_ASSERT(!work_directory.empty());
}

#if defined(BOOST_POSIX_API)
template<class Path>
void
basic_context<Path>::apply() const
{
    if (!chroot.empty())
        if (::chroot(chroot.c_str()) == -1)
            BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("chroot(2) failed");
    if (egid != ::getegid())
        if (::setegid(egid) == -1)
            BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("setegid(2) failed");
    if (euid != ::geteuid())
        if (::seteuid(euid) == -1)
            BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("seteuid(2) failed");
    BOOST_ASSERT(!work_directory.empty());
    if (::chdir(work_directory.c_str()) == -1)
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("chdir(2) failed");
}
#endif // #if defined(BOOST_POSIX_API)

} // namespace processes
} // namespace boost

#endif
