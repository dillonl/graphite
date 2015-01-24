// Boost.Process
//
// Copyright (c) 2006 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

/// \file boost/process/detail/pipe.hpp
///
/// Includes the declaration of the pipe class. This file is for
/// internal usage only and must not be included by the library user.

#ifndef BOOST_PROCESSES_PIPE_HPP_170038
#define BOOST_PROCESSES_PIPE_HPP_170038

#include "boost/processes/config.hpp"
#include "boost/processes/pipe_end.hpp"

namespace boost {
namespace processes {
namespace detail {

/// \brief Simple RAII model for anonymous pipes.
///
/// The pipe class is a simple RAII model for anonymous pipes. It
/// provides a portable constructor that allocates a new %pipe and creates
/// a pipe object that owns the two file handles associated to it: the
/// read end and the write end.
///
class BOOST_PROCESSES_DECL pipe
{
    /// \brief The %pipe's read end.
    ///
    pipe_end read_end_;

    /// \brief The %pipe's write end.
    ///
    pipe_end write_end_;

public:
    /// \brief Creates a new %pipe.
    ///
    /// The default pipe constructor allocates a new anonymous %pipe
    /// and assigns its ownership to the created pipe object.
    ///
    /// \throw system_error If the anonymous %pipe creation fails.
    ///
    pipe();

    /// \brief Returns the %pipe's read end.
    ///
    /// \return A reference to the %pipe's read end file handle.
    ///
    pipe_end& rend()
    {
        return read_end_;
    }

    /// \brief Returns the %pipe's write end.
    ///
    /// \return A reference to the %pipe's write end file handle.
    ///
    pipe_end& wend()
    {
        return write_end_;
    }
};

} // namespace detail
} // namespace processes
} // namespace boost

#ifdef BOOST_PROCESSES_HEADER_ONLY
#include "boost/processes/detail/impl/pipe.hpp"
#endif

#endif
