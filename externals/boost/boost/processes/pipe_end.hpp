// Boost.Process
//
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

#ifndef BOOST_PROCESSES_PIPE_END_HPP_114037
#define BOOST_PROCESSES_PIPE_END_HPP_114037

#include "boost/processes/config.hpp"
#include "boost/processes/detail/handle.hpp"
#include <ios>

namespace boost {
namespace processes {

class BOOST_PROCESSES_DECL pipe_end: public detail::handle {
public:
    pipe_end() {} // never throws

    explicit pipe_end(native_type native)
    : detail::handle(native)
    {}

    std::streamsize read(char* s, std::streamsize n);
    std::streamsize write(const char* s, std::streamsize n);
};

} // namespace processes {
} // namespace boost {

#ifdef BOOST_PROCESSES_HEADER_ONLY
#include "boost/processes/detail/impl/pipe_end.hpp"
#endif

#endif
