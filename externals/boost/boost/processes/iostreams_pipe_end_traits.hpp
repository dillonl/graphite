// Boost.Process
//
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

#ifndef BOOST_PROCESSES_IOSTREAMS_PIPE_END_TRAITS_HPP_105121
#define BOOST_PROCESSES_IOSTREAMS_PIPE_END_TRAITS_HPP_105121

#include "boost/processes/pipe_end.hpp"
#include "boost/iostreams/categories.hpp"
#include "boost/iostreams/traits_fwd.hpp"

namespace boost {
namespace iostreams {

template<>
struct char_type_of<boost::processes::pipe_end>
{ typedef char type; };

template<>
struct category_of<boost::processes::pipe_end>
{ typedef bidirectional_device_tag type; };

} // namespace iostreams {
} // namespace boost {

#endif
