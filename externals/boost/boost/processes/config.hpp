// Boost.Process
//
// Copyright (c) 2006 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Copyright John Maddock 1998
//     We are using quotes from
//     'Guidelines for Authors of Boost Libraries Containing Separate Source',
//     http://www.boost.org/development/separate_compilation.html
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

/// \file boost/process/config.hpp
///
/// Defines macros that are used by the library's code to determine the
/// operating system it is running under and the features it supports.

#ifndef BOOST_PROCESSES_CONFIG_HPP_170912
#define BOOST_PROCESSES_CONFIG_HPP_170912

#include "boost/config.hpp"
#include "boost/system/config.hpp"

#ifdef BOOST_PROCESSES_HEADER_ONLY
    #define BOOST_PROCESSES_INLINE_IF_HEADER_ONLY inline
#else
    #define BOOST_PROCESSES_INLINE_IF_HEADER_ONLY
    #ifdef BOOST_HAS_DECLSPEC // defined in config system
        // we need to import/export our code only if the user has specifically
        // asked for it by defining either BOOST_ALL_DYN_LINK if they want all
        // boost libraries to be dynamically linked, or BOOST_PROCESS_DYN_LINK
        // if they want just this one to be dynamically liked
        #if defined(BOOST_ALL_DYN_LINK) || defined(BOOST_PROCESSES_DYN_LINK)
            // export if this is our own source, otherwise import
            #ifdef BOOST_PROCESSES_SOURCE
                #define BOOST_PROCESSES_DECL __declspec(dllexport)
            #else
                #define BOOST_PROCESSES_DECL __declspec(dllimport)
            #endif  // BOOST_WHATEVER_SOURCE
        #endif  // DYN_LINK
    #endif  // BOOST_HAS_DECLSPEC

    // Automatically link to the correct build variant where possible
    #if !defined(BOOST_ALL_NO_LIB) \
        && !defined(BOOST_PROCESSES_NO_LIB) && !defined(BOOST_PROCESSES_SOURCE)
        // Set the name of our library,
        // this will get undef'ed by auto_link.hpp once it's done with it
        #define BOOST_LIB_NAME boost_process
        // If we're importing code from a dll, then tell auto_link.hpp about it
        #if defined(BOOST_ALL_DYN_LINK) || defined(BOOST_PROCESSES_DYN_LINK)
            #define BOOST_DYN_LINK
        #endif
        // And include the header that does the work
        #include "boost/config/auto_link.hpp"
    #endif // auto-linking disabled

#endif // BOOST_PROCESSES_HEADER_ONLY

// if BOOST_WHATEVER_DECL isn't defined yet define it now:
#ifndef BOOST_PROCESSES_DECL
#define BOOST_PROCESSES_DECL
#endif

#endif
