// Boost.Process
//
// Copyright (c) 2006, 2007 Julio M. Merino Vidal.
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

#ifndef BOOST_PROCESSES_CHILD_HPP_145715
#define BOOST_PROCESSES_CHILD_HPP_145715

#include "boost/processes/child.hpp"

namespace boost {
namespace processes {

BOOST_PROCESSES_INLINE_IF_HEADER_ONLY
const
status child::wait()
{
#if defined(BOOST_POSIX_API)
    int s;
    if (::waitpid(pid_, &s, 0) == -1)
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("waitpid(2) failed");
    return status(s);
#elif defined(BOOST_WINDOWS_API)
    using namespace detail::windows_api;
    DWORD code;
    // XXX This loop should go away in favour of a passive wait.
    do {
        GetExitCodeProcess(handle_.get(), &code);
        Sleep(500);
    } while (code == STILL_ACTIVE_);
    WaitForSingleObject(handle_.get(), 0);
    return status(code);
#endif // #elif defined(BOOST_WINDOWS_API)
}

BOOST_PROCESSES_INLINE_IF_HEADER_ONLY
void
child::terminate(bool force) const
{
#if defined(BOOST_POSIX_API)
    if (::kill(pid_, force ? SIGKILL : SIGTERM) == -1)
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("kill(2) failed");
#elif defined(BOOST_WINDOWS_API)
    if (detail::windows_api::TerminateProcess(handle_.get(), EXIT_FAILURE) == 0)
        BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("TerminateProcess() failed");
#endif
}

BOOST_PROCESSES_INLINE_IF_HEADER_ONLY
bool
child::terminate_process_tree() const
{
#if defined(BOOST_POSIX_API)
	pid_t pid = pid_;
	if (::killpg(getpgid(pid), SIGTERM) == -1)
	{
		BOOST_PROCESSES_THROW_LAST_SYSTEM_ERROR("kill(2) failed");
		return false;
	}
	else
	{
		std::cout << "pid:" << pid << " killed!!!" << std::endl;
		return true;
	}

#elif defined(BOOST_WINDOWS_API)
	detail::windows_api::DWORD processId = detail::windows_api::GetProcessId(handle_.get());
	std::string pid_str =  boost::lexical_cast<std::string>(processId);
			
	std::string cmd = "taskkill /F /T /PID " + pid_str;
	::system( cmd.c_str() );
	return true;
#endif
}


#if defined(BOOST_POSIX_API)

BOOST_PROCESSES_INLINE_IF_HEADER_ONLY pipe_end
child::get_pipe_end(int desc) const
{
    pe_map::const_iterator iter = pes_.find(desc);
    BOOST_ASSERT(iter != pes_.end());
    return iter->second;
}

#endif // #if defined(BOOST_POSIX_API)

} // namespace processes
} // namespace boost

#endif
