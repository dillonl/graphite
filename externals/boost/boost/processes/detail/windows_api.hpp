// Boost.Process
//
// Copyright 2008 Ilya Sokolov
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt.)

#ifndef BOOST_PROCESSES_WINDOWS_API_HPP_181027
#define BOOST_PROCESSES_WINDOWS_API_HPP_181027

#include "boost/config.hpp"

#if !defined(BOOST_WINDOWS)
#error this file is for windows only
#endif

#ifdef BOOST_USE_WINDOWS_H
#include <windows.h>
#endif

// workaround a bug in WinBase.h
#ifdef GetEnvironmentStrings
#undef GetEnvironmentStrings
#endif

namespace boost {
namespace processes {
namespace detail {
namespace windows_api {

#ifdef BOOST_USE_WINDOWS_H

    using ::DWORD;
    using ::HANDLE;
    using ::PROCESS_INFORMATION;
    using ::STARTUPINFOA;

    static const unsigned long CREATE_NEW_ = CREATE_NEW;
    static const unsigned long DUPLICATE_SAME_ACCESS_ =
        DUPLICATE_SAME_ACCESS;
    static const unsigned long ERROR_BROKEN_PIPE_ = ERROR_BROKEN_PIPE;
    static const unsigned long FILE_ATTRIBUTE_NORMAL_ =
        FILE_ATTRIBUTE_NORMAL;
    static const unsigned long GENERIC_READ_ = GENERIC_READ;
    static const unsigned long GENERIC_WRITE_ = GENERIC_WRITE;
    static const unsigned long HANDLE_FLAG_INHERIT_ = HANDLE_FLAG_INHERIT;
    static void* const INVALID_HANDLE_VALUE_ = INVALID_HANDLE_VALUE;
    static const unsigned long MAX_PATH_ = MAX_PATH;
    static const unsigned long OPEN_EXISTING_ = OPEN_EXISTING;
    static const unsigned long STD_ERROR_HANDLE_ = STD_ERROR_HANDLE;
    static const unsigned long STD_INPUT_HANDLE_ = STD_INPUT_HANDLE;
    static const unsigned long STD_OUTPUT_HANDLE_ = STD_OUTPUT_HANDLE;
    static const unsigned long STARTF_USESTDHANDLES_ = STARTF_USESTDHANDLES;
    static const unsigned long STILL_ACTIVE_ = STILL_ACTIVE;

    using ::CloseHandle;
    using ::CreateFileA;
    using ::CreatePipe;
    using ::CreateProcessA;
    using ::DuplicateHandle;
    using ::FreeEnvironmentStringsA;
    using ::GetCurrentDirectoryA;
    using ::GetCurrentProcess;
    using ::GetCurrentProcessId;
    using ::GetEnvironmentStrings;
    using ::GetExitCodeProcess;
    using ::GetSystemDirectoryA;
    using ::GetLastError;
    using ::GetStdHandle;
    using ::ReadFile;
    using ::SearchPathA;
    using ::SetHandleInformation;
    using ::Sleep;
    using ::TerminateProcess;
    using ::WaitForSingleObject;
    using ::WriteFile;
	using ::OpenProcess; //added by rtao, 01/25/2010
	using ::GetProcessId;
	

	static const unsigned long PROCESS_ALL_ACCESS_ = PROCESS_ALL_ACCESS;
	

#else // #ifdef BOOST_USE_WINDOWS_H

    typedef unsigned long DWORD;
    typedef void* HANDLE;

    struct PROCESS_INFORMATION
    {
        void* hProcess;
        void* hThread;
        unsigned int dwProcessId;
        unsigned int dwThreadId;
    };

    struct STARTUPINFOA
    {
        unsigned long cb;
        char* lpReserved;
        char* lpDesktop;
        char* lpTitle;
        unsigned long dwX;
        unsigned long dwY;
        unsigned long dwXSize;
        unsigned long dwYSize;
        unsigned long dwXCountChars;
        unsigned long dwYCountChars;
        unsigned long dwFillAttribute;
        unsigned long dwFlags;
        unsigned short wShowWindow;
        unsigned short cbReserved2;
        unsigned char* lpReserved2;
        void* hStdInput;
        void* hStdOutput;
        void* hStdError;
    };

    static const unsigned long CREATE_NEW_ = 1;
    static const unsigned long DUPLICATE_SAME_ACCESS_ = 2;
    static const unsigned long ERROR_BROKEN_PIPE_ = 109;
    static const unsigned long FILE_ATTRIBUTE_NORMAL_ = 0x00000080ul;
    static const unsigned long GENERIC_READ_ = 0x80000000ul;
    static const unsigned long GENERIC_WRITE_ = 0x40000000ul;
    static const unsigned long HANDLE_FLAG_INHERIT_ = 1;
    static void* const INVALID_HANDLE_VALUE_ = (void*)-1;
    static const unsigned long MAX_PATH_ = 260;
    static const unsigned long OPEN_EXISTING_ = 3;
    static const unsigned long STD_ERROR_HANDLE_ = (unsigned long)-12;
    static const unsigned long STD_INPUT_HANDLE_ = (unsigned long)-10;
    static const unsigned long STD_OUTPUT_HANDLE_ = (unsigned long)-11;
    static const unsigned long STARTF_USESTDHANDLES_ = 0x100ul;
    static const unsigned long STILL_ACTIVE_ = 0x103ul;

    extern "C" __declspec(dllimport) int __stdcall
    CloseHandle(void*);

    extern "C" __declspec(dllimport) void* __stdcall
    CreateFileA(const char*, unsigned long, unsigned long,
                void*, unsigned long, unsigned long, void*);

    extern "C" __declspec(dllimport) int __stdcall
    CreatePipe(void*, void*, void*, unsigned long);

    extern "C" __declspec(dllimport) int __stdcall
    CreateProcessA(const char*, char*, void*,
        void*, int, unsigned long, void*, const char*,
        STARTUPINFOA*, PROCESS_INFORMATION*);

    extern "C" __declspec(dllimport) int __stdcall
    DuplicateHandle(void*, void*, void*, void*,
                    unsigned int, int, unsigned int);

    extern "C" __declspec(dllimport) int __stdcall
    FreeEnvironmentStringsA(char*);

    extern "C" __declspec(dllimport) unsigned int __stdcall
    GetCurrentDirectoryA(unsigned int, char*);

    extern "C" __declspec(dllimport) void* __stdcall
    GetCurrentProcess();

    extern "C" __declspec(dllimport) unsigned long __stdcall
    GetCurrentProcessId();

    extern "C" __declspec(dllimport) char* __stdcall
    GetEnvironmentStrings();

    extern "C" __declspec(dllimport) int __stdcall
    GetExitCodeProcess(void*, unsigned long*);

    extern "C" __declspec(dllimport) unsigned int __stdcall
    GetSystemDirectoryA(char*, unsigned int);

    extern "C" __declspec(dllimport) unsigned long __stdcall
    GetLastError();

    extern "C" __declspec(dllimport) void* __stdcall
    GetStdHandle(unsigned long);

    extern "C" __declspec(dllimport) int __stdcall
    ReadFile(void*, const void*, unsigned long, unsigned long*, void*);

    extern "C" __declspec(dllimport) unsigned long __stdcall
    SearchPathA(const char*, const char*, const char*,
                unsigned long, char*, char**);

    extern "C" __declspec(dllimport) int __stdcall
    SetHandleInformation(void*, unsigned long, unsigned long);

    extern "C" __declspec(dllimport) void __stdcall
    Sleep(unsigned long);

    extern "C" __declspec(dllimport) int __stdcall
    TerminateProcess(void*, unsigned int);

    extern "C" __declspec(dllimport) unsigned long __stdcall
    WaitForSingleObject(void*, unsigned long);

    extern "C" __declspec(dllimport) int __stdcall
    WriteFile(void*, const void*, unsigned long, unsigned long*, void*);
	

	/*This for terminate process tree, added by rtao*/
	extern "C" __declspec(dllimport) void* __stdcall
    OpenProcess(unsigned long, bool, unsigned long);
    
    extern "C" __declspec(dllimport) unsigned long __stdcall
    GetProcessId(void* );

	static const unsigned long PROCESS_ALL_ACCESS_ = 0x1fffff;


#endif // #else // #ifdef BOOST_USE_WINDOWS_H

} // namespace windows_api {
} // namespace detail {
} // namespace processes {
} // namespace boost {

#endif
