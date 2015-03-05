#ifndef GWIZ_THREADPOOL_HPP
#define GWIZ_THREADPOOL_HPP

#include <boost/noncopyable.hpp>
#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>

#include <thread>
#include <memory>

namespace gwiz
{
	class ThreadPool : private boost::noncopyable
	{
	public:

		static ThreadPool* Instance()
		{
			static ThreadPool* s_threadpool = NULL; // lazy initialization
			if (s_threadpool == NULL)
			{
				s_threadpool = new ThreadPool();
			}
			return s_threadpool;
		}

		void postJob(boost::function< void() > funct)
		{
			this->m_io_service.dispatch(std::move(funct));
		}

		void joinAll()
		{
			m_io_service.stop();
			this->m_threadgroup.join_all();
		}

		void startIOService()
		{
			m_io_service.run();
		}

	private:
		ThreadPool() : m_work(m_io_service)
		{
			initializeThreadPool();
		}

		~ThreadPool()
		{
			m_io_service.stop();
			m_threadgroup.join_all();
		}

		void initializeThreadPool()
		{
			for(uint32_t i = 0; i < std::thread::hardware_concurrency(); ++i)
			// for(uint32_t i = 0; i < 1; ++i)
			{
				m_threadgroup.create_thread(boost::bind(&boost::asio::io_service::run, &m_io_service));
			}
			// m_io_service.run();
		}

		static ThreadPool* s_threadpool;

		boost::asio::io_service m_io_service;
		boost::thread_group m_threadgroup;
		boost::asio::io_service::work m_work;
	};
}

#endif //GWIZ_THREADPOOL_HPP
