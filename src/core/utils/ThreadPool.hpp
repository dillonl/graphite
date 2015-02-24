#ifndef GWIZ_THREADPOOL_HPP
#define GWIZ_THREADPOOL_HPP

#include <boost/noncopyable.hpp>
#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>

#include <memory>

namespace gwiz
{
	class ThreadPool : private boost::noncopyable
	{
	public:
		static void PostJob(boost::function< void() > funct)
		{
			static ThreadPool s_threadpool;
			s_threadpool.postThreadPoolJob(funct);
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
			{
				m_threadgroup.create_thread(boost::bind(&boost::asio::io_service::run, &m_io_service));
			}
		}

		void postThreadPoolJob(boost::function< void() > funct)
		{
			this->m_io_service.post(funct);
		}

		static ThreadPool s_threadpool;

		boost::asio::io_service m_io_service;
		boost::thread_group m_threadgroup;
		boost::asio::io_service::work m_work;
	};
}

#endif //GWIZ_THREADPOOL_HPP
