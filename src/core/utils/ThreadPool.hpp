#ifndef GWIZ_THREADPOOL_HPP
#define GWIZ_THREADPOOL_HPP

#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>
#include <boost/bind.hpp>


/*
#include <boost/noncopyable.hpp>
#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>

#include <thread>
#include <memory>
*/

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

		// void postJob(boost::function< void() > funct)
		// {
		// }

		void joinAll()
		{
			this->m_condition.notify_all();
			for (std::thread& worker : this->m_workers)
			{
				worker.join();
			}
		}

		// template< class F, class... Args >
		// auto enqueue(F&& funct, Args&& args)->std::future< typename std::result_of< F(Args...) >::type >
		template<class F, class... Args>
		auto enqueue(F&& funct, Args&&... args)
			-> std::future<typename std::result_of<F(Args...)>::type>
		{
			using return_type = typename std::result_of< F(Args...) >::type;

			auto task = std::make_shared< std::packaged_task< return_type() > >(
				std::bind(std::forward< F >(funct), std::forward< Args >(args)...)
				);
			std::future< return_type > res = task->get_future();
			{
				std::unique_lock< std::mutex > lock(this->m_tasks_mutex);
				if (this->m_stopped)
				{
					throw std::runtime_error("enqueue on stopped ThreadPool");
				}
				this->m_tasks.emplace([task](){(*task)();});
			}
			this->m_condition.notify_one();
			return res;
		}

	private:

		ThreadPool() :
			m_stopped(false)
		{
			init();
		}

		~ThreadPool()
		{
			{
				std::unique_lock< std::mutex > lock(this->m_tasks_mutex);
				this->m_stopped = true;
			}
			joinAll();
		}

		void init()
		{
			size_t numberOfThreads = std::thread::hardware_concurrency();
			for (size_t i = 0; i < numberOfThreads; ++i)
			{
				this->m_workers.emplace_back(
					[this]
					{
						for (;;)
						{
							std::function< void() > task;
							{
								std::unique_lock< std::mutex > lock(this->m_tasks_mutex);
								this->m_condition.wait(lock,
													   [this]{ return this->m_stopped || !this->m_tasks.empty(); });
								if (this->m_stopped && this->m_tasks.empty()) { return; }
								task = std::move(this->m_tasks.front());
								this->m_tasks.pop();
							}
							task();
						}
					}
					);
			}
		}

		static ThreadPool* s_threadpool;
		std::vector< std::thread > m_workers;
		std::queue< std::function< void() > > m_tasks;

		std::mutex m_tasks_mutex;
		std::condition_variable m_condition;
		bool m_stopped;
    };
	/*
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
	*/
}

#endif //GWIZ_THREADPOOL_HPP
