#ifndef GRAPHITE_THREADPOOL_HPP
#define GRAPHITE_THREADPOOL_HPP

#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>

#include "core/util/Noncopyable.hpp"

namespace graphite
{

    class ThreadPool : private Noncopyable
    {
	public:
		static ThreadPool* Instance()
		{
			static ThreadPool* s_threadpool = new ThreadPool(); // lazy initialization
			return s_threadpool;
		}

		uint32_t getThreadCount()
		{
			return m_thread_count;
		}

		void joinAll()
		{
			stop();
			start();
		}

		void start()
		{
			if (this->m_stopped)
			{
				{
					std::unique_lock<std::mutex> lock(m_tasks_mutex);
					this->m_stopped = false;
					this->m_workers.clear();
				}
				init();
			}
		}

		void stop()
		{
			{
				std::unique_lock<std::mutex> lock(m_tasks_mutex);
				this->m_stopped = true;
			}
			this->m_condition.notify_all();
			for (std::thread& worker : this->m_workers)
			{
				worker.join();
			}
		}

		template<class F, class... Args>
		auto enqueue(F&& funct, Args&&... args)
			-> std::shared_ptr< std::future< typename std::result_of< F(Args...) >::type > >
		{
			using return_type = typename std::result_of< F(Args...) >::type;

			auto task = std::make_shared< std::packaged_task< return_type() > >(
				std::bind(std::forward< F >(funct), std::forward< Args >(args)...)
				);
			std::shared_ptr< std::future< return_type > > res = std::make_shared< std::future< return_type > >(task->get_future());
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

		void setThreadCount(uint32_t threadCount)
		{
			this->m_thread_count = threadCount;
			stop();
			start();
		}

		int getTaskCount()
		{
			std::unique_lock< std::mutex > lock(this->m_tasks_mutex);
			return this->m_tasks.size();
		}

	private:

		ThreadPool() :
			m_thread_count(std::thread::hardware_concurrency() * 2),
			m_stopped(false)
		{
			init();
		}

		~ThreadPool()
		{
			stop();
		}

		void init()
		{
			for (size_t i = 0; i < this->m_thread_count; ++i)
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
		uint32_t m_thread_count;
    };

}

#endif //GRAPHITE_THREADPOOL_HPP
