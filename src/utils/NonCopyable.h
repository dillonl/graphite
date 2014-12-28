#ifndef GWIZ_NON_COPYABLE_H
#define GWIZ_NON_COPYABLE_H

namespace gwiz
{

//  Private copy constructor and copy assignment ensure classes derived from
//  class noncopyable cannot be copied.

//  Contributed by Dave Abrahams
//  This code has been repurposed from boost's noncopyable class

	namespace noncopyable_  // protection from unintended ADL
	{
		class noncopyable
		{
		protected:
			noncopyable() {}
			virtual ~noncopyable() {}
			noncopyable( const noncopyable& ) = delete;
			noncopyable& operator=( const noncopyable& ) = delete;
		};
	}

	typedef noncopyable_::noncopyable noncopyable;

} // namespace gwiz

#endif //GWIZ_NON_COPYABLE_H
