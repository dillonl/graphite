#ifndef GRAPHITE_SRC_UTILS_NONCOPYABLE_HPP
#define GRAPHITE_SRC_UTILS_NONCOPYABLE_HPP

namespace graphite
{
	class Noncopyable
	{
	public:
		Noncopyable( const Noncopyable& noncopyable) = delete;
		Noncopyable& operator=( const Noncopyable& noncopyable) = delete;

		Noncopyable() = default;

	};
}

#endif //GRAPHITE_SRC_UTILS_NONCOPYABLE_HPP
