#ifndef GRAPHITE_ALLELEMETADATA_H
#define GRAPHITE_ALLELEMETADATA_H

#include <boost/noncopyable.hpp>

#include <memory>

namespace graphite
{
	class AlleleMetaData : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< AlleleMetaData > SharedPtr;
	    AlleleMetaData(uint16_t paddingPrefix, uint16_t paddingSuffix) :
		    m_padding_prefix(paddingPrefix), m_padding_suffix(paddingSuffix)
		{
		}
		~AlleleMetaData() {}

		uint16_t getPaddingPrefix() { return this->m_padding_prefix; }
		uint16_t getPaddingSuffix() { return this->m_padding_suffix; }

	private:
		uint16_t m_padding_prefix;
		uint16_t m_padding_suffix;
	};
}

#endif //GRAPHITE_ALLELEMETADATA_H
