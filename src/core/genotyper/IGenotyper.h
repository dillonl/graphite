#ifndef GWIZ_GENOTYPER_H
#define GWIZ_GENOTYPER_H

#include <memory>

#include <boost/noncopyable.hpp>

namespace gwiz
{
	class IGenotyper : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IGenotyper > SharedPtr;
	    IGenotyper(IGenotyperVariant::SharedPtr genotyperVariant) :
		    m_genotyper_variant(genotyperVariant)
		{
		}

		virtual ~IGenotyper()
		{
		}

	protected:
		std::vector< IGenotyperVariant::SharedPtr > m_genotyper_variant;
	};
}

#endif //GWIZ_GENOTYPER_H
