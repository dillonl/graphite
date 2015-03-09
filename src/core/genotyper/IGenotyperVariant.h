#ifndef GWIZ_GENOTYPERVARIANT_H
#define GWIZ_GENOTYPERVARIANT_H

#include <memory>

#include <boost/noncopyable.hpp>

namespace gwiz
{
	class IGenotyperVariant : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IGenotyperVariant > SharedPtr;
		IGenotyperVariant()
		{
		}

		virtual ~IGenotyperVariant()
		{
		}

		void addAllele(IGenotyperAllele::SharedPtr allele) = 0;

	protected:
		std::vector< IGenotypeAllele::SharedPtr > m_alleles;
	};
}

#endif //GWIZ_GENOTYPERVARIANT_H