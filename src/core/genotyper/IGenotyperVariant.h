#ifndef GWIZ_GENOTYPERVARIANT_H
#define GWIZ_GENOTYPERVARIANT_H

#include <memory>

#include <boost/noncopyable.hpp>

#include "GenotyperAllele.hpp"

namespace gwiz
{
	class IGenotyperVariant : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IGenotyperVariant > SharedPtr;
	IGenotyperVariant(position pos) :
		    m_position(pos)
		{
		}

		virtual ~IGenotyperVariant()
		{
		}

		void addAllele(GenotyperAllele::SharedPtr allele)
		{
			m_alleles.push_back(allele);
		}

		std::list< GenotyperAllele::SharedPtr > getGenotyperAlleles() { return this->m_alleles; }

		position getPosition() { return this->m_position; }

	protected:
		position m_position;
		std::list< GenotyperAllele::SharedPtr > m_alleles;
	};
}

#endif //GWIZ_GENOTYPERVARIANT_H
