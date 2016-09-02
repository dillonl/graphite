#ifndef GRAPHITE_GENOTYPERVARIANT_H
#define GRAPHITE_GENOTYPERVARIANT_H

#include <memory>

#include "GenotyperAllele.hpp"
#include "core/util/Noncopyable.hpp"

namespace graphite
{
	class IGenotyperVariant : private Noncopyable
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

#endif //GRAPHITE_GENOTYPERVARIANT_H
