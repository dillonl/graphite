#ifndef GRAPHITE_IMAPPING_H
#define GRAPHITE_IMAPPING_H

#include "core/alignment/IAlignment.h"
#include "core/allele/IAllele.h"
#include "core/mapping/MappingAlignmentInfo.h"

#include <memory>
#include <vector>

#include <boost/noncopyable.hpp>

namespace graphite
{
	class IMapping : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IMapping > SharedPtr;
		IMapping() {}
		~IMapping() {}

		virtual MappingAlignmentInfo::SharedPtr getMappingAlignmentInfo(IAllele::SharedPtr allelePtr, std::shared_ptr< IAdjudicator > adjudicatorPtr) = 0;
		virtual int getMappingScore() = 0;
		virtual IAlignment::SharedPtr getAlignmentPtr() = 0;
		virtual std::vector< IAllele::SharedPtr > getAllelePtrs() = 0;
		virtual position getPosition() = 0;
		/* virtual MappingAlignment::SharedPtr getGSSWAlignmentPtrFromAllelePtr(IAllele::SharedPtr allelePtr) = 0; */
	};
}

#endif //GRAPHITE_IMAPPING_H
