#ifndef GRAPHITE_MAPPINGMANAGER_H
#define GRAPHITE_MAPPINGMANAGER_H

#include "IMapping.h"
#include "core/adjudicator/IAdjudicator.h"

#include <unordered_map>
#include <mutex>

#include <boost/noncopyable.hpp>

namespace graphite
{
	class MappingManager : private boost::noncopyable
	{
	public:
		static MappingManager* Instance();

		/*
		 * Checks the mapping's mappingScore and only
		 * adds the alignmentPtr if the passed in
		 * mappingPtr's mapping score is larger.
		 */
		void registerMapping(IMapping::SharedPtr mappingPtr);
		void evaluateAlignmentMappings(IAdjudicator::SharedPtr adjudicatorPtr);
		void clearRegisteredMappings();
	private:
		MappingManager();
		~MappingManager();

		static MappingManager* s_instance;
		std::mutex m_alignment_mapping_map_mutex;
		std::unordered_map< IAlignment::SharedPtr, IMapping::SharedPtr > m_alignment_mapping_map;
		std::vector< IMapping::SharedPtr > m_mappings;
	};
}

#endif //GRAPHITE_MAPPINGMANAGER_H
