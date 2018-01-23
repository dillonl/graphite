#ifndef GRAPHITE_IALIGNMENT_H
#define GRAPHITE_IALIGNMENT_H

#include "core/util/Types.h"
#include "core/util/Noncopyable.hpp"
#include "core/allele/IAllele.h"
#include "core/sample/Sample.h"

#include "api/BamAlignment.h"

#include <memory>
#include <unordered_map>
#include <mutex>
#include <vector>
#include <string>

namespace graphite
{

	class IMapping;
	class IAlignment : private Noncopyable
	{
	public:
		typedef std::shared_ptr< IAlignment > SharedPtr;

	    IAlignment() : m_mapping_mutex(new std::recursive_mutex()) {}
		virtual ~IAlignment() { delete this->m_mapping_mutex; }

		virtual const char* getSequence() = 0;
		virtual const position getPosition() = 0;
		virtual const size_t getLength() = 0;
		virtual const std::string getID() { return ""; }
		virtual const bool isFirstMate() { return false;}
		virtual const bool isMapped() { return false; }
		virtual const bool isReverseStrand() { return false; }
		virtual const bool isDuplicate() { return false; }
		virtual const uint16_t getOriginalMapQuality() { return 0; }
		virtual std::vector< std::shared_ptr< IMapping > > getMappingPtrs() { return this->m_mapping_ptrs; }
		virtual void addMapping(std::shared_ptr< IMapping > mappingPtr)
		{
			std::lock_guard< std::recursive_mutex > r_lock(*this->m_mapping_mutex);
			m_mapping_ptrs.emplace_back(mappingPtr);
		}
		std::recursive_mutex* getMappingMutex() { return this->m_mapping_mutex; }
		const Sample::SharedPtr getSample() { return m_sample_ptr; }
		virtual const std::shared_ptr< BamTools::BamAlignment > getBamToolsAlignmentPtr() = 0;


	protected:
		std::mutex m_mutex;
		std::unordered_map< uint32_t, std::string > m_mapped_variants_information;
		std::vector< std::shared_ptr< IMapping > > m_mapping_ptrs;
		std::recursive_mutex* m_mapping_mutex;
		Sample::SharedPtr m_sample_ptr;
	};
}

#endif //GRAPHITE_IALIGNMENT_H
