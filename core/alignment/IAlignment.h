#ifndef GRAPHITE_IALIGNMENT_H
#define GRAPHITE_IALIGNMENT_H

#include <boost/noncopyable.hpp>

#include <memory>
#include <unordered_map>
#include <mutex>
#include <vector>
#include <string>

#include "core/util/Types.h"

namespace graphite
{

	class IMapping;
	class IAlignment : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IAlignment > SharedPtr;

	    IAlignment() : m_mapped(false), m_score(0), m_mapping_mutex(new std::recursive_mutex()) {}
		virtual ~IAlignment() { delete this->m_mapping_mutex; }

		virtual const char* getSequence() = 0;
		virtual const position getPosition() = 0;
		virtual const size_t getLength() = 0;
		virtual const std::string getID() { return ""; }
		virtual const bool isFirstMate() { return false;}
		virtual const bool isMapped() { return false; }
		virtual const bool isReverseStrand() { return false; }
		virtual const uint16_t getOriginalMapQuality() { return 0; }
		virtual std::weak_ptr< IMapping > getMapping() { return this->m_mapping_wptr; }
		virtual void setMapping(std::shared_ptr< IMapping > mappingPtr)
		{
			std::lock_guard< std::recursive_mutex > r_lock(*this->m_mapping_mutex);
			this->m_mapping_wptr = mappingPtr;
		}
		std::recursive_mutex* getMappingMutex() { return this->m_mapping_mutex; }

		/*
		 * the variantInformation is the variantID and the allele that was matched
		 */
		void setMappingInformation(uint32_t score, std::vector< std::tuple< uint32_t, std::string > >& variantInformation)
		{
			std::lock_guard< std::mutex > lock(this->m_mutex);
			if (score < this->m_score)
			{
				return;
			}
			this->m_mapped = true;
			this->m_score = score;
			this->m_mapped_variants_information.clear();
			for (const auto& variantInfo : variantInformation)
			{
				this->m_mapped_variants_information[std::get< 0 >(variantInfo)] = std::get< 1 >(variantInfo);
			}
		}

		int32_t getVariantMappingScore(const uint32_t variantID)
		{
			std::lock_guard< std::mutex > lock(this->m_mutex);
			auto variantIDIter = this->m_mapped_variants_information.find(variantID);
			int32_t score = (variantIDIter == this->m_mapped_variants_information.end()) ? -1 : this->m_score;
			return score;
		}

		std::string getVariantAllele(uint32_t variantID)
		{
			std::lock_guard< std::mutex > lock(this->m_mutex);
			auto variantIDIter = this->m_mapped_variants_information.find(variantID);
			std::string variantAllele = (variantIDIter == this->m_mapped_variants_information.end()) ? "" : variantIDIter->second;
			return variantAllele;
		}

		void setMapped(bool mapped) { this->m_mapped = mapped; }
		bool getMapped() { return this->m_mapped; }

	protected:
		std::mutex m_mutex;
		bool m_mapped;
		uint32_t m_score;
		std::unordered_map< uint32_t, std::string > m_mapped_variants_information;
		std::weak_ptr< IMapping > m_mapping_wptr;
		std::recursive_mutex* m_mapping_mutex;
	};
}

#endif //GRAPHITE_IALIGNMENT_H
