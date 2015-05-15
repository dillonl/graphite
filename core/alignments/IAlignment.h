#ifndef GWIZ_IALIGNMENT_H
#define GWIZ_IALIGNMENT_H

#include <boost/noncopyable.hpp>

#include <memory>
#include <unordered_map>
#include <mutex>

#include "core/utils/Types.h"
#include "core/variants/Variant.h"

namespace gwiz
{

	class IAlignment : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IAlignment > SharedPtr;

	    IAlignment() : m_mapped(false) {}
		virtual ~IAlignment() {}

		virtual const char* getSequence() = 0;
		virtual const position getPosition() = 0;
		virtual const size_t getLength() = 0;
		virtual const std::string getID() { return ""; }
		virtual const bool isFirstMate() { return false;}
		virtual const bool isMapped() { return false; }
		virtual const bool isReverseStrand() { return false; }
		virtual const uint16_t getOriginalMapQuality() { return 0; }

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
	};
}

#endif //GWIZ_IALIGNMENT_H
