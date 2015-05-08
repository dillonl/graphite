#ifndef GWIZ_IALIGNMENT_H
#define GWIZ_IALIGNMENT_H

#include <boost/noncopyable.hpp>

#include <memory>
#include <unordered_map>

#include "core/utils/Types.h"

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
		virtual const std::string getID() { std::lock_guard< std::mutex > lock(this->m_mutex); return ""; }
		virtual const bool isFirstMate() { std::lock_guard< std::mutex > lock(this->m_mutex); return false;}
		virtual const bool isMapped() { std::lock_guard< std::mutex > lock(this->m_mutex); return false; }
		virtual const bool isReverseStrand() { std::lock_guard< std::mutex > lock(this->m_mutex); return false; }
		virtual const uint16_t getOriginalMapQuality() { std::lock_guard< std::mutex > lock(this->m_mutex); return 0; }
		void setMappingScoreAndVariants(uint32_t score, std::unordered_map< uint32_t, int32_t >& variantsIDAndScoreMap)
		{
			std::lock_guard< std::mutex > lock(this->m_mutex);
			if (this->m_score < score)
			{
				this->m_score = score;
				this->m_variants_id_and_score_map = variantsIDAndScoreMap;
			}
		}

		int32_t getVariantMappingScore(uint32_t variantID) { std::lock_guard< std::mutex > lock(this->m_mutex); return (this->m_variants_id_and_score_map.find(variantID) == this->m_variants_id_and_score_map.end()) ? -1 : this->m_variants_id_and_score_map[variantID]; }
		void setMapped(bool mapped) { std::lock_guard< std::mutex > lock(this->m_mutex); this->m_mapped = mapped; }
		bool getMapped() { std::lock_guard< std::mutex > lock(this->m_mutex); return this->m_mapped; }

	protected:
		std::mutex m_mutex;
		bool m_mapped;
		uint32_t m_score;
		std::unordered_map< uint32_t, int32_t > m_variants_id_and_score_map;
	};
}

#endif //GWIZ_IALIGNMENT_H
