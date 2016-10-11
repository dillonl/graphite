#ifndef GRAPHITE_HTSLIBALIGNMENT_H
#define GRAPHITE_HTSLIBALIGNMENT_H

#include "IAlignment.h"

namespace graphite
{
	class HTSLibAlignment : public IAlignment
	{
	public:
		typedef std::shared_ptr< HTSLibAlignment > SharedPtr;
	    HTSLibAlignment();
		~HTSLibAlignment();

		const char* getSequence() override;
		const position getPosition() override { return m_position; }
		const size_t getLength() override;
		const bool isReverseStrand() override { return m_reverse_strand; }
		std::recursive_mutex* getMappingMutex() { return this->m_mapping_mutex; }
		const std::shared_ptr< Sample > getSample() { return m_sample_ptr; }

		const std::string getID() { return m_id; }
		const void setSequence(char* seq, uint32_t len) override;
		const void setPosition(position pos) { m_position = pos; }
		const void setSample(std::shared_ptr< Sample > samplePtr) { m_sample_ptr = samplePtr; }
		const void setName(char* name, bool isFirstMate)
		{
				m_first_mate = isFirstMate;
				m_id = std::string(name) + std::to_string(isFirstMate);
		}
		const void setDuplicate(bool isDup) { m_duplicate = isDup; }
		const void setIsReverseStrand(bool isReverse) { m_reverse_strand = isReverse; }
		const void setRefID(uint32_t refID) { m_ref_id = refID; }
		const void setFilePosition(uint64_t filePos) { m_file_position = filePos; }
		const void removeSequence() override;
		const void incrementReferenceCount() override;

	private:
		std::mutex m_sequence_mutex;
		uint16_t m_sequence_counter;
		std::string m_sequence_string;
		uint32_t m_ref_id;
		position m_position;
		std::string m_id;
		bool m_first_mate;
		bool m_mapped;
		bool m_reverse_strand;
		bool m_duplicate;
		uint16_t m_original_map_quality;
		uint64_t m_file_position;
	};
}

#endif //GRAPHITE_HTSLIBALIGNMENT_H
