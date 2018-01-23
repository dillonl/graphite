#ifndef GRAPHITE_BAMFILEWRITER_H
#define GRAPHITE_BAMFILEWRITER_H

#include "core/util/Noncopyable.hpp"

#include "api/BamWriter.h"
#include "api/BamAlignment.h"
#include "api/SamHeader.h"

#include <memory>
#include <mutex>
#include <vector>
#include <unordered_map>

namespace graphite
{
	class BamFileWriter : private Noncopyable
	{
	public:
		typedef std::shared_ptr< BamFileWriter > SharedPtr;
		BamFileWriter(const std::string& path, BamTools::SamHeader header, BamTools::RefVector& refVec);
		~BamFileWriter();

        void writeAlignment(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, const std::string& ref);
		void close();

	private:
		bool m_open;
		std::string m_path;
		std::string m_tmp_path;
		BamTools::BamWriter m_tmp_bam_writer;
		BamTools::SamHeader m_header;
		BamTools::RefVector m_ref_vec;
		std::unordered_map< std::string, uint32_t > m_ref_map;
		std::unordered_map< std::string, uint32_t > m_ref_id_map; // stores the idx for the region for fast lookup
		std::vector< std::string > m_ref;
		std::mutex m_writer_lock;
		std::mutex m_file_open_lock;

	};
}
#endif //GRAPHITE_BAMFILEWRITER_H
