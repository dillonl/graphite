#ifndef GRAPHITE_VCFREADER_H
#define GRAPHITE_VCFREADER_H

#include "core2/util/Noncopyable.hpp"
#include "core2/util/gzstream.h"
#include "core2/region/Region.h"
#include "core2/sample/Sample.h"
#include "Variant.h"
#include "VCFWriter.h"

#include <memory>

#include <istream>


namespace graphite
{
	class VCFReader : private Noncopyable
	{
	public:
		typedef std::shared_ptr< VCFReader > SharedPtr;
		VCFReader(const std::string& filename, std::vector< graphite::Sample::SharedPtr >& bamSamplePtrs, Region::SharedPtr regionPtr, VCFWriter::SharedPtr vcfWriter);
		~VCFReader();
		bool getNextVariants(std::vector< Variant::SharedPtr >& variantPtrs, uint32_t spacing);

	private:
		void openFile();
		void getRegionsFromVCF(std::vector< Region::SharedPtr >& regionPtrs);
		void processHeader(std::vector< graphite::Sample::SharedPtr >& bamSamplePtrs);
		Variant::SharedPtr getNextVariant();
		void setRegion(Region::SharedPtr regionPtr);

		std::string setSamplePtrs(const std::string& columnLine, std::vector< graphite::Sample::SharedPtr >& bamSamplePtrs);

		inline bool getNextLine(std::string& line)
		{
			return std::getline(*this->m_file_stream_ptr, line);
		}

		Variant::SharedPtr m_preloaded_variant;
		VCFWriter::SharedPtr m_vcf_writer;
		Region::SharedPtr m_region_ptr;
		std::string m_filename;
		std::shared_ptr< std::istream > m_file_stream_ptr;
        std::unordered_map< std::string, Sample::SharedPtr > m_sample_ptrs_map;
	};
}

#endif // GRAPHITE_VCFREADER_H
