#ifndef GRAPHITE_VCFFILEREADER_H
#define GRAPHITE_VCFFILEREADER_H

#include "core/variant/IVariantList.h"
#include "core/file/IFile.h"

#include "core/region/Region.h"
#include "core/util/SharedCreator.hpp"
#include "core/util/Noncopyable.hpp"

#include "Variant.h"

#include <list>
#include <tuple>
#include <map>
#include <thread>
#include <mutex>
#include <future>

namespace graphite
{
	class VCFFileReader : private Noncopyable
	{
    public:
		typedef std::shared_ptr<VCFFileReader> SharedPtr;
		typedef std::weak_ptr<VCFFileReader> WeakPtr;
		~VCFFileReader();

		inline static VCFFileReader::SharedPtr CreateVCFFileReader(const std::string& path, const std::vector< Sample::SharedPtr >& samplePtrs, IReference::SharedPtr referencePtr, uint32_t maxAllowedAlleleSize)
		{
			auto vcfFileReaderPtr = std::shared_ptr< VCFFileReader >(new VCFFileReader(path, samplePtrs, referencePtr, maxAllowedAlleleSize));
			return vcfFileReaderPtr;
		}

		std::string getFilePath();

		uint32_t getID() { return this->m_id; }
		std::vector< IVariant::SharedPtr > getVariantsInRegion(Region::SharedPtr regionPtr);

		static std::vector< Region::SharedPtr > GetAllRegionsInVCF(const std::string& vcfPath);

		VCFFileReader(const std::string& path, const std::vector< Sample::SharedPtr >& samplePtrs);
		VCFHeader::SharedPtr getVCFHeader() { return this->m_vcf_header; }

	protected:
		VCFFileReader(const std::string& path, const std::vector< Sample::SharedPtr >& samplePtrs, IReference::SharedPtr referencePtr, uint32_t maxAllowedAlleleSize);
	private:
		void Open();
        void readHeader();
		static position getPositionFromLine(const char* line);
		void setFileReader(const std::string& path);
		uint64_t findRegionStartPosition(Region::SharedPtr regionPtr);
		uint64_t getPositionFromFile(uint64_t seekPosition, uint64_t endSeekPosition, std::shared_ptr< std::atomic< bool > > posFound, Region::SharedPtr regionPtr, std::string path);

		VCFHeader::SharedPtr m_vcf_header;
		std::string m_path;
		VCFFileReader::WeakPtr m_this_wk_ptr;
		IFile::SharedPtr m_file_ptr;
		uint32_t m_id;
		IReference::SharedPtr m_reference_ptr;
		uint32_t m_max_allowed_allele_size;

		std::mutex m_region_mutex;

	};
} // end namespace graphite

#endif  //GRAPHITE_VCFFILEREADER_H
