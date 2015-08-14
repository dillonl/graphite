#ifndef GRAPHITE_VCFFILEREADER_H
#define GRAPHITE_VCFFILEREADER_H

#include "core/variant/IVariantList.h"
#include "core/file/IFile.h"

#include "core/parser/VCFParser.hpp"
#include "core/parser/ChromParser.hpp"
#include "core/region/Region.h"
#include "core/util/SharedCreator.hpp"

#include "Variant.h"

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

#include <list>
#include <tuple>
#include <map>
#include <thread>
#include <mutex>
#include <future>

namespace graphite
{
	class VCFFileReader : private boost::noncopyable
	{
    public:
		typedef std::shared_ptr<VCFFileReader> SharedPtr;
		typedef std::weak_ptr<VCFFileReader> WeakPtr;
		~VCFFileReader();

		inline static VCFFileReader::SharedPtr CreateVCFFileReader(const std::string& path)
		{
			auto vcfFileReaderPtr = std::shared_ptr< VCFFileReader >(new VCFFileReader(path));
			return vcfFileReaderPtr;
		}

		std::string getFilePath();

		uint32_t getID() { return this->m_id; }
		std::vector< IVariant::SharedPtr > getVariantsInRegion(Region::SharedPtr regionPtr);
	protected:
		VCFFileReader(const std::string& path);
	private:

		void Open();
        void readHeader();
		position getPositionFromLine(const char* line);
		void setFileReader(const std::string& path);

		std::string m_path;
		VCFFileReader::WeakPtr m_this_wk_ptr;
		IFile::SharedPtr m_file_ptr;
		VariantParser< const char* > m_vcf_parser;
		uint32_t m_id;

		std::mutex m_region_mutex;

	};
} // end namespace graphite

#endif  //GRAPHITE_VCFFILEREADER_H
