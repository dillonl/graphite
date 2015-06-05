#ifndef GWIZ_VCFFILEREADER_H
#define GWIZ_VCFFILEREADER_H

#include "core/variant/IVariantList.h"
#include "core/file/IFile.h"

#include "core/parser/VCFParser.hpp"
#include "core/parser/ChromParser.hpp"
#include "core/region/Region.h"

#include "Variant.h"

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

#include <list>
#include <tuple>
#include <map>
#include <thread>
#include <mutex>
#include <future>

namespace gwiz
{
	class VCFFileReader : private boost::noncopyable
	{
    public:
		typedef std::shared_ptr<VCFFileReader> SharedPtr;
		VCFFileReader(const std::string& path);
		~VCFFileReader();

		std::vector< IVariant::SharedPtr > getVariantsInRegion(Region::SharedPtr regionPtr);
	private:
		void Open();
        void readHeader();
		position getPositionFromLine(const char* line);
		void setFileReader(const std::string& path);

		IFile::SharedPtr m_file_ptr;
		VariantParser< const char* > m_vcf_parser;

		std::mutex m_region_mutex;
	};
} // end namespace gwiz

#endif  //GWIZ_VCFFILEREADER_H
