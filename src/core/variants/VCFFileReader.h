#ifndef GWIZ_VCFFILEREADER_H
#define GWIZ_VCFFILEREADER_H

#include "core/utils/file/ASCIIFileReader.h"
#include "core/variants/IVariantList.h"

#include "core/variants/VCFParser.hpp"
#include "core/variants/ChromParser.hpp"
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
	class VCFFileReader : public ASCIIFileReader, public IVariantList
	{
    public:
		typedef std::shared_ptr<VCFFileReader> SharedPtr;
		VCFFileReader(const std::string& path, Region::SharedPtr region);
		VCFFileReader(const std::string& path);
		~VCFFileReader();

		bool getNextVariant(Variant::SharedPtr& variant) override
		{
			const char* line = getNextLine();
			if (line != NULL)
			{
				variant = variant->BuildVariant(line, m_vcf_parser); // static function call to BuildVariant
				return true;
			}
			else
			{
				variant = NULL;
				return false;
			}
		}

		void Open() override;
		void Open(Region::SharedPtr region);
	protected:


	private:
		size_t getHeaderSize();
        void readHeader();
		position getPositionFromLine(const char* line);
		bool setRegionPositions(Region::SharedPtr regionPtr, const char* startLine, const char* endLine);
		bool registerRegion(Region::SharedPtr region);

		VariantParser< const char* > m_vcf_parser;

		std::mutex m_region_mutex;
	};
} // end namespace gwiz

#endif  //GWIZ_VCFFILEREADER_H
