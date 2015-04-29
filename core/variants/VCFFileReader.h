#ifndef GWIZ_VCFFILEREADER_H
#define GWIZ_VCFFILEREADER_H

#include "core/variants/IVariantList.h"
#include "core/utils/file/IFile.h"

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
	/* class VCFFileReader : public ASCIIFileReader, public IVariantList */
	class VCFFileReader : public IVariantList
	{
    public:
		typedef std::shared_ptr<VCFFileReader> SharedPtr;
		VCFFileReader(const std::string& path, Region::SharedPtr region);
		VCFFileReader(const std::string& path);
		~VCFFileReader();

		bool getNextVariant(Variant::SharedPtr& variant) override
		{
			/* const char* line = getNextLine(); */
			std::string line;
			if (m_file_ptr->getNextLine(line))
			{
				variant = variant->BuildVariant(line.c_str(), m_vcf_parser); // static function call to BuildVariant
				return true;
			}
			else
			{
				variant = NULL;
				return false;
			}
		}

		void Open();
		/* void Open(Region::SharedPtr region); */
		void rewind() override;
		IVariantList::SharedPtr getVariantsInRegion(Region::SharedPtr region) override { return nullptr; }
		size_t getCount() override { throw "Not implemented"; }

		void printToVCF(std::ostream& out) override;
	protected:
		void setFileReader(const std::string& path);
		IFile::SharedPtr m_file_ptr;
		Region::SharedPtr m_region_ptr;

	private:
		/* size_t getHeaderSize(); */
        void readHeader();
		position getPositionFromLine(const char* line);
		/* bool setRegionPositions(Region::SharedPtr regionPtr, const char* startLine, const char* endLine); */
		/* bool registerRegion(Region::SharedPtr region); */

		VariantParser< const char* > m_vcf_parser;

		std::mutex m_region_mutex;
	};
} // end namespace gwiz

#endif  //GWIZ_VCFFILEREADER_H
