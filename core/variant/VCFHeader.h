#ifndef GRAPHITE_VCFHEADER_H
#define GRAPHITE_VCFHEADER_H

#include "IHeader.h"
#include "core/util/Types.h"

#include <memory>
#include <vector>
#include <string>
#include <unordered_set>

namespace graphite
{
	class IAlignmentReader;
	class VCFHeader : public IHeader
	{
	public:
		typedef std::shared_ptr< VCFHeader > SharedPtr;
		VCFHeader(const std::vector< std::string >& vcfHeaderLines);
		~VCFHeader();

		std::string getHeader() override;
		void registerReferencePath(const std::string& referencePath);
		void registerActiveSample(std::shared_ptr< Sample > samplePtr) override;
		std::vector< std::string > getSampleNames() override { return m_sample_names; }
		int32_t getColumnPosition(const std::string& columnTitle) override;
		void setColumn(const std::string& column) override;
		std::vector< std::string > getColumnNames() override;
		bool isActiveSampleColumnName(const std::string& headerName) override;
	private:
		void setColumns(const std::string& headerString);
		std::string getColumnsString();
		void addFormatToHeader();

		std::string m_reference_path;
		std::string m_columns_line;
		std::vector< std::string > m_lines;
		std::vector< std::string > m_columns;
	    std::vector< std::string > m_sample_names;
		std::vector< AlleleCountType > m_allele_count_type_registry;
		std::unordered_set< std::string > m_active_sample_name;
	};
}

#endif //GRAPHITE_VCFHEADER_H
