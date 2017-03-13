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
		bool registerNewSample(std::shared_ptr< Sample > samplePtr) override;
		std::vector< std::string > getSampleNames() override;
		int32_t getColumnPosition(const std::string& columnTitle) override;
		void setColumn(const std::string& column) override;
		std::vector< std::string > getColumnNames() override;
		/* bool isActiveSampleColumnName(const std::string& headerName) override; */
		bool isNewSampleColumnName(const std::string& headerName) override;
		bool isSampleColumnName(const std::string& headerName) override;
	private:
		void setColumns(const std::string& headerString);
		std::string getColumnsString();
		void addFormatToHeader();

		std::string m_reference_path;
		std::string m_columns_line;
		std::vector< std::string > m_lines;
		std::vector< std::string > m_columns;
		std::vector< AlleleCountType > m_allele_count_type_registry;
		std::unordered_set< std::string > m_new_sample_names;
		std::unordered_set< std::string > m_sample_names;
		std::vector< std::string > m_sample_names_by_column_order;
	};
}

#endif //GRAPHITE_VCFHEADER_H
