#ifndef GRAPHITE_VCFHEADER_H
#define GRAPHITE_VCFHEADER_H

#include "core/sample/SampleManager.h"
#include "core/util/Types.h"
#include "core/util/Noncopyable.hpp"

#include <memory>
#include <vector>
#include <string>
#include <unordered_set>

namespace graphite
{
	class IAlignmentReader;
	class VCFHeader : private Noncopyable
	{
	public:
		typedef std::shared_ptr< VCFHeader > SharedPtr;
		VCFHeader(const std::vector< std::string >& vcfHeaderLines);
		~VCFHeader();

		std::string getHeader();
		void registerReferencePath(const std::string& referencePath);
		void registerActiveSample(SampleManager::SharedPtr sampleManagerPtr);
		std::vector< std::string > getSampleNames();
		int32_t getColumnPosition(const std::string& columnTitle);
		void setColumn(const std::string& column);
		std::vector< std::string > getColumnNames();
		bool isActiveSampleColumnName(const std::string& headerName);
		bool isSampleColumnName(const std::string& headerName);
	private:
		void setColumns(const std::string& headerString);
		std::string getColumnsString();
		void addFormatToHeader();

		std::string m_reference_path;
		std::string m_columns_line;
		std::vector< std::string > m_lines;
		std::vector< std::string > m_columns;
		std::vector< AlleleCountType > m_allele_count_type_registry;
		std::unordered_set< std::string > m_active_sample_names;
		std::unordered_set< std::string > m_sample_names;
		std::vector< std::string > m_sample_names_by_column_order;
	};
}

#endif //GRAPHITE_VCFHEADER_H
