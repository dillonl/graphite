#ifndef GRAPHITE_VCFHEADER_H
#define GRAPHITE_VCFHEADER_H

#include "IHeader.h"

#include <memory>
#include <vector>
#include <string>

namespace graphite
{
	class IAlignmentReader;
	class VCFHeader : public IHeader
	{
	public:
		typedef std::shared_ptr< VCFHeader > SharedPtr;
		VCFHeader();
		~VCFHeader();

		void addHeaderLine(const std::string& headerLine) override;
		std::string getHeader() override;
		void registerReferencePath(const std::string& referencePath);
		void registerSample(std::shared_ptr< Sample > samplePtr) override;
		std::vector< std::shared_ptr< Sample > > getSamplePtrs() override { return m_sample_ptrs; }
	private:
		std::string m_reference_path;
		std::vector< std::string > m_header_lines;
	    std::vector< std::shared_ptr< Sample > > m_sample_ptrs;
	};
}

#endif //GRAPHITE_VCFHEADER_H
