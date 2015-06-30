#ifndef GWIZ_ALLELEMETADATA_H
#define GWIZ_ALLELEMETADATA_H

#include <boost/noncopyable.hpp>

namespace gwiz
{
	class VCFFileReader;
	class IAllele;
	class AlleleMetaData : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< AlleleMetaData > SharedPtr;
	    AlleleMetaData(std::shared_ptr< VCFFileReader > vcfFileReaderPtr, uint16_t paddingPrefix, uint16_t paddingSuffix) :
		    m_vcf_file_reader_ptr(vcfFileReaderPtr), m_padding_prefix(paddingPrefix), m_padding_suffix(paddingSuffix)
		{
		}
		~AlleleMetaData() {}

		std::shared_ptr< VCFFileReader > getVCFFileReaderWeakPtr() { return this->m_vcf_file_reader_ptr; }
		uint16_t getPaddingPrefix() { return this->m_padding_prefix; }
		uint16_t getPaddingSuffix() { return this->m_padding_suffix; }

	private:
		std::shared_ptr< VCFFileReader > m_vcf_file_reader_ptr;
		uint16_t m_padding_prefix;
		uint16_t m_padding_suffix;
	};
}

#endif //GWIZ_ALLELEMETADATA_H
