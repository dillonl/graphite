#ifndef GWIZ_VCFALLELE_H
#define GWIZ_VCFALLELE_H

#include "Allele.h"

namespace gwiz
{
	class VCFAllele : public Allele
	{
	public:
		VCFAllele() {}
		~Allele() {}

		/*
		  We use a createAllelePtr static method so we can create a weak pointer
		  to the allelemetadata object. That way the allelemetadata has a way
		  to trace back to the allele and similarely we will have similar functionality
		  on the variant creation.
		*/
		static Allele::SharedPtr createAllelePtr(std::shared_ptr< Sequence > sequencePtr, std::shared_ptr< VCFFileReader > vcfFileReaderPtr, uint16_t paddingPrefix, uint16_t paddingSuffix)
		{
			auto allelePtr = std::shared_ptr< Allele >(new Allele(sequencePtr));
			allelePtr->m_allele_meta_data_ptrs.emplace_back(std::make_shared< AlleleMetaData >(vcfFileReaderPtr, paddingPrefix, paddingSuffix));
			return allelePtr;
		}

		std::vector< std::shared_ptr< AlleleMetaData > > getAlleleMetaDataPtrs() override { return m_allele_meta_data_ptrs; }

		void addAlleleMetaData(std::shared_ptr< VCFFileReader > vcfReaderPtr, uint16_t paddingPrefix, uint16_t paddingSuffix) override
		{
			this->m_allele_meta_data_ptrs.emplace_back(std::make_shared< AlleleMetaData >(vcfReaderPtr, paddingPrefix, paddingSuffix));
		}

	protected:
		std::vector< AlleleMetaData::SharedPtr > m_allele_meta_data_ptrs;

	};
}

#endif //GWIZ_VCFALLELE_H
