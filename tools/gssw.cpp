#include "core/alignments/BamAlignmentReader.h"
#include "core/alignments/BamAlignmentReaderManager.h"
#include "core/alignments/BamAlignmentReaderPreloadManager.h"
#include "core/variants/VCFFileReader.h"
#include "core/variants/IVariant.h"
#include "core/variants/VariantListVCFPreloaded.h"
#include "core/reference/FastaReference.h"
#include "core/utils/Parameters.h"
#include "gssw/graph/GraphManager.h"
#include "gssw/graph/GSSWAdjudicator.h"


int main(int argc, char** argv)
{
		gwiz::Parameters::Instance()->setParams(argc, argv);
		std::string fastaPath = gwiz::Parameters::Instance()->getParameterValue("fasta");
		std::string vcfPath = gwiz::Parameters::Instance()->getParameterValue("vcf");
		std::string bamPath = gwiz::Parameters::Instance()->getParameterValue("bam");

		gwiz::Region::SharedPtr regionPtr = std::make_shared< gwiz::Region >("20");
		auto fastaReferencePtr = std::make_shared< gwiz::FastaReference >(fastaPath, regionPtr);
		auto vcfFileReaderPtr = std::make_shared< gwiz::VariantListVCFPreloaded >(vcfPath, regionPtr);
		auto bamAlignmentReaderPreloadManager = std::make_shared< gwiz::BamAlignmentReaderPreloadManager >(bamPath, regionPtr);
		vcfFileReaderPtr->loadVariantsFromFile();

		std::cout << "Finished loading BAM and VCF" << std::endl;

		auto gsswAdjudicator = std::make_shared< gwiz::gssw::GSSWAdjudicator >();
		auto gsswGraphManager = std::make_shared< gwiz::gssw::GraphManager >(fastaReferencePtr, vcfFileReaderPtr, bamAlignmentReaderPreloadManager, gsswAdjudicator);
		auto variantList = gsswGraphManager->buildGraphs(fastaReferencePtr->getRegion(), 3000, 1000, 100);

		std::cout << "starting to print vcf" << std::endl;

		std::ofstream outVCF;
		outVCF.open("output.vcf", std::ios::out);
		variantList->printToVCF(outVCF);
		outVCF.close();
		return 0;
}
