/*
#include "core/alignment/BamAlignmentReader.h"
#include "core/file/IFile.h"

#include "VisualizationToolKit.h"

namespace graphite
{
	VisualizationToolKit::VisualizationToolKit(const std::string& outputPath, const std::vector< std::string >& bamPaths, int matchValue) :
		m_output_path(outputPath), m_match_value(matchValue)
	{
		initializeKit(bamPaths);
	}

	VisualizationToolKit::~VisualizationToolKit()
	{

	}

	void VisualizationToolKit::initializeKit(const std::vector< std::string >& bamPaths)
	{
		for (auto bamPath : bamPaths)
		{
			auto samplePtrs = BamAlignmentReader::GetBamReaderSamples(bamPath);
			auto bamHeader = BamAlignmentReader::GetBamReaderHeader(bamPath);
			auto bamRefVec = BamAlignmentReader::GetBamReaderRefVector(bamPath);
			std::string fastaPath = m_output_path + "/" + IFile::getBaseNameWithoutExtension(bamPath) + ".graphite.fasta";
			std::string bamOutputPath = m_output_path + "/" + IFile::getBaseNameWithoutExtension(bamPath) + ".graphite.bam";
			auto fastaPtr = std::make_shared< FastaFileWriter >(fastaPath);
			auto bamWriterPtr = std::make_shared< BamFileWriter >(bamOutputPath, bamHeader, bamRefVec);
			m_bam_writer_ptrs.emplace_back(bamWriterPtr);
			m_fasta_writer_ptrs.emplace_back(fastaPtr);
			for (auto samplePtr : samplePtrs)
			{
				m_fasta_files.emplace(samplePtr->getReadgroup(), fastaPtr);
				m_bam_files.emplace(samplePtr->getReadgroup(), bamWriterPtr);
			}
		}
	}

	void VisualizationToolKit::setAlignmentAndMapping(IAlignment::SharedPtr alignmentPtr, GSSWGraph::SharedPtr gsswGraphPtr, GSSWMapping::SharedPtr refMapping, GSSWMapping::SharedPtr altMapping)
	{
		std::string refTracebackSequence = "";
		std::string refTracebackID = std::to_string(gsswGraphPtr->getStartPosition()) + "-" + std::to_string(gsswGraphPtr->getEndPosition());
		std::string altTracebackSequence = "";
		std::string altTracebackID = std::to_string(gsswGraphPtr->getStartPosition()) + "-" + std::to_string(gsswGraphPtr->getEndPosition());

		auto readGroup = alignmentPtr->getSample()->getReadgroup();
		auto fastaIter = m_fasta_files.find(readGroup);
		if (fastaIter != m_fasta_files.end())
		{
			FastaFileWriter::SharedPtr fastaWriterPtr = fastaIter->second;
			refMapping->setTracebackSequenceAndID(refTracebackSequence, refTracebackID);
			altMapping->setTracebackSequenceAndID(altTracebackSequence, altTracebackID);
			fastaWriterPtr->writeSequence(refTracebackID, refTracebackSequence);
			uint32_t altCounter = 0;
			for (auto allelePtr : altMapping->getAllelePtrs())
			{
				auto variantPtr = allelePtr->getVariantWPtr().lock();
				if (variantPtr != nullptr)
				{
					altCounter++;
					std::string tmpAltTracebackID = altTracebackID + "->" + std::to_string(variantPtr->getPosition()) + ":" + std::to_string(altCounter);
					fastaWriterPtr->writeSequence(tmpAltTracebackID, altTracebackSequence);
				}
			}
		}
		auto bamtoolsReaderIter = m_bam_files.find(readGroup);
		if (bamtoolsReaderIter != m_bam_files.end())
		{
			std::shared_ptr< BamTools::BamAlignment > bamtoolsAlignmentPtr = alignmentPtr->getBamToolsAlignmentPtr();
			auto bamtoolsWriterPtr = bamtoolsReaderIter->second;

			uint32_t swPercent = ((refMapping->getMappingScore() / (double)(alignmentPtr->getLength() * this->m_match_value)) * 100);
			std::string aName = bamtoolsAlignmentPtr->Name;
			bamtoolsAlignmentPtr->Name = aName + ":" + std::to_string(swPercent);
			bamtoolsAlignmentPtr->CigarData = cigarToBamToolsCigar(refMapping->getCigarData());
			bamtoolsAlignmentPtr->Position = refMapping->getAlignmentMappedPosition();
			bamtoolsWriterPtr->writeAlignment(bamtoolsAlignmentPtr, refTracebackID);

			uint32_t altCounter = 0;
			for (auto allelePtr : altMapping->getAllelePtrs())
			{
				auto variantPtr = allelePtr->getVariantWPtr().lock();
				if (variantPtr != nullptr)
				{
					altCounter++;
					std::string tmpAltTracebackID = altTracebackID + "->" +  std::to_string(variantPtr->getPosition()) + ":" + std::to_string(altCounter);

					swPercent = ((altMapping->getMappingScore() / (double)(alignmentPtr->getLength() * this->m_match_value)) * 100);
					bamtoolsAlignmentPtr->Name = aName + ":" + std::to_string(swPercent);
					bamtoolsAlignmentPtr->CigarData = cigarToBamToolsCigar(altMapping->getCigarData());
					bamtoolsAlignmentPtr->Position = altMapping->getAlignmentMappedPosition();
					bamtoolsWriterPtr->writeAlignment(bamtoolsAlignmentPtr, tmpAltTracebackID);
				}
			}
		}
    }

	std::vector< BamTools::CigarOp > VisualizationToolKit::cigarToBamToolsCigar(const std::vector< std::tuple< char, uint32_t > >& cigar)
	{
		std::vector< BamTools::CigarOp > cigarData;
		for (auto c : cigar)
		{
			BamTools::CigarOp op(std::get< 0 >(c), std::get< 1 >(c));
			cigarData.emplace_back(op);
		}
		return cigarData;
	}

	void VisualizationToolKit::closeResources()
	{
		for (auto fastaWriterPtr : m_fasta_writer_ptrs)
		{
			fastaWriterPtr->close();
		}
		for (auto bamWriterPtr : m_bam_writer_ptrs)
		{
			bamWriterPtr->close();
		}
	}
}
*/
