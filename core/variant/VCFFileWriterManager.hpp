#ifndef GRAPHITE_VCFFILEWRITERMANAGER_HPP
#define GRAPHITE_VCFFILEWRITERMANAGER_HPP

#include "VCFFileWriter.h"

#include <unordered_map>

namespace graphite
{
	class VCFFileWriterManager : private Noncopyable
	{
	public:
		static VCFFileWriterManager* Instance()
		{
			static VCFFileWriterManager* s_instance = new VCFFileWriterManager();
			return s_instance;
		}

		void addVCFFileWritersForVCFs(const std::vector< std::string >& originalVCFPaths, const std::string& outputDirectory)
		{
			for (auto originalVCFPath : originalVCFPaths)
			{
				addVCFFileWriterForVCF(originalVCFPath, outputDirectory, FileType::ASCII); // hard-coded for now
			}
		}

		VCFFileWriter::SharedPtr getVCFFileWriter(const std::string& originalVCFPath)
		{
			auto iter = m_vcf_writers.find(originalVCFPath);
			if (iter != m_vcf_writers.end())
			{
				return iter->second;
			}
			return nullptr;
		}

		void closeAllVCFFileWriters()
		{
			for (auto iter : m_vcf_writers)
			{
				iter.second->close();
			}
		}

	private:
		VCFFileWriterManager(){}
		~VCFFileWriterManager(){}

		void addVCFFileWriterForVCF(const std::string& originalVCFPath, const std::string& outputDirectory, FileType fileType)
		{
			std::string originalVcfFileName = originalVCFPath.substr(originalVCFPath.find_last_of("/") + 1);
			std::string outputVCFFilePath = outputDirectory + "/" + originalVcfFileName;
			uint32_t counter = 1;
			while (graphite::IFile::fileExists(outputVCFFilePath, false))
			{
				std::string extension = originalVCFPath.substr(originalVCFPath.find_last_of(".") + 1);
				std::string fileNameWithoutExtension = originalVcfFileName.substr(0, originalVcfFileName.find_last_of("."));
				outputVCFFilePath = outputDirectory + "/" + fileNameWithoutExtension + "." + std::to_string(counter) + "." + extension;
				++counter;
			}
			auto vcfFileWriterPtr = std::make_shared< VCFFileWriter >(outputVCFFilePath, originalVCFPath, fileType);
			m_vcf_writers.emplace(originalVCFPath, vcfFileWriterPtr);
		}

		std::unordered_map< std::string, VCFFileWriter::SharedPtr > m_vcf_writers;
	};
}

#endif //GRAPHITE_VCFFILEWRITERMANAGER_HPP
