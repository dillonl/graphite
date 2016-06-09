#ifndef GRAPHITE_S3FILETESTS_HPP
#define GRAPHITE_S3FILETESTS_HPP

#include "TestConfig.h"

#include "core/alignment/SamtoolsAlignmentReader.h"
#include "core/alignment/SampleManager.hpp"


#include <vector>

namespace
{
	namespace adj_test
	{
		using namespace graphite;

		TEST(S3FileTest, fileTest)
		{
			std::string path = "http://s3.amazonaws.com/iobio/NA12878/NA12878.autsome.bam";
			auto samplePtrs = graphite::SamtoolsAlignmentReader::GetBamReaderSamples(path);
			for (auto samplePtr : samplePtrs)
			{
				graphite::SampleManager::Instance()->addSamplePtr(samplePtr);
			}
			auto samtoolsReader = std::make_shared< SamtoolsAlignmentReader >(path);

		}
	}
}
#endif
