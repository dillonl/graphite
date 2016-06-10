#ifndef GRAPHITE_SAMTOOLSALIGNMENTREADER_HPP
#define GRAPHITE_SAMTOOLSALIGNMENTREADER_HPP

#include "core/region/Region.h"
#include "IAlignment.h"
#include "IAlignmentList.h"
#include "IAlignmentReader.h"
#include "Sample.hpp"

#include <boost/noncopyable.hpp>

namespace graphite
{
	class SamtoolsAlignmentReader : public IAlignmentReader
	{
	public:
        typedef std::shared_ptr< SamtoolsAlignmentReader > SharedPtr;
        SamtoolsAlignmentReader(const std::string& path);
		~SamtoolsAlignmentReader();

        std::vector< IAlignment::SharedPtr > loadAlignmentsInRegion(Region::SharedPtr regionPtr, bool excludeDuplicateReads = false) override;
		void loadAlignmentSequencesInRegion(Region::SharedPtr regionPtr, std::shared_ptr< std::unordered_map< std::string, IAlignment::SharedPtr > > nameAlignmentPtrsMap);

		static std::vector< Sample::SharedPtr > GetBamReaderSamples(const std::string& bamPath);
		static position GetLastPositionInBam(const std::string& bamPath, Region::SharedPtr regionPtr);

	private:
		static void InitReader(const std::string& path);

		static std::unordered_map< std::string, position > s_region_last_positions;
		static std::unordered_map< std::string, uint32_t > s_region_map;
		static std::mutex s_region_ids_lock;
		static std::mutex s_region_last_positions_lock;
		std::string m_path;
	};
}

#endif //GRAPHITE_SAMTOOLSALIGNMENTREADER_HPP
