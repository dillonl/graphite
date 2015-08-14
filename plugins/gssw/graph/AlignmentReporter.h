#ifndef GRAPHITE_GSSW_ALIGNMENTREPORTER_H
#define GRAPHITE_GSSW_ALIGNMENTREPORTER_H

#include "AlignmentReport.h"

#include <boost/noncopyable.hpp>
#include <memory>
#include <mutex>

namespace graphite
{
namespace gssw
{
	class AlignmentReporter : boost::noncopyable
	{
	public:
		static AlignmentReporter* Instance();

		void printAlignmentReportsToStream(std::ostream& out);
		void addAlignmentReport(AlignmentReport::SharedPtr alignmentReportPtr);

	private:
		AlignmentReporter() {}
		~AlignmentReporter() {}

		static AlignmentReporter* s_instance;

		std::mutex m_alignment_reports_mutex;
		std::list< AlignmentReport::SharedPtr > m_alignment_reports;
	};
}
}

#endif //GRAPHITE_GSSW_ALIGNMENTREPORTER_H
