#ifndef GRAPHITE_ADJUDICATOR_ALIGNMENTREPORTER_H
#define GRAPHITE_ADJUDICATOR_ALIGNMENTREPORTER_H

#include "AlignmentReport.h"

#include <memory>
#include <mutex>
#include <list>

namespace graphite
{
namespace adjudicator
{
	class AlignmentReporter : Noncopyable
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
