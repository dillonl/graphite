#include "AlignmentReporter.h"


namespace graphite
{
namespace adjudicator
{
	AlignmentReporter* AlignmentReporter::s_instance = nullptr;


	AlignmentReporter* AlignmentReporter::Instance()
	{
		if (s_instance == nullptr)
		{
			s_instance = new AlignmentReporter();
		}
		return s_instance;
	}

	void AlignmentReporter::printAlignmentReportsToStream(std::ostream& out)
	{
		std::lock_guard< std::mutex > lock(this->m_alignment_reports_mutex); // unlocks when out of scope
		while (!this->m_alignment_reports.empty())
		{
			auto alignmentReportPtr = this->m_alignment_reports.back();
			this->m_alignment_reports.pop_back();
			auto alignmentReportString = alignmentReportPtr->toString();
			out.write(alignmentReportString.c_str(), alignmentReportString.size());
		}
	}

	void AlignmentReporter::addAlignmentReport(AlignmentReport::SharedPtr alignmentReportPtr)
	{
		std::lock_guard< std::mutex > lock(this->m_alignment_reports_mutex); // unlocks when out of scope
		this->m_alignment_reports.emplace_front(alignmentReportPtr);
	}

}
}
