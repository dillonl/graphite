#ifndef GWIZ_GSSW_ALIGNMENTREPORTER_H
#define GWIZ_GSSW_ALIGNMENTREPORTER_H

#include <boost/noncopyable.hpp>
#include <memory>

namespace gwiz
{
namespace gssw
{
	class AlignmentReporter : boost::noncopyable
	{
	public:
		static AlignmentReporter* Instance();

		void AddAlignmentEvidence();

	private:
		AlignmentReporter() {}
		~AlignmentReporter() {}

		static AlignmentReporter* s_instance;
	};
}
}

#endif //GWIZ_GSSW_ALIGNMENTREPORTER_H
