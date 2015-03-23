#include "AlignmentReporter.h"


namespace gwiz
{
namespace gssw
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

	void AlignmentReporter::AddAlignmentEvidence(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr, IAlignment::SharedPtr alignmentPtr, std::shared_ptr< gssw_graph_mapping > graphMappingSharePtr)
	{

	}

}
}
