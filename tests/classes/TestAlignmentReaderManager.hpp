#ifndef GWIZ_TESTING_TESTALIGNMENTREADERMANAGER
#define GWIZ_TESTING_TESTALIGNMENTREADERMANAGER

namespace gwiz
{
namespace testing
{
	class TestAlignmentReaderManager : public IAlignmentReaderManager
	{
	public:
		typedef std::shared_ptr< TestAlignmentReaderManager > SharedPtr;
		TestAlignmentReaderManager() {}
		~TestAlignmentReaderManager() {}

		void addAlignments(IAlignmentReader::SharedPtr alignmentReader)
		{
			this->m_alignments = alignmentReader;
		}

		IAlignmentReader::SharedPtr generateAlignmentReader()
		{
			return this->m_alignments;
		}

	private:
		IAlignmentReader::SharedPtr m_alignments;

	};
}
}

#endif //GWIZ_TESTING_TESTALIGNMENTREADERMANAGER
