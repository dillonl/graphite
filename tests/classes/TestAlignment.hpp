#ifndef GWIZ_TESTING_TESTALIGNMENT_HPP
#define GWIZ_TESTING_TESTALIGNMENT_HPP

namespace gwiz
{
namespace testing
{
	class TestAlignment : public IAlignment
	{
	public:
		typedef std::shared_ptr< TestAlignment > SharedPtr;
		TestAlignment(const position pos, const std::string& sequence) :
			m_position(pos), m_sequence(sequence)
		{
		}

		~TestAlignment() {}

		virtual const char* getSequence() override { return this->m_sequence.c_str(); };
		virtual const position getPosition() override { return this->m_position; }
		virtual const size_t getLength() override { return m_sequence.size(); }

	private:
		position m_position;
		std::string m_sequence;
	};
}
}

#endif //GWIZ_TESTING_TESTALIGNMENT_HPP
