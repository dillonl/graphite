#ifndef GRAPHITE_IHEADER_H
#define GRAPHITE_IHEADER_H

#include <memory>
#include <vector>
#include <string>

#include <boost/noncopyable.hpp>

namespace graphite
{
	class Sample;
	class IHeader : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IHeader > SharedPtr;
		IHeader() {}
		virtual ~IHeader() {}

		virtual void addHeaderLine(const std::string& headerLine) = 0;
		virtual std::string getHeader() = 0;
		virtual void registerSample(std::shared_ptr< Sample > samplePtr) = 0;
		virtual std::vector< std::shared_ptr< Sample > > getSamplePtrs() = 0;
	private:
	};
}

#endif //GRAPHITE_IHEADER_H
