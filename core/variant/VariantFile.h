#ifndef GWIZ_VARIANT_FILE_H
#define GWIZ_VARIANT_FILE_H

#include <boost/noncopyable.hpp>

namespace gwiz
{
	class VariantFile : private boost::noncopyable
	{
	public:
		VariantFile(const std::string& file_path);
		~VariantFile();

	private:
		VariantFile(const VariantFile& other) = delete; // no reason to use copy const
		VariantFile& operator=(const VariantFile& other) = delete; // no reason to use = operator

		std::string m_file_name;
	};
}

#endif //GWIZ_VARIANT_FILE_H
