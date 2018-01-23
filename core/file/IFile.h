#ifndef GRAPHITE_IFILE_H
#define GRAPHITE_IFILE_H

#include <string>
#include <memory>
#include <iostream>
#include <fstream>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "core/util/Noncopyable.hpp"

/* #include <boost/filesystem.hpp> */

namespace graphite
{
	class IFile : private Noncopyable
	{
	public:
		typedef std::shared_ptr<IFile> SharedPtr;
	    IFile(const std::string& path)
			: m_file_path(path), m_opened(false)//, m_file_size(0)
		{
		}
		virtual ~IFile(){}

		virtual void Open() = 0;
		virtual void Close() = 0;
		/* virtual const char* getNextLine() = 0; */
		virtual bool getNextLine(std::string& line) = 0;
		virtual void setFilePosition(uint64_t pos) = 0;

		static bool fileExists(const std::string& name, bool exitOnFailure)
		{
			std::ifstream f(name.c_str());
			if (!f)
			{
				if (exitOnFailure)
				{
					std::cout << "File not found: " << name << std::endl;
					exit(EXIT_FAILURE);
				}
				return false;
			}
			return true;
		}

		static bool folderExists(const std::string& path, bool exitOnFailure)
		{
			if(!path.empty())
			{
				struct stat sb;

				if (stat(path.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
				{
					return true;
				}
			}
			if (exitOnFailure)
			{
                std::cout << "Folder not found: " << path << std::endl;
				exit(EXIT_FAILURE);
			}
			return false;
		}

		static std::string getBaseName(const std::string& path)
		{
			return path.substr(path.find_last_of("/") + 1);
		}

		static std::string getBaseNameWithoutExtension(const std::string& path)
		{
			return path.substr(path.find_last_of("/") + 1, path.find_last_of("."));
		}

		static std::string getExtension(const std::string& path)
		{
			return path.substr(path.find_last_of(".") + 1);
		}

	protected:
		inline bool fileExists(const std::string& filename)
		{
			/* return boost::filesystem::exists(name); */
			std::ifstream ifile(filename.c_str());
			return (bool)ifile;
		}

		std::string m_file_path;
		bool m_opened;
		/* size_t m_file_size; */
	};
} // end namespace graphite

#endif // GRAPHITE_IFILE_H
