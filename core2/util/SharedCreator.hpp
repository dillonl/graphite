#ifndef GRAPHITE_SHAREDCREATOR_H
#define GRAPHITE_SHAREDCREATOR_H

#include <memory>

namespace graphite
{
	template <class T>
	class SharedCreator
	{
		// private structs so the constructor and destructor can be private
		// forces use of the CreateVCFFileReader function
	private:
		struct _private_const {};
		struct _private_destr { void operator()(T* ptr) const { delete ptr; } };

	public:
		template < typename ...ArgsT >
		static std::shared_ptr< T > ConstructShared(...)
		{
			// auto sharedPtr = std::make_shared< T >({params});
			// return sharedPtr;
			return nullptr;
		}
	};
}

#endif //GRAPHITE_SHAREDCREATOR_H
