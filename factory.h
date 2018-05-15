#include <map>
#include "base_classes.h"

namespace factory
{

	//************ FACTORY ***************
	template<typename T> Residue * createT()
	{
		return new T;
	}
	///
	typedef std::map<std::string, Residue*> map_type;
	///
	class ResidueFactory {
	public:
		inline static Residue* getInstance(const std::string & resname) {
			map_type::iterator it = ResidueMap->find(resname);
			if(it == ResidueMap->end())
				return NULL;
			return it->second;
		}

		static map_type * ResidueMap;
	};
	///
	template<typename T> class Registrar{ 
	public:
		Registrar(const std::string & resname)
		{ 
				ResidueFactory::ResidueMap->insert(std::pair<std::string,Residue*>(resname,createT<T>()));
		}
	};
	///
	//map_type * ResidueFactory::ResidueMap = new map_type;
	// ********** END OF FACTORY **********/

}