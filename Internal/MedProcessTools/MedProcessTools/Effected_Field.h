#ifndef __EFFECTED_FIELD__H__
#define __EFFECTED_FIELD__H__

#include <string>
using namespace std;
typedef enum {
	PREDICTION = 0,
	NUMERIC_ATTRIBUTE = 1,
	STRING_ATTRIBUTE = 2,
	JSON_DATA = 3,
	FEATURE = 4,

	LAST_UNDEFINED = 5
} Field_Type;

class Effected_Field {
public:
	Field_Type field = Field_Type::LAST_UNDEFINED;
	string value_name = "";
	Effected_Field() {}
	Effected_Field(Field_Type f, const string &val) {
		field = f;
		value_name = val;
	}

	bool operator==(const Effected_Field& other) const
	{
		if (this->field == other.field && this->value_name == other.value_name) return true;
		else return false;
	}

	struct HashFunction
	{
		size_t operator()(const Effected_Field& other) const
		{
			size_t xHash = std::hash<int>()(other.field);
			size_t yHash = std::hash<string>()(other.value_name) << 1;
			return xHash ^ yHash;
		}
	};
};

#endif
